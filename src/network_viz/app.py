from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st

from .core import find_connected_genes, load_network_from_dataframes
from .plotting import (
    draw_highlighted_network,
    draw_interactive_network_plotly,
    draw_interactive_network_plotly_3d,
    draw_zoomed_network,
)
from .reports import generate_network_statistics_df


def _parse_genes(raw_text: str) -> list[str]:
    genes = [gene.strip() for gene in raw_text.replace("\n", ",").split(",") if gene.strip()]
    seen: set[str] = set()
    unique_genes: list[str] = []
    for gene in genes:
        if gene not in seen:
            unique_genes.append(gene)
            seen.add(gene)
    return unique_genes


def _render_results(payload: dict[str, Any]) -> None:
    graph = payload["graph"]
    found_targets = payload["found_targets"]
    connected_genes = payload["connected_genes"]
    module_color = payload["module_color"]
    show_only_target_labels = payload["show_only_target_labels"]

    # --- Header Metrics ---
    st.markdown("### :material/analytics: Network Overview")
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric("Total Nodes", graph.number_of_nodes())
    with m2:
        st.metric("Total Edges", graph.number_of_edges())
    with m3:
        st.metric("Targets Found", len(found_targets))
    with m4:
        st.metric("Connected Genes", len(connected_genes))

    st.divider()

    # --- Organize with Tabs ---
    tab_viz, tab_interactive, tab_data, tab_stats, tab_export = st.tabs([
        ":material/hub: Visualizers",
        ":material/touch_app: Interactive Map",
        ":material/account_tree: Connected Genes",
        ":material/insights: Global Stats",
        ":material/file_download: Exports"
    ])

    # Visualizers Tab
    with tab_viz:
        col_full, col_zoom = st.columns(2)

        # Generate figures
        with st.spinner("Generating network diagrams..."):
            fig, pos_full = draw_highlighted_network(
                graph,
                found_targets,
                connected_genes,
                module_color=module_color,
                show_only_target_labels=show_only_target_labels,
            )
            zoom_fig = draw_zoomed_network(
                graph,
                found_targets,
                connected_genes,
                module_color=module_color,
                pos_full=pos_full,
            )

        with col_full:
            st.markdown("#### Full Filtered Network")
            st.caption("Global topology with target genes highlighted.")
            st.pyplot(fig, width='stretch')

        with col_zoom:
            st.markdown("#### Neighborhood Zoom")
            if zoom_fig is not None:
                st.caption("Direct neighbors (1-hop) of discovered targets.")
                st.pyplot(zoom_fig, width='stretch')
            else:
                st.warning("Zoom unavailable: No target genes were found.")
 
    # Interactive Tab
    with tab_interactive:
        st.markdown("#### Interactive Exploration")
        col_ctrl, _ = st.columns([1, 3])
        with col_ctrl:
            viz_mode = st.radio("Visualization Mode", ["2D Map", "3D Sphere"], horizontal=True, label_visibility="collapsed")
        
        st.caption(f"Showcasing {viz_mode}. Pan, zoom, and hover for metadata.")
        
        with st.spinner(f"Preparing {viz_mode}..."):
            if viz_mode == "2D Map":
                plotly_fig = draw_interactive_network_plotly(
                    graph,
                    found_targets,
                    connected_genes,
                    module_color=module_color
                )
            else:
                plotly_fig = draw_interactive_network_plotly_3d(
                    graph,
                    found_targets,
                    connected_genes,
                    module_color=module_color
                )
            st.plotly_chart(plotly_fig, width='stretch')

    # Connected Genes Tab
    with tab_data:
        if connected_genes:
            rows = []
            for gene in connected_genes:
                connected_to = [target for target in found_targets if target in graph and gene in graph.neighbors(target)]
                rows.append(
                    {
                        "gene": gene,
                        "degree": graph.degree(gene),
                        "connected_to_targets": ", ".join(connected_to),
                    }
                )
            connected_df = pd.DataFrame(rows).sort_values("degree", ascending=False)
            st.markdown("#### Target Neighborhood Analysis")
            st.dataframe(connected_df, width='stretch', hide_index=True)

            st.download_button(
                "Download CSV",
                connected_df.to_csv(index=False).encode("utf-8"),
                file_name="connected_genes.csv",
                mime="text/csv",
                key="dl_connected",
                icon=":material/download:"
            )
        else:
            st.info("No connections found.")

    # Global Stats Tab
    with tab_stats:
        stats_df = generate_network_statistics_df(graph)
        st.markdown("#### Full Network Statistics")
        st.dataframe(stats_df, width='stretch', hide_index=True)

        st.download_button(
            "Download CSV",
            stats_df.to_csv(index=False).encode("utf-8"),
            file_name="network_statistics.csv",
            mime="text/csv",
            key="dl_stats",
            icon=":material/download:"
        )

    # Exports Tab (PDFs)
    with tab_export:
        st.markdown("#### Download Publication-Ready Figures")
        st.info("Static plots with transparent background halos (High-quality PDF).")

        c1, c2 = st.columns(2)
        with c1:
            pdf_buffer = io.BytesIO()
            fig.savefig(pdf_buffer, format="pdf", bbox_inches="tight", pad_inches=0.1)
            st.download_button(
                "Full Network (PDF)",
                data=pdf_buffer.getvalue(),
                file_name="network_full.pdf",
                mime="application/pdf",
                width='stretch',
                icon=":material/picture_as_pdf:"
            )

        with c2:
            if zoom_fig is not None:
                zoom_pdf_buffer = io.BytesIO()
                zoom_fig.savefig(zoom_pdf_buffer, format="pdf", bbox_inches="tight", pad_inches=0.1)
                st.download_button(
                    "Zoomed View (PDF)",
                    data=zoom_pdf_buffer.getvalue(),
                    file_name="network_zoom.pdf",
                    mime="application/pdf",
                    width='stretch',
                    icon=":material/picture_as_pdf:"
                )


def run_app() -> None:
    st.set_page_config(
        page_title="Network Viz Pro",
        page_icon=":material/hub:",
        layout="wide"
    )

    st.markdown("""
        <style>
        .main { background-color: #f8f9fa; }
        .stMetric { background-color: #ffffff; padding: 15px; border-radius: 10px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); }
        </style>
    """, unsafe_allow_html=True)

    st.title(":material/hub: Network Module Visualization")
    st.markdown("---")

    if "analysis_payload" not in st.session_state:
        st.session_state["analysis_payload"] = None

    with st.sidebar:
        st.header(":material/settings: Configuration")

        with st.expander(":material/folder: Data Files", expanded=True):
            nodes_file = st.file_uploader("Nodes (TSV)", type=["txt", "tsv"])
            edges_file = st.file_uploader("Edges (TSV)", type=["txt", "tsv"])

        # Auto-detect columns and colors if files are uploaded
        node_id_col = "nodeName"
        node_label_col = "altName"
        from_node_col = "fromNode"
        to_node_col = "toNode"
        weight_col = "weight"
        suggested_color = "red"

        if nodes_file and edges_file:
            try:
                # Preview headers
                n_df_head = pd.read_csv(nodes_file, sep="\t", nrows=1)
                e_df_head = pd.read_csv(edges_file, sep="\t", nrows=1)
                
                with st.expander(":material/list_alt: Mapping", expanded=True):
                    st.write("**Nodes Mapping**")
                    cols_n = list(n_df_head.columns)
                    node_id_col = st.selectbox("Unique ID Column", cols_n, index=cols_n.index("nodeName") if "nodeName" in cols_n else 0)
                    node_label_col = st.selectbox("Display Name Column", cols_n, index=cols_n.index("altName") if "altName" in cols_n else 0)
                    
                    st.write("**Edges Mapping**")
                    cols_e = list(e_df_head.columns)
                    from_node_col = st.selectbox("From Node Column", cols_e, index=cols_e.index("fromNode") if "fromNode" in cols_e else 0)
                    to_node_col = st.selectbox("To Node Column", cols_e, index=cols_e.index("toNode") if "toNode" in cols_e else 0)
                    weight_col = st.selectbox("Weight Column", cols_e, index=cols_e.index("weight") if "weight" in cols_e else 0)

                # Attempt to guess color
                # If there is a column with "Attr", check for values
                attr_cols = [c for c in cols_n if "Attr" in c]
                if attr_cols:
                    # Read a few more rows to see a common value
                    n_df_peek = pd.read_csv(nodes_file, sep="\t", nrows=10)
                    val = n_df_peek[attr_cols[0]].iloc[0]
                    if isinstance(val, str) and len(val) > 2:
                        suggested_color = val
            except:
                pass

        with st.expander(":material/tune: Filters", expanded=True):
            weight = st.slider("Weight threshold", 0.0, 1.0, 0.2, 0.01)
            min_degree = st.number_input("Min degree", 1, 1000, 1)

        with st.expander(":material/palette: Appearance", expanded=False):
            module_color = st.text_input("Module color", value=suggested_color)
            show_only_target_labels = st.checkbox("Show only target labels", value=False)

        with st.expander(":material/gps_fixed: Targets", expanded=True):
            genes_text = st.text_area("Gene names (comma/newline)", placeholder="e.g. HLM1, THO2")

        run_btn = st.button("Run Analysis", width='stretch', type="primary", icon=":material/rocket_launch:")

    if run_btn:
        if not nodes_file or not edges_file:
            st.error("Please upload both Nodes and Edges files.")
        else:
            with st.spinner("Processing network data..."):
                try:
                    # Reset buffers for reading
                    nodes_file.seek(0)
                    edges_file.seek(0)
                    
                    nodes_df = pd.read_csv(nodes_file, sep="\t")
                    edges_df = pd.read_csv(edges_file, sep="\t")
                    genes = _parse_genes(genes_text)

                    graph = load_network_from_dataframes(
                        nodes_df,
                        edges_df,
                        weight_threshold=weight,
                        min_degree=int(min_degree),
                        node_id_col=node_id_col,
                        node_label_col=node_label_col,
                        from_node_col=from_node_col,
                        to_node_col=to_node_col,
                        weight_col=weight_col
                    )
                    found_targets, connected_genes = find_connected_genes(graph, genes)
                    st.session_state["analysis_payload"] = {
                        "graph": graph,
                        "found_targets": found_targets,
                        "connected_genes": connected_genes,
                        "module_color": module_color,
                        "show_only_target_labels": show_only_target_labels,
                    }
                except Exception as exc:
                    st.session_state["analysis_payload"] = None
                    st.exception(exc)

    payload = st.session_state.get("analysis_payload")
    if payload:
        _render_results(payload)
    else:
        st.info(":material/info: Welcome! Use the sidebar to upload files and define targets.")
