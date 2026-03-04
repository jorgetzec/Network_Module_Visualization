from __future__ import annotations

import io
from typing import Any

import pandas as pd
import streamlit as st

from .core import find_connected_genes, load_network_from_dataframes
from .plotting import draw_highlighted_network, draw_zoomed_network
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

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Nodes", graph.number_of_nodes())
    col2.metric("Edges", graph.number_of_edges())
    col3.metric("Found targets", len(found_targets))
    col4.metric("Connected genes", len(connected_genes))
    st.caption(
        "These metrics summarize the filtered network after applying weight threshold and min degree."
    )

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

    st.subheader("Network visualizers")
    st.markdown(
        "Use these two panels together: the left panel shows global context, and the right panel focuses on the local neighborhood of target genes."
    )
    col_full, col_zoom = st.columns(2)
    with col_full:
        st.markdown("**Full network**")
        st.caption(
            "Complete filtered network. Target genes and their direct links are highlighted, while the rest of the graph remains in the background."
        )
        st.pyplot(fig, use_container_width=True)
    with col_zoom:
        st.markdown("**Target subnetwork (zoom)**")
        if zoom_fig is not None:
            st.caption(
                "Zoomed 1-hop subnetwork: only target genes and their direct neighbors. This is useful to inspect local structure."
            )
            st.pyplot(zoom_fig, use_container_width=True)
        else:
            st.info("Zoom view unavailable. Ensure at least one target gene exists in the network.")

    if connected_genes:
        rows = []
        for gene in connected_genes:
            connected_to = [target for target in found_targets if target in graph and gene in graph.neighbors(target)]
            rows.append(
                {
                    "gene": gene,
                    "degree": graph.degree(gene),
                    "connected_to_targets": ",".join(connected_to),
                }
            )
        connected_df = pd.DataFrame(rows).sort_values("degree", ascending=False)
        st.subheader("Connected genes table (target subnetwork)")
        st.caption(
            "Rows are genes directly connected to at least one target gene (1-hop neighborhood). "
            "`degree` is computed on the full filtered network."
        )
        st.dataframe(connected_df, use_container_width=True)
        st.download_button(
            "Download connected_genes.csv",
            connected_df.to_csv(index=False).encode("utf-8"),
            file_name="connected_genes.csv",
            mime="text/csv",
        )

    stats_df = generate_network_statistics_df(graph)
    st.subheader("Connectivity summary table (full network)")
    st.caption(
        "This table includes all genes in the filtered network, sorted by degree. "
        "`connections` lists each gene's neighbors in the full filtered graph."
    )
    st.dataframe(stats_df, use_container_width=True)
    st.download_button(
        "Download network_statistics.csv",
        stats_df.to_csv(index=False).encode("utf-8"),
        file_name="network_statistics.csv",
        mime="text/csv",
    )

    pdf_buffer = io.BytesIO()
    fig.savefig(pdf_buffer, format="pdf", bbox_inches="tight", pad_inches=0.1)
    pdf_buffer.seek(0)
    st.download_button(
        "Download figure PDF",
        data=pdf_buffer,
        file_name="network_highlighted.pdf",
        mime="application/pdf",
    )

    if zoom_fig is not None:
        zoom_pdf_buffer = io.BytesIO()
        zoom_fig.savefig(zoom_pdf_buffer, format="pdf", bbox_inches="tight", pad_inches=0.1)
        zoom_pdf_buffer.seek(0)
        st.download_button(
            "Download zoom figure PDF",
            data=zoom_pdf_buffer,
            file_name="network_zoomed.pdf",
            mime="application/pdf",
        )


def run_app() -> None:
    st.set_page_config(page_title="Network Module Visualization", layout="wide")
    st.title("Network Module Visualization")
    st.caption("Visualize WGCNA/Cytoscape/VisANT networks and highlight target genes.")

    if "analysis_payload" not in st.session_state:
        st.session_state["analysis_payload"] = None

    with st.sidebar:
        st.header("Configuration")
        with st.form("network_form"):
            nodes_file = st.file_uploader("Nodes file (TSV)", type=["txt", "tsv"])
            edges_file = st.file_uploader("Edges file (TSV)", type=["txt", "tsv"])
            weight = st.slider("Weight threshold", min_value=0.0, max_value=1.0, value=0.2, step=0.01)
            min_degree = st.number_input("Min degree", min_value=1, max_value=1000, value=1, step=1)
            module_color = st.text_input("Module color", value="red")
            show_only_target_labels = st.checkbox("Show only target labels", value=False)
            genes_text = st.text_area("Target genes (comma or newline separated)", value="")
            run_btn = st.form_submit_button("Run")

    if run_btn:
        if not nodes_file or not edges_file:
            st.error("You must upload nodes and edges files.")
        else:
            try:
                nodes_df = pd.read_csv(nodes_file, sep="\t")
                edges_df = pd.read_csv(edges_file, sep="\t")
                genes = _parse_genes(genes_text)

                graph = load_network_from_dataframes(
                    nodes_df,
                    edges_df,
                    weight_threshold=weight,
                    min_degree=int(min_degree),
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
        st.info("Upload files and click Run.")
