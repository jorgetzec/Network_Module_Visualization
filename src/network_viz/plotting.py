from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import networkx as nx
import plotly.graph_objects as go


def build_suffix(weight_threshold: float, module_color: str) -> str:
    current_date = datetime.now().strftime("%Y%m%d")
    return f"_w{weight_threshold}_{module_color}_{current_date}"


def draw_highlighted_network(
    graph: nx.Graph,
    found_targets: list[str],
    connected_genes: list[str],
    module_color: str,
    show_only_target_labels: bool = False,
) -> tuple[plt.Figure, dict[str, tuple[float, float]]]:
    fig = plt.figure(figsize=(15, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    try:
        pos = nx.kamada_kawai_layout(graph)
    except ModuleNotFoundError:
        # Fallback when scipy is not installed in the runtime environment.
        pos = nx.spring_layout(graph, seed=42)

    other_nodes = [node for node in graph.nodes() if node not in found_targets and node not in connected_genes]
    connected_nodes = [node for node in connected_genes if node in graph]

    other_edges = [(u, v) for u, v in graph.edges() if not (u in found_targets or v in found_targets)]
    other_edge_widths = [10 * graph[u][v]["weight"] for u, v in other_edges]
    if other_edges:
        nx.draw_networkx_edges(graph, pos, edgelist=other_edges, width=other_edge_widths, alpha=0.6, edge_color="silver")

    if other_nodes:
        other_node_sizes = [50 + 5 * graph.degree(node) for node in other_nodes]
        nx.draw_networkx_nodes(graph, pos, nodelist=other_nodes, node_size=other_node_sizes, node_color="gray", alpha=0.8)

    target_edges = []
    target_edge_widths = []
    for target in found_targets:
        if target not in graph:
            continue
        for neighbor in graph.neighbors(target):
            edge_tuple = (target, neighbor)
            if edge_tuple in target_edges:
                continue
            target_edges.append(edge_tuple)
            target_edge_widths.append(max(3.0, 20 * graph[target][neighbor]["weight"]))

    if target_edges:
        edge_color = "black" if module_color == "black" else module_color
        ax = plt.gca()
        for i, (u, v) in enumerate(target_edges):
            ax.plot(
                [pos[u][0], pos[v][0]],
                [pos[u][1], pos[v][1]],
                color=edge_color,
                linewidth=target_edge_widths[i],
                alpha=1.0,
                solid_capstyle="round",
                zorder=3,
            )

    if connected_nodes:
        connected_node_sizes = [50 + 5 * graph.degree(node) for node in connected_nodes]
        ax = plt.gca()
        ax.scatter(
            [pos[node][0] for node in connected_nodes],
            [pos[node][1] for node in connected_nodes],
            s=connected_node_sizes,
            c="white",
            edgecolors=module_color,
            linewidths=2,
            alpha=0.9,
            zorder=4,
        )

    if found_targets:
        target_node_sizes = [50 + 5 * graph.degree(node) for node in found_targets if node in graph]
        valid_targets = [node for node in found_targets if node in graph]
        ax = plt.gca()
        ax.scatter(
            [pos[node][0] for node in valid_targets],
            [pos[node][1] for node in valid_targets],
            s=target_node_sizes,
            c=module_color,
            alpha=1.0,
            zorder=5,
        )

    if show_only_target_labels:
        labels_to_show = set(found_targets) | set(connected_nodes)
        label_dict = {node: node for node in labels_to_show if node in graph}
    else:
        label_dict = {node: node for node in graph.nodes()}

    # Improved Label Drawing: added halo (path effect) for readability without opaque boxes
    labels = nx.draw_networkx_labels(
        graph,
        pos,
        labels=label_dict,
        font_size=6,
        font_color="black",
    )

    # Use a white halo (stroke) for better visibility against backgrounds
    halo = [path_effects.withStroke(linewidth=3, foreground="white", alpha=0.8)]

    for node, t in labels.items():
        t.set_zorder(10)
        t.set_path_effects(halo)

        if node in found_targets:
            t.set_fontsize(12)
            t.set_color(module_color if module_color != "white" else "black")
            t.set_weight("bold")
            # Stronger halo for targets
            t.set_path_effects([path_effects.withStroke(linewidth=4, foreground="white", alpha=1.0)])
        elif node in connected_nodes:
            t.set_fontsize(8)

    target_str = ", ".join(found_targets) if found_targets else "No target genes"
    plt.title(f"Gene network highlighting {target_str} ({module_color} module)", fontsize=16, fontweight="bold", pad=20)
    plt.axis("off")
    return fig, pos


def draw_zoomed_network(
    graph: nx.Graph,
    found_targets: list[str],
    connected_genes: list[str],
    module_color: str,
    pos_full: dict[str, tuple[float, float]] | None = None,
) -> plt.Figure | None:
    if not found_targets:
        return None

    nodes_to_include = set(found_targets)
    for target in found_targets:
        if target in graph:
            nodes_to_include.update(graph.neighbors(target))

    nodes_to_include = [node for node in nodes_to_include if node in graph]
    if len(nodes_to_include) < 2:
        return None

    graph_zoom = graph.subgraph(nodes_to_include).copy()

    if pos_full is None:
        try:
            pos_full = nx.kamada_kawai_layout(graph)
        except ModuleNotFoundError:
            pos_full = nx.spring_layout(graph, seed=42)

    pos_zoom = {node: pos_full[node] for node in nodes_to_include if node in pos_full}
    if not pos_zoom:
        return None

    x_coords = [pos[0] for pos in pos_zoom.values()]
    y_coords = [pos[1] for pos in pos_zoom.values()]
    center_x = sum(x_coords) / len(x_coords)
    center_y = sum(y_coords) / len(y_coords)
    scale_factor = 3.0

    pos_scaled: dict[str, tuple[float, float]] = {}
    for node, (x, y) in pos_zoom.items():
        pos_scaled[node] = ((x - center_x) * scale_factor + center_x, (y - center_y) * scale_factor + center_y)

    connected_nodes = [node for node in connected_genes if node in graph_zoom]
    other_nodes = [node for node in graph_zoom.nodes() if node not in found_targets and node not in connected_nodes]

    fig = plt.figure(figsize=(15, 10))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

    other_edges = [(u, v) for u, v in graph_zoom.edges() if not (u in found_targets or v in found_targets)]
    other_edge_widths = [10 * graph_zoom[u][v].get("weight", 0.1) for u, v in other_edges]
    if other_edges:
        nx.draw_networkx_edges(
            graph_zoom,
            pos_scaled,
            edgelist=other_edges,
            width=other_edge_widths,
            alpha=0.3,
            edge_color="silver",
        )

    if other_nodes:
        other_node_sizes = [(50 + 5 * graph.degree(node)) * 1.3 for node in other_nodes]
        nx.draw_networkx_nodes(
            graph_zoom,
            pos_scaled,
            nodelist=other_nodes,
            node_size=other_node_sizes,
            node_color="gray",
            alpha=0.8,
        )

    target_edges = []
    target_edge_widths = []
    for target in found_targets:
        if target not in graph_zoom:
            continue
        for neighbor in graph_zoom.neighbors(target):
            edge_tuple = (target, neighbor)
            if edge_tuple in target_edges:
                continue
            target_edges.append(edge_tuple)
            target_edge_widths.append(max(4.0, 30 * graph_zoom[target][neighbor].get("weight", 0.1)))

    if target_edges:
        edge_color = "black" if module_color == "black" else module_color
        ax = plt.gca()
        for i, (u, v) in enumerate(target_edges):
            ax.plot(
                [pos_scaled[u][0], pos_scaled[v][0]],
                [pos_scaled[u][1], pos_scaled[v][1]],
                color=edge_color,
                linewidth=target_edge_widths[i],
                alpha=1.0,
                solid_capstyle="round",
                zorder=3,
            )

    if connected_nodes:
        connected_node_sizes = [(50 + 5 * graph.degree(node)) * 1.3 for node in connected_nodes]
        ax = plt.gca()
        ax.scatter(
            [pos_scaled[node][0] for node in connected_nodes],
            [pos_scaled[node][1] for node in connected_nodes],
            s=connected_node_sizes,
            c="white",
            edgecolors=module_color,
            linewidths=2,
            alpha=0.9,
            zorder=4,
        )

    valid_targets = [target for target in found_targets if target in graph_zoom]
    if valid_targets:
        target_node_sizes = [(50 + 5 * graph.degree(node)) * 1.3 for node in valid_targets]
        ax = plt.gca()
        ax.scatter(
            [pos_scaled[node][0] for node in valid_targets],
            [pos_scaled[node][1] for node in valid_targets],
            s=target_node_sizes,
            c=module_color,
            alpha=1.0,
            zorder=5,
        )

    label_dict = {node: node for node in graph_zoom.nodes()}
    labels = nx.draw_networkx_labels(
        graph_zoom,
        pos_scaled,
        labels=label_dict,
        font_size=10,
        font_color="black"
    )

    halo = [path_effects.withStroke(linewidth=4, foreground="white", alpha=0.9)]

    for node, t in labels.items():
        t.set_zorder(10)
        t.set_path_effects(halo)

        if node in valid_targets:
            t.set_fontsize(16)
            t.set_color(module_color if module_color != "white" else "black")
            t.set_weight("bold")
            t.set_path_effects([path_effects.withStroke(linewidth=5, foreground="white", alpha=1.0)])
        elif node in connected_nodes:
            t.set_fontsize(12)

    target_str = ", ".join(valid_targets)
    plt.title(f"Zoomed view: {target_str} and neighbors ({module_color} module)", fontsize=20, fontweight="bold", pad=25)
    plt.axis("off")
    return fig


def save_network_figure(
    fig: plt.Figure, output_prefix: str, weight_threshold: float, module_color: str, highlighted: bool = True
) -> Path:
    suffix = build_suffix(weight_threshold, module_color)
    filename = f"{output_prefix}_network_highlighted{suffix}.pdf" if highlighted else f"{output_prefix}_network{suffix}.pdf"
    out_path = Path(filename)
    fig.savefig(out_path, bbox_inches="tight", pad_inches=0.1)
    return out_path


def draw_interactive_network_plotly(
    graph: nx.Graph,
    found_targets: list[str],
    connected_genes: list[str],
    module_color: str,
) -> go.Figure:
    """Creates an enhanced interactive Plotly 2D visualization for Streamlit."""
    try:
        pos = nx.kamada_kawai_layout(graph)
    except:
        pos = nx.spring_layout(graph, seed=42)

    # Separate edges for highlighting
    target_edge_x, target_edge_y = [], []
    normal_edge_x, normal_edge_y = [], []

    for edge in graph.edges():
        u, v = edge
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        
        if u in found_targets or v in found_targets:
            target_edge_x.extend([x0, x1, None])
            target_edge_y.extend([y0, y1, None])
        else:
            normal_edge_x.extend([x0, x1, None])
            normal_edge_y.extend([y0, y1, None])

    # Color safety
    common_css_colors = {'red', 'blue', 'green', 'black', 'white', 'gray', 'yellow', 'cyan', 'magenta', 'orange', 'purple', 'salmon', 'turquoise', 'brown', 'gold', 'silver', 'darksalmon', 'lightsalmon'}
    def is_safe(c: str) -> bool:
        return isinstance(c, str) and (c.startswith('#') or c.lower() in common_css_colors)
    
    safe_module_color = module_color if is_safe(module_color) else 'black'

    normal_edge_trace = go.Scatter(
        x=normal_edge_x, y=normal_edge_y,
        line=dict(width=0.5, color='rgba(200,200,200,0.5)'),
        hoverinfo='none',
        mode='lines',
        name='Normal'
    )

    target_edge_trace = go.Scatter(
        x=target_edge_x, y=target_edge_y,
        line=dict(width=2, color=safe_module_color),
        hoverinfo='none',
        mode='lines',
        name='High-Weight'
    )

    node_x, node_y, node_text, node_color, node_size, node_border = [], [], [], [], [], []

    for node in graph.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
        
        attrs = graph.nodes[node]
        label = attrs.get("display_label", node)
        degree = graph.degree(node)
        
        hover = f"<b>ID:</b> {node}<br><b>Label:</b> {label}<br><b>Degree:</b> {degree}"
        for k, v in attrs.items():
            if k not in ["display_label", "pos"]:
                hover += f"<br><b>{k}:</b> {v}"
        node_text.append(hover)
        
        if node in found_targets:
            node_color.append(safe_module_color)
            node_size.append(30)
            node_border.append(3)
        elif node in connected_genes:
            node_color.append("white")
            node_size.append(20)
            node_border.append(2)
        else:
            node_color.append("#D3D3D3")
            node_size.append(12)
            node_border.append(1)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text' if len(graph.nodes()) < 40 else 'markers',
        hoverinfo='text',
        text=[graph.nodes[n].get("display_label", n) for n in graph.nodes()],
        textposition="top center",
        marker=dict(
            showscale=False,
            color=node_color,
            size=node_size,
            line=dict(color='rgba(50,50,50,0.8)', width=node_border)
        )
    )
    node_trace.hovertext = node_text

    fig = go.Figure(data=[normal_edge_trace, target_edge_trace, node_trace],
                 layout=go.Layout(
                    title={'text': f'Interactive Network ({len(found_targets)} targets)', 'font': {'size': 16}},
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20, l=5, r=5, t=40),
                    height=800,  # Make it more square-like
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)'
                ))
    
    return fig


def draw_interactive_network_plotly_3d(
    graph: nx.Graph,
    found_targets: list[str],
    connected_genes: list[str],
    module_color: str,
) -> go.Figure:
    """Creates an enhanced interactive Plotly 3D visualization for Streamlit."""
    # Compute 3D layout
    try:
        pos = nx.kamada_kawai_layout(graph, dim=3)
    except:
        pos = nx.spring_layout(graph, dim=3, seed=42)

    # Separate edges for highlighting
    target_edge_x, target_edge_y, target_edge_z = [], [], []
    normal_edge_x, normal_edge_y, normal_edge_z = [], [], []

    for edge in graph.edges():
        u, v = edge
        x0, y0, z0 = pos[u]
        x1, y1, z1 = pos[v]
        
        if u in found_targets or v in found_targets:
            target_edge_x.extend([x0, x1, None])
            target_edge_y.extend([y0, y1, None])
            target_edge_z.extend([z0, z1, None])
        else:
            normal_edge_x.extend([x0, x1, None])
            normal_edge_y.extend([y0, y1, None])
            normal_edge_z.extend([z0, z1, None])

    # Color safety for Plotly
    common_css_colors = {'red', 'blue', 'green', 'black', 'white', 'gray', 'yellow', 'cyan', 'magenta', 'orange', 'purple', 'salmon', 'turquoise', 'brown', 'gold', 'silver'}
    def is_safe(c: str) -> bool:
        return isinstance(c, str) and (c.startswith('#') or c.lower() in common_css_colors)
    
    safe_module_color = module_color if is_safe(module_color) else 'black'

    # Normal edges: light gray for background
    normal_edge_trace = go.Scatter3d(
        x=normal_edge_x, y=normal_edge_y, z=normal_edge_z,
        line=dict(width=2, color='rgba(200,200,200,0.5)'),
        hoverinfo='none',
        mode='lines',
        name='Other Connections'
    )

    # Highlighted edges: very thick and module-colored
    target_edge_trace = go.Scatter3d(
        x=target_edge_x, y=target_edge_y, z=target_edge_z,
        line=dict(width=10, color=safe_module_color), # Even thicker
        hoverinfo='none',
        mode='lines',
        name='Target Connections'
    )

    node_x, node_y, node_z, node_text, node_color, node_size = [], [], [], [], [], []

    for node in graph.nodes():
        x, y, z = pos[node]
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)
        
        attrs = graph.nodes[node]
        label = attrs.get("display_label", node)
        hover = f"<b>ID:</b> {node}<br><b>Label:</b> {label}<br><b>Degree:</b> {graph.degree(node)}"
        for k, v in attrs.items():
            if k not in ["display_label", "pos"]:
                hover += f"<br><b>{k}:</b> {v}"
        node_text.append(hover)
        
        if node in found_targets:
            node_color.append(safe_module_color)
            node_size.append(35)  # Hero size for targets
        elif node in connected_genes:
            node_color.append("white")
            node_size.append(25)  # Prominent neighbors
        else:
            node_color.append("#D3D3D3")
            node_size.append(15)   # Background nodes

    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        name='Genes',
        marker=dict(
            showscale=False,
            color=node_color,
            size=node_size,
            line=dict(color='rgba(0,0,0,0.8)', width=1)
        )
    )

    fig = go.Figure(data=[normal_edge_trace, target_edge_trace, node_trace],
                 layout=go.Layout(
                    title={'text': f'Interactive 3D Network ({len(found_targets)} targets)', 'font': {'size': 16}},
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=0, l=0, r=0, t=40),
                    height=800, # Square-like height
                    scene=dict(
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title='', visible=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title='', visible=False),
                        zaxis=dict(showgrid=False, zeroline=False, showticklabels=False, title='', visible=False),
                        bgcolor='rgba(0,0,0,0)'
                    ),
                    paper_bgcolor='rgba(0,0,0,0)',
                    plot_bgcolor='rgba(0,0,0,0)'
                ))
    
    return fig
