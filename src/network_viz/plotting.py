from __future__ import annotations

from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx


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
    labels = nx.draw_networkx_labels(graph, pos, labels=label_dict, font_size=4, font_color="black")

    for label in labels.values():
        label.set_zorder(6)
    for target in found_targets:
        if target in labels:
            labels[target].set_fontsize(10)
            labels[target].set_color("black")
            labels[target].set_weight("bold")
    for connected in connected_nodes:
        if connected in labels:
            labels[connected].set_fontsize(6)

    target_str = ", ".join(found_targets) if found_targets else "No target genes"
    plt.title(f"Gene network highlighting {target_str} ({module_color} module)", fontsize=14)
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
    labels = nx.draw_networkx_labels(graph_zoom, pos_scaled, labels=label_dict, font_size=8, font_color="black")
    for label in labels.values():
        label.set_zorder(6)
    for target in valid_targets:
        if target in labels:
            labels[target].set_fontsize(14)
            labels[target].set_color("black")
            labels[target].set_weight("bold")
    for connected in connected_nodes:
        if connected in labels:
            labels[connected].set_fontsize(10)

    target_str = ", ".join(valid_targets)
    plt.title(f"Zoomed view: {target_str} and neighbors ({module_color} module)", fontsize=16, fontweight="bold")
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
