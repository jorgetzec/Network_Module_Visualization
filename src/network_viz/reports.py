from __future__ import annotations

from datetime import datetime
from pathlib import Path

import networkx as nx
import pandas as pd


def build_suffix(weight_threshold: float, module_color: str) -> str:
    current_date = datetime.now().strftime("%Y%m%d")
    return f"_w{weight_threshold}_{module_color}_{current_date}"


def generate_network_statistics_df(graph: nx.Graph) -> pd.DataFrame:
    degree_dict = dict(graph.degree())
    rows = []
    for gene, degree in sorted(degree_dict.items(), key=lambda x: x[1], reverse=True):
        neighbors = list(graph.neighbors(gene))
        rows.append(
            {
                "gene": gene,
                "degree": degree,
                "connections": ",".join(neighbors),
            }
        )
    return pd.DataFrame(rows)


def write_network_statistics(
    graph: nx.Graph, output_prefix: str, weight_threshold: float, module_color: str
) -> Path:
    suffix = build_suffix(weight_threshold, module_color)
    out_path = Path(f"{output_prefix}_network_statistics{suffix}.txt")

    lines = [
        "=== NETWORK STATISTICS ===",
        f"Weight threshold: {weight_threshold}",
        f"Module color: {module_color}",
        f"Analysis date: {datetime.now().strftime('%Y%m%d')}",
        f"Number of nodes: {graph.number_of_nodes()}",
        f"Number of edges: {graph.number_of_edges()}",
    ]
    if graph.number_of_nodes() > 0:
        avg_degree = sum(dict(graph.degree()).values()) / graph.number_of_nodes()
        lines.append(f"Average connections per node: {avg_degree:.2f}")
        lines.append("")
        lines.append("=== ALL GENES SORTED BY CONNECTIONS (DESCENDING) ===")
        lines.append("Gene\tDegree\tConnections")
        lines.append("-" * 80)
        for _, row in generate_network_statistics_df(graph).iterrows():
            lines.append(f"{row['gene']}\t{row['degree']}\t{row['connections']}")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path


def write_connected_genes_report(
    graph: nx.Graph,
    found_targets: list[str],
    connected_genes: list[str],
    output_prefix: str,
    weight_threshold: float,
    module_color: str,
) -> Path | None:
    if not found_targets:
        return None

    suffix = build_suffix(weight_threshold, module_color)
    out_path = Path(f"{output_prefix}_connected_genes{suffix}.txt")
    lines = [
        "=== CONNECTED GENES INFORMATION ===",
        f"Weight threshold: {weight_threshold}",
        f"Module color: {module_color}",
        f"Analysis date: {datetime.now().strftime('%Y%m%d')}",
        f"Target genes found in network: {', '.join(found_targets)}",
        "",
        f"Total unique genes connected to ANY target gene: {len(connected_genes)}",
        "",
        "Gene\tDegree\tConnected to target(s)",
        "-" * 80,
    ]

    for gene in connected_genes:
        if gene not in graph:
            continue
        degree = graph.degree(gene)
        connected_to = [target for target in found_targets if target in graph and gene in graph.neighbors(target)]
        lines.append(f"{gene}\t{degree}\t{','.join(connected_to)}")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    return out_path
