from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import pandas as pd


@dataclass
class AnalysisResult:
    graph: nx.Graph
    found_targets: list[str]
    connected_genes: list[str]


def load_genes(genes_arg: str | None = None, genes_file_arg: str | None = None) -> list[str]:
    genes: list[str] = []

    if genes_arg:
        genes.extend([gene.strip() for gene in genes_arg.split(",") if gene.strip()])

    if genes_file_arg:
        genes_file = Path(genes_file_arg)
        if not genes_file.exists():
            raise FileNotFoundError(f"Genes file not found: {genes_file_arg}")
        file_genes = [line.strip() for line in genes_file.read_text(encoding="utf-8").splitlines() if line.strip()]
        genes.extend(file_genes)

    unique_genes: list[str] = []
    seen: set[str] = set()
    for gene in genes:
        if gene not in seen:
            unique_genes.append(gene)
            seen.add(gene)
    return unique_genes


def load_network(nodes_file: str, edges_file: str, weight_threshold: float, min_degree: int = 1) -> nx.Graph:
    nodes_path = Path(nodes_file)
    edges_path = Path(edges_file)
    if not nodes_path.exists():
        raise FileNotFoundError(f"Nodes file not found: {nodes_file}")
    if not edges_path.exists():
        raise FileNotFoundError(f"Edges file not found: {edges_file}")

    nodes_df = pd.read_csv(nodes_path, sep="\t")
    edges_df = pd.read_csv(edges_path, sep="\t")
    return load_network_from_dataframes(nodes_df, edges_df, weight_threshold, min_degree=min_degree)


def load_network_from_dataframes(
    nodes_df: pd.DataFrame, edges_df: pd.DataFrame, weight_threshold: float, min_degree: int = 1
) -> nx.Graph:
    required_nodes_cols = {"nodeName", "altName"}
    required_edges_cols = {"fromNode", "toNode", "weight"}
    missing_nodes = required_nodes_cols - set(nodes_df.columns)
    missing_edges = required_edges_cols - set(edges_df.columns)
    if missing_nodes:
        raise ValueError(f"Missing required node columns: {sorted(missing_nodes)}")
    if missing_edges:
        raise ValueError(f"Missing required edge columns: {sorted(missing_edges)}")

    edges_df = edges_df[edges_df["weight"] >= weight_threshold]

    graph = nx.Graph()
    node_attr_col = "nodeAttr[nodesPresent, ]"

    for _, row in nodes_df.iterrows():
        graph.add_node(
            row["nodeName"],
            altName=row["altName"],
            nodeAttr=row[node_attr_col] if node_attr_col in nodes_df.columns else None,
        )

    for _, row in edges_df.iterrows():
        graph.add_edge(row["fromNode"], row["toNode"], weight=float(row["weight"]))

    if min_degree > 1:
        nodes_with_min_degree = [node for node, degree in dict(graph.degree()).items() if degree >= min_degree]
        graph = graph.subgraph(nodes_with_min_degree).copy()

    return graph


def find_connected_genes(graph: nx.Graph, target_genes: list[str]) -> tuple[list[str], list[str]]:
    connected_genes: set[str] = set()
    found_targets: list[str] = []
    for target in target_genes:
        if target in graph:
            found_targets.append(target)
            connected_genes.update(graph.neighbors(target))
    return found_targets, sorted(connected_genes)


def run_analysis(
    nodes_file: str,
    edges_file: str,
    weight_threshold: float,
    target_genes: list[str],
    min_degree: int = 1,
) -> AnalysisResult:
    graph = load_network(nodes_file, edges_file, weight_threshold=weight_threshold, min_degree=min_degree)
    found_targets, connected_genes = find_connected_genes(graph, target_genes)
    return AnalysisResult(graph=graph, found_targets=found_targets, connected_genes=connected_genes)
