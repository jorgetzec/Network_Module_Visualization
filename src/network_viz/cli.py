from __future__ import annotations

import argparse

from .core import load_genes, run_analysis
from .plotting import draw_highlighted_network, save_network_figure
from .reports import write_connected_genes_report, write_network_statistics


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Network analysis tool for gene modules with customizable highlighting",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  uv run network-viz --nodes data/nodes.txt --edges data/edges.txt --weight 0.17 --genes GPD2,GPD3 --module green
  uv run network-viz --nodes data/nodes.txt --edges data/edges.txt --weight 0.2 --genes-file genes_list.txt --module blue
        """,
    )
    parser.add_argument("--nodes", required=True, help="Path to nodes file (tab-separated)")
    parser.add_argument("--edges", required=True, help="Path to edges file (tab-separated)")
    parser.add_argument("--weight", type=float, required=True, help="Weight threshold for filtering edges")
    parser.add_argument("--genes", help="Comma-separated list of genes to highlight")
    parser.add_argument("--genes-file", help="Path to file containing genes to highlight (one per line)")
    parser.add_argument("--module", default="red", help="Color for highlighting genes and connections")
    parser.add_argument("--output-prefix", default="network_analysis", help="Prefix for output files")
    parser.add_argument("--min-degree", type=int, default=1, help="Minimum degree for node filtering")
    parser.add_argument(
        "--show-only-target-labels",
        action="store_true",
        help="Show labels only for target genes and their connected genes",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    target_genes = load_genes(args.genes, args.genes_file)
    analysis = run_analysis(
        nodes_file=args.nodes,
        edges_file=args.edges,
        weight_threshold=args.weight,
        target_genes=target_genes,
        min_degree=args.min_degree,
    )

    stats_file = write_network_statistics(analysis.graph, args.output_prefix, args.weight, args.module)
    connected_file = write_connected_genes_report(
        analysis.graph,
        analysis.found_targets,
        analysis.connected_genes,
        args.output_prefix,
        args.weight,
        args.module,
    )

    fig, _ = draw_highlighted_network(
        analysis.graph,
        analysis.found_targets,
        analysis.connected_genes,
        module_color=args.module,
        show_only_target_labels=args.show_only_target_labels,
    )
    network_file = save_network_figure(fig, args.output_prefix, args.weight, args.module, highlighted=True)

    print(f"Nodes: {analysis.graph.number_of_nodes()}")
    print(f"Edges: {analysis.graph.number_of_edges()}")
    print(f"Found targets: {analysis.found_targets}")
    print(f"Connected genes: {len(analysis.connected_genes)}")
    print(f"Network statistics: {stats_file}")
    if connected_file:
        print(f"Connected genes report: {connected_file}")
    print(f"Network figure: {network_file}")


if __name__ == "__main__":
    main()
