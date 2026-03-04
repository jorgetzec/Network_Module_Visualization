# Network Module Visualization

[![Python](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
[![Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://networkmodule-visualization.streamlit.app/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Network Module Visualization** is a dual-interface tool designed for analyzing and visualizing gene co-expression networks (WGCNA, Cytoscape, or VisANT-like exports). It features a robust 3D interactive engine and automated reporting.

- **Python CLI** for reproducible batch analysis.
* **Streamlit Web App** for interactive 2D and 3D exploration.
- **Targeted Analysis**: Highlight specific genes and their neighborhoods.
- **Docker Ready**: One-command deployment.

---

## Features

- **Multi-dimensional Visualization**: Toggle between **2D Maps** and **3D Spherical** interactive views.
- **Robust Filtering**: Filter by edge weight and minimum node degree iteratively.
- **High-Quality Exports**: Generate publication-ready PDFs and detailed CSV reports.
- **Smart Color Detection**: Automatic suggestion of module colors based on node attributes.
- **Deep Zoom**: Focus on the 1-hop neighborhood of your target genes.

## CLI Usage

The project includes a powerful command-line interface for batch processing.

### Command Reference

| Argument | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--nodes` | string | **Required** | Path to the nodes tab-separated file. |
| `--edges` | string | **Required** | Path to the edges tab-separated file. |
| `--weight` | float | **Required** | Weight threshold for filtering edges (0.0 to 1.0). |
| `--genes` | string | None | Comma-separated list of gene IDs to highlight. |
| `--genes-file` | string | None | Path to a file containing gene IDs (one per line). |
| `--module` | string | `red` | Color name or hex code for highlighting. |
| `--min-degree` | int | `1` | Minimum degree required for a node to be kept. |
| `--output-prefix` | string | `network_analysis` | Prefix for generated reports and plots. |
| `--show-only-target-labels` | flag | `False` | If set, skip labels for non-related genes. |

### Example
```bash
uv run network-viz --nodes data/nodes.txt --edges data/edges.txt --weight 0.2 --genes HLM1,THO2 --module salmon
```

---

## Quick Start (Docker)

Recommended for a clean environment without local dependencies.

1. **Build the image**:
   ```bash
   docker build -t network-viz .
   ```

2. **Run the Streamlit app**:
   ```bash
   docker run --rm -p 8501:8501 network-viz
   ```
   *Access it at: [http://localhost:8501](http://localhost:8501)*

3. **Run the CLI**:
   ```bash
   docker run --rm -v "${PWD}/data:/data" --entrypoint uv network-viz run network-viz --nodes /data/nodes.txt --edges /data/edges.txt --weight 0.15
   ```

---

## Local Development

Requires [uv](https://github.com/astral-sh/uv).

```bash
# Install dependencies
uv sync

# Run the app
uv run streamlit run streamlit_app.py

# Run tests
uv run pytest
```

---

## Data Formats

### Edges File (TSV)
Required columns: `fromNode`, `toNode`, `weight`.
Optional: `direction`, `fromAltName`, `toAltName`.

### Nodes File (TSV)
Required columns: `nodeName`, `altName`.
Supports custom attributes for automatic color detection (e.g., `nodeAttr[nodesPresent, ]`).

---

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details (or assume MIT for open-source use).

Created with focus on biological network clarity and interactive depth.
