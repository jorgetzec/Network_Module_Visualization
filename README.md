# Network Module Visualization

Network Module Visualization is a dual-interface tool for gene co-expression networks (WGCNA/Cytoscape/VisANT-like exports).

- Python CLI for reproducible batch analysis
- Streamlit web app for interactive exploration
- Full-network and target-centered zoom visualizers

## Features

- Load `nodes` and `edges` tab-separated files
- Filter by edge weight and minimum node degree
- Highlight target genes and their direct neighbors
- Inspect two views:
  - Full filtered network
  - Zoomed 1-hop target subnetwork
- Export analysis outputs (`CSV`, `PDF`)

## Expected Input Format

### `edges` file (TSV)

Required columns:

- `fromNode`
- `toNode`
- `weight`

Optional columns:

- `direction`
- `fromAltName`
- `toAltName`

### `nodes` file (TSV)

Required columns:

- `nodeName`
- `altName`

Optional column:

- `nodeAttr[nodesPresent, ]`

## Quick Start (Docker-first)

This is the recommended way if you want to avoid creating a local `.venv`.

1. Generate or refresh lockfile (when dependencies change):

```bash
uv lock
```

2. Build image:

```bash
docker build -t network-viz .
```

3. Run Streamlit app:

```bash
docker run --rm -p 8501:8501 network-viz
```

4. Open:

```text
http://localhost:8501
```

## Run CLI in Docker

Example using local `data/` mounted as `/data`:

```bash
docker run --rm -v "${PWD}/data:/data" --entrypoint uv network-viz run network-viz --nodes /data/nodes-MODULE_salmon.txt --edges /data/edges-MODULE_salmon.txt --weight 0.2 --genes HLM1,THO2 --module salmon
```

## Local Development with `uv`

```bash
uv sync
uv run network-viz --help
uv run streamlit run streamlit_app.py
```

By default, `uv sync` creates `.venv` in the project folder.

To place the environment outside the repository (PowerShell):

```powershell
$env:UV_PROJECT_ENVIRONMENT="$env:USERPROFILE\.virtualenvs\network-viz"
uv sync
```

## Streamlit Upload Limit

Configured in `.streamlit/config.toml`:

- `server.maxUploadSize = 1024` (MB)

## Project Structure

```text
.
├── Dockerfile
├── pyproject.toml
├── streamlit_app.py
└── src/network_viz/
    ├── app.py
    ├── cli.py
    ├── core.py
    ├── plotting.py
    └── reports.py
```
