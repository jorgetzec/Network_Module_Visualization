FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim

WORKDIR /app

COPY pyproject.toml /app/
COPY uv.lock /app/
COPY README.md /app/
COPY .streamlit /app/.streamlit
COPY src /app/src
COPY streamlit_app.py /app/

RUN uv sync --frozen && uv pip install --python /app/.venv/bin/python scipy

EXPOSE 8501

CMD ["uv", "run", "streamlit", "run", "streamlit_app.py", "--server.address=0.0.0.0", "--server.port=8501"]
