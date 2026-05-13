# syntax=docker/dockerfile:1.7

# ---- Builder: install deps into a wheel cache --------------------------------
FROM python:3.11-slim AS builder

ENV PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libxrender1 \
        libxext6 \
        libfontconfig1 \
        poppler-utils \
        libcairo2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY requirements-finder.txt .
RUN pip install --prefix=/install -r requirements-finder.txt

# ---- Runtime ----------------------------------------------------------------
FROM python:3.11-slim AS runtime

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false \
    STREAMLIT_SERVER_ENABLE_CORS=false \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=true

RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 \
        libxext6 \
        libfontconfig1 \
        poppler-utils \
        libcairo2 \
        ca-certificates \
        curl \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --system app && useradd --system --gid app --home /home/app --create-home app

COPY --from=builder /install /usr/local

WORKDIR /app
COPY . /app
RUN chown -R app:app /app
USER app

EXPOSE 8501

# Streamlit ships its own /healthz endpoint at /_stcore/health.
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health || exit 1

CMD ["python", "-m", "streamlit", "run", "streamlit_finder.py", \
     "--server.port=8501", "--server.address=0.0.0.0"]
