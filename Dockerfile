# =============================================================================
# gnomix-1000g — Docker image
# =============================================================================
# Build:
#   docker build -t gnomix-1000g .
#
# Run (downloads, intermediate files and outputs land in ./gnomix-1000g-data).
# Make the folder first: if Docker makes it, root owns it and --user cannot write.
#   mkdir -p gnomix-1000g-data
#   docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)/gnomix-1000g-data:/data" gnomix-1000g
#
# Settings: -e THREADS=8 -e CHROMS="1-22" -e RUN_TRIO_BENCHMARK=0 (see README)
# =============================================================================

# Python 3.9: the pretrained Gnomix pickles need scikit-learn 1.0.1 and
# xgboost 1.1.1 (see requirements.txt).
# Base image pinned by digest (python:3.9-slim-bookworm, October 2026) for reproducible builds.
FROM python:3.9-slim-bookworm@sha256:a02e9c5406c416c504d6c9a1a306ff4080c3173f1008d192f953bd20382a2d5c

LABEL org.opencontainers.image.source="https://github.com/human-genomics/gnomix-1000g"
LABEL org.opencontainers.image.description="Local ancestry (pretrained Gnomix + Gnofix phase correction) for all 3,202 high-coverage 1000 Genomes samples"
LABEL org.opencontainers.image.licenses="MIT"

# Java runs Beagle for the trio phasing benchmark. libgomp1 is needed by xgboost.
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates default-jre-headless libgomp1 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY lib.sh setup_tools.sh /app/
RUN bash /app/setup_tools.sh /opt/tools

COPY main.sh download.sh /app/
COPY gnomix1000g /app/gnomix1000g
COPY reference /app/reference

ENV GNOMIX1000G_PYTHON=python3 \
    PYTHONPATH=/app \
    GNOMIX1000G_TOOLS=/opt/tools \
    GNOMIX1000G_DATA_DIR=/data \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONWARNINGS=ignore \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

WORKDIR /data
ENTRYPOINT ["bash", "/app/main.sh"]
