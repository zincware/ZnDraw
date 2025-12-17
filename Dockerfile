# Multi-stage Dockerfile for ZnDraw - Optimized for minimal size
#
# Stage 1: Build wheel with correct version
# - SETUPTOOLS_SCM_PRETEND_VERSION overrides dirty version detection
# - hatch_build.py builds frontend during wheel build
#
# Stage 2: Minimal runtime - just pip install the wheel

# ============================================================================
# Stage 1: Build wheel with uv (includes frontend build via hatch_build.py)
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS builder

# Version override for Docker builds (avoids dirty version detection)
# See: https://github.com/pypa/setuptools-scm/issues/77
ARG SETUPTOOLS_SCM_PRETEND_VERSION
ENV SETUPTOOLS_SCM_PRETEND_VERSION=${SETUPTOOLS_SCM_PRETEND_VERSION}

# Install git (for hatch-vcs fallback) and bun (for frontend build)
RUN apt-get update && \
    apt-get install -y --no-install-recommends git curl unzip && \
    curl -fsSL https://bun.sh/install | bash && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/root/.bun/bin:$PATH"

WORKDIR /build

# Copy source and build wheel
COPY . /src
RUN --mount=type=cache,target=/root/.cache/uv \
    cd /src && uv build --wheel --out-dir /build/dist

# ============================================================================
# Stage 2: Minimal runtime image
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS runtime

# Install only runtime dependencies for RDKit Draw module
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libx11-6 \
    libexpat1 \
    libfreetype6 \
    libpng16-16 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Create non-root user
RUN groupadd --system --gid 999 appuser && \
    useradd --system --gid 999 --uid 999 --create-home appuser && \
    mkdir -p /app/data /tmp/zndraw_uploads && \
    chown -R appuser:appuser /app /tmp/zndraw_uploads

WORKDIR /app

# Copy wheel and install
COPY --from=builder /build/dist/*.whl /tmp/
RUN WHEEL=$(ls /tmp/*.whl) && uv pip install --system --no-cache-dir "${WHEEL}[full]" && rm /tmp/*.whl

# Set environment
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    ZNDRAW_STORAGE_PATH=/app/data/zndraw-data \
    ZNDRAW_UPLOAD_TEMP=/tmp/zndraw_uploads \
    ZNDRAW_SERVER_HOST=0.0.0.0 \
    ZNDRAW_REDIS_URL=redis://redis:6379

USER appuser

EXPOSE 5000

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:5000').read()" || exit 1

CMD ["zndraw", "--no-browser"]
