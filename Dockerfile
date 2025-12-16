# Multi-stage Dockerfile for ZnDraw - Optimized for minimal size
#
# Stage 1: Build everything (Python + Frontend) with uv
# - hatch_build.py automatically builds frontend during `uv sync`
# - hatch-vcs gets version from git tags
# - Frontend version is injected via VITE_APP_VERSION
#
# Stage 2: Minimal runtime (no build tools)

# ============================================================================
# Stage 1: Build with uv (includes frontend build via hatch_build.py)
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS builder

# Install git (for hatch-vcs) and bun (for frontend build)
RUN apt-get update && \
    apt-get install -y --no-install-recommends git curl unzip && \
    curl -fsSL https://bun.sh/install | bash && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/root/.bun/bin:$PATH"
ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy

WORKDIR /app

# Copy dependency files for layer caching
COPY pyproject.toml uv.lock* hatch_build.py ./

# Sync dependencies without installing project (for caching)
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-install-project --extra full

# Copy all source code
COPY src/ ./src/
COPY app/ ./app/
COPY LICENSE README.md ./

# Install project with frontend build (bind-mount .git for version detection)
# hatch_build.py will:
# 1. Get version from git via hatch-vcs
# 2. Build frontend with bun, injecting VITE_APP_VERSION
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=source=.git,target=.git,type=bind \
    uv sync --frozen --no-dev --extra full

# ============================================================================
# Stage 2: Minimal runtime image
FROM python:3.12-slim-bookworm AS runtime

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

# Copy only the virtual environment and application from builder
COPY --from=builder --chown=appuser:appuser /app/.venv /app/.venv
COPY --from=builder --chown=appuser:appuser /app/src /app/src

# Set environment
ENV PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    ZNDRAW_STORAGE_PATH=/app/data/zndraw-data \
    ZNDRAW_UPLOAD_TEMP=/tmp/zndraw_uploads \
    ZNDRAW_SERVER_HOST=0.0.0.0 \
    ZNDRAW_REDIS_URL=redis://redis:6379 \
    OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

USER appuser

EXPOSE 5000

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:5000').read()" || exit 1

CMD ["zndraw", "--no-browser"]
