# Multi-stage Dockerfile for ZnDraw - Optimized for minimal size
# Stage 1: Build frontend with Bun
FROM oven/bun:1-slim AS frontend-builder

WORKDIR /build

# Copy frontend package files first for layer caching
COPY app/package.json app/bun.lock* ./

# Install frontend dependencies
RUN bun install --frozen-lockfile

# Copy frontend source
COPY app/ ./

# Build frontend (outputs to /src/zndraw/static)
# Skip tsc - vite handles transpilation
RUN bun vite build

# ============================================================================
# Stage 2: Build Python environment with uv
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS python-builder

# Install git for hatch-vcs version detection
RUN apt-get update && apt-get install -y --no-install-recommends git && rm -rf /var/lib/apt/lists/*

ENV UV_COMPILE_BYTECODE=1 UV_LINK_MODE=copy

WORKDIR /app

# Copy dependency files for layer caching
COPY pyproject.toml uv.lock* ./

# Sync dependencies without installing project (for caching)
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-install-project

# Copy source code
COPY src/ ./src/
COPY LICENSE README.md hatch_build.py ./

# Copy built frontend from frontend-builder
COPY --from=frontend-builder /src/zndraw/static/ ./src/zndraw/static/

# Install project (bind-mount .git for hatch-vcs version detection)
RUN --mount=type=cache,target=/root/.cache/uv \
    --mount=source=.git,target=.git,type=bind \
    uv sync --frozen --no-dev

# ============================================================================
# Stage 3: Minimal runtime image (no uv, no bun, no build tools)
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
COPY --from=python-builder --chown=appuser:appuser /app/.venv /app/.venv
COPY --from=python-builder --chown=appuser:appuser /app/src /app/src

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
