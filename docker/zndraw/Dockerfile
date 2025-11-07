# Multi-stage Dockerfile for ZnDraw Application
# Stage 1: Build frontend with Bun
FROM oven/bun:1 AS frontend-builder

WORKDIR /build

# Copy frontend package files
COPY app/package.json app/bun.lock* ./

# Install frontend dependencies
RUN bun install --frozen-lockfile

# Copy frontend source and build configuration
COPY app/ ./

# Build frontend (outputs to /src/zndraw/static per vite.config.ts)
# Skip tsc type checking for faster builds - vite will handle transpilation
RUN bun run vite build

# ============================================================================
# Stage 2: Python application with uv
FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS app

# Install runtime dependencies for RDKit Draw module
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libx11-6 \
    libexpat1 \
    libfreetype6 \
    libpng16-16 \
    && rm -rf /var/lib/apt/lists/*

# Create app user for security
RUN useradd -m -u 1000 -s /bin/bash appuser && \
    mkdir -p /app/data /tmp/zndraw_uploads && \
    chown -R appuser:appuser /app /tmp/zndraw_uploads

# Set working directory
WORKDIR /app

# Copy Python dependency files first (for better layer caching)
COPY --chown=appuser:appuser pyproject.toml uv.lock* LICENSE README.md ./

# Install dependencies into virtual environment as appuser
USER appuser
RUN --mount=type=cache,target=/home/appuser/.cache/uv,uid=1000,gid=1000 \
    uv sync --frozen --no-dev --no-install-project

# Copy application source
COPY --chown=appuser:appuser src/ ./src/

# Copy built frontend from frontend-builder stage
# Vite builds to /src/zndraw/static (see app/vite.config.ts)
COPY --from=frontend-builder --chown=appuser:appuser /src/zndraw/static/ ./src/zndraw/static/

# Install the project itself
RUN --mount=type=cache,target=/home/appuser/.cache/uv,uid=1000,gid=1000 \
    uv sync --frozen --no-dev

# Set environment variables
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    # Storage paths
    ZNDRAW_STORAGE_PATH=/app/data/zndraw-data \
    ZNDRAW_UPLOAD_TEMP=/tmp/zndraw_uploads \
    # Server configuration
    ZNDRAW_SERVER_HOST=0.0.0.0 \
    # Redis URL (override in docker-compose)
    ZNDRAW_REDIS_URL=redis://redis:6379 \
    # Fix for Celery on Apple Silicon hosts
    OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:5000').read()" || exit 1

# Default command: Use zndraw CLI with eventlet server
# Configuration via environment variables (see docker-compose.yaml)
CMD ["uv", "run", "zndraw", "--no-browser"]
