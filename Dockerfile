# Multi-stage Dockerfile for ZnDraw
# Stage 1: Build frontend with Bun
FROM oven/bun:1 AS frontend-builder

WORKDIR /build

# Copy frontend package files
COPY app/package.json app/bun.lock* ./

# Install frontend dependencies
RUN bun install --frozen-lockfile

# Copy frontend source and build configuration
COPY app/ ./

# Build frontend (outputs to dist/)
# Skip tsc type checking for faster builds - vite will handle transpilation
RUN bun run vite build

# ============================================================================
# Stage 2: Build Python application with uv
FROM ghcr.io/astral-sh/uv:python3.12-bookworm AS builder

# Set working directory
WORKDIR /app

# Copy Python dependency files first (for better layer caching)
COPY pyproject.toml uv.lock* LICENSE README.md ./

# Install dependencies into virtual environment
# Use --no-install-project to cache dependencies separately from source code
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev --no-install-project

# Copy application source
COPY src/ ./src/

# Copy built frontend from previous stage
# Vite outputs to /src/zndraw/static/ in the builder container
COPY --from=frontend-builder /src/zndraw/static/ ./src/zndraw/static/

# Install the project itself
RUN --mount=type=cache,target=/root/.cache/uv \
    uv sync --frozen --no-dev

# ============================================================================
# Stage 3: Production runtime
FROM python:3.12-slim

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
RUN useradd -m -u 1000 -s /bin/bash appuser

# Set working directory
WORKDIR /app

# Copy virtual environment from builder
COPY --from=builder --chown=appuser:appuser /app/.venv /app/.venv

# Copy application code from builder
COPY --from=builder --chown=appuser:appuser /app/src /app/src

# Copy pyproject.toml for metadata
COPY --from=builder --chown=appuser:appuser /app/pyproject.toml /app/pyproject.toml

# Set environment variables
ENV VIRTUAL_ENV=/app/.venv \
    PATH="/app/.venv/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    # Storage paths
    ZNDRAW_STORAGE_PATH=/app/data/zndraw-data.zarr \
    ZNDRAW_UPLOAD_TEMP=/tmp/zndraw_uploads \
    # Server configuration
    ZNDRAW_SERVER_HOST=0.0.0.0 \
    # Redis URL (override in docker-compose)
    ZNDRAW_REDIS_URL=redis://redis:6379 \
    # Fix for Celery on Apple Silicon hosts
    OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES

# Create data directories
RUN mkdir -p /app/data /tmp/zndraw_uploads && \
    chown -R appuser:appuser /app/data /tmp/zndraw_uploads

# Switch to non-root user
USER appuser

# Expose port
EXPOSE 5000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:5000').read()" || exit 1

# Default command (can be overridden in docker-compose)
# Use --no-celery to prevent embedded worker (dedicated celery-worker container handles tasks)
# Use --host with the Docker service name so celery tasks can connect via internal network
CMD ["zndraw", "--host", "zndraw", "--port", "5000", "--no-celery"]
