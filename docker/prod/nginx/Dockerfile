# Nginx Dockerfile for static file serving and reverse proxy
# This builds the frontend and serves it via Nginx

# Stage 1: Build frontend with Bun (reused from zndraw)
FROM oven/bun:1 AS frontend-builder

WORKDIR /build

# Copy frontend package files
COPY app/package.json app/bun.lock* ./

# Install frontend dependencies
RUN bun install --frozen-lockfile

# Copy frontend source and build configuration
COPY app/ ./

# Build frontend (outputs to /src/zndraw/static per vite.config.ts)
RUN bun run vite build

# ============================================================================
# Stage 2: Nginx production server
FROM nginx:1.25-alpine

# Copy custom nginx configuration
COPY docker/nginx/nginx.conf /etc/nginx/nginx.conf

# Copy built frontend static files from frontend-builder stage
# Vite builds to /src/zndraw/static (see app/vite.config.ts)
COPY --from=frontend-builder /src/zndraw/static /usr/share/nginx/html

# Create nginx user and set permissions
RUN chown -R nginx:nginx /usr/share/nginx/html && \
    chmod -R 755 /usr/share/nginx/html

# Expose HTTP port
EXPOSE 80

# Health check
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD wget --no-verbose --tries=1 --spider http://localhost/health || exit 1

# Start nginx in foreground
CMD ["nginx", "-g", "daemon off;"]
