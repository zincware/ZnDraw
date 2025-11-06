# Docker Production Deployment

This directory contains the production-ready Docker configuration for ZnDraw with separate Dockerfiles for clarity and maintainability.

## Directory Structure

```
docker/
├── README.md              # This file
├── docker-compose.yaml    # Production compose configuration
├── zndraw/
│   └── Dockerfile        # ZnDraw application (Python + Gunicorn)
└── nginx/
    ├── Dockerfile        # Nginx reverse proxy + static files
    └── nginx.conf        # Nginx configuration
```

## Architecture

```
Internet → Nginx (port 80)
             ├─ Static files (served directly)
             ├─ /api/* → Load balanced across ZnDraw replicas
             └─ /socket.io/* → Sticky session routing
                    ↓
            ZnDraw App × 3 replicas
                    ↓
                Redis (pub/sub + Celery broker)
                    ↓
            Celery Workers × 2 replicas
```

## Why Separate Dockerfiles?

**Benefits:**
- **Clarity**: Each service has its own clear, focused Dockerfile
- **No Build Confusion**: No need for `target:` specifications
- **Easier Maintenance**: Changes to one service don't affect others
- **CI/CD Friendly**: Can build and deploy services independently
- **Better Caching**: Docker can cache each service's layers independently

**Trade-off:**
- Frontend is built twice (once for nginx, once for zndraw app)
- This is acceptable for production where build time matters less than clarity

## Quick Start

```bash
cd docker

# Build and start all services
docker compose build
docker compose up -d

# Check status (should see 6 containers)
docker compose ps

# View logs
docker compose logs -f

# Stop all services
docker compose down
```

## Access

- **Application**: http://localhost
- **Health Check**: http://localhost/health
- **Redis** (for monitoring): localhost:6380

## Configuration

### Resource Limits

Default configuration (for 128 core, 500GB server):

| Service | Replicas | Cores/Container | Memory/Container | Total Cores | Total Memory |
|---------|----------|-----------------|------------------|-------------|--------------|
| ZnDraw  | 3        | 16              | 32GB             | 48          | 96GB         |
| Celery  | 2        | 16              | 32GB             | 32          | 64GB         |
| Redis   | 1        | 4               | 8GB              | 4           | 8GB          |
| Nginx   | 1        | 2               | 2GB              | 2           | 2GB          |
| **Total** | **7**   | -               | -                | **86**      | **170GB**    |

Leaves ~40 cores and 330GB for system and other processes.

### Scaling

Edit `docker-compose.yaml`:

```yaml
zndraw:
  deploy:
    replicas: 5  # Increase for more capacity

celery-worker:
  deploy:
    replicas: 4  # Increase for more background processing
```

Then restart:
```bash
docker compose up -d
```

### Environment Variables

Key variables in `docker-compose.yaml`:

```yaml
# Security
- FLASK_SECRET_KEY=change-this-in-production

# Scaling
- GUNICORN_WORKERS=1      # Always 1 for threaded mode
- GUNICORN_THREADS=100    # Adjust per container capacity

# Features
- ZNDRAW_SIMGEN_ENABLED=true
- ZNDRAW_FILE_BROWSER_ENABLED=true
```

## Dockerfile Details

### zndraw/Dockerfile

**Stages:**
1. `frontend-builder`: Builds React frontend with Bun
2. `app`: Python application with uv, includes built frontend

**Key Features:**
- Multi-stage build for efficiency
- Non-root user (appuser) for security
- Virtual environment with uv
- Production extras (gunicorn, eventlet)
- Health checks built-in

### nginx/Dockerfile

**Stages:**
1. `frontend-builder`: Builds React frontend with Bun (same as zndraw)
2. Production: Nginx with static files and configuration

**Key Features:**
- Lightweight alpine-based nginx
- Custom configuration
- Static file serving
- Reverse proxy with load balancing
- Sticky sessions for WebSockets

## Build Process

Both Dockerfiles use the same frontend builder stage, so the frontend is built twice. This is intentional for clarity:

```
zndraw build: Bun → Build Frontend → Copy to /app → Python App
nginx build:  Bun → Build Frontend → Copy to /usr/share/nginx/html
```

## Development vs Production

**This is for production.** For development, use the root `docker-compose.yaml` or run services directly:

```bash
# Development (root directory)
docker compose up

# Production (docker directory)
cd docker
docker compose up
```

## Monitoring

```bash
# View resource usage
docker stats

# View specific service logs
docker compose logs -f zndraw
docker compose logs -f celery-worker
docker compose logs -f nginx

# Check health status
docker compose ps

# Inspect specific container
docker inspect $(docker compose ps -q zndraw | head -1)
```

## Troubleshooting

### Services Not Starting

```bash
# Check logs
docker compose logs

# Check individual service
docker compose logs zndraw
docker compose logs nginx
```

### Nginx Can't Connect to Backend

```bash
# Verify zndraw is running and healthy
docker compose ps | grep zndraw

# Check nginx logs
docker compose logs nginx | grep "upstream"
```

### Celery Workers Not Starting

```bash
# Check celery logs
docker compose logs celery-worker

# Verify PATH is set correctly
docker compose exec celery-worker sh -c "echo \$PATH"
```

### Frontend Not Loading

```bash
# Check if static files exist in nginx
docker compose exec nginx ls -la /usr/share/nginx/html

# Check nginx error logs
docker compose logs nginx | grep "error"
```

## Production Checklist

Before deploying to production:

- [ ] Change `FLASK_SECRET_KEY` to a random value
- [ ] Set `ZNDRAW_ADMIN_USERNAME` and `ZNDRAW_ADMIN_PASSWORD` if needed
- [ ] Review and adjust resource limits for your server
- [ ] Review and adjust replica counts for your load
- [ ] Configure SSL/TLS (external load balancer or nginx)
- [ ] Set up log aggregation (e.g., ELK, Grafana Loki)
- [ ] Set up monitoring and alerting (e.g., Prometheus + Grafana)
- [ ] Configure automated backups for data volumes
- [ ] Test failover by stopping individual containers
- [ ] Load test under expected peak conditions
- [ ] Document your specific configuration

## Data Storage

Production deployment uses Docker named volumes for data persistence:

```yaml
volumes:
  redis-data:    # Redis persistence
  zndraw-data:   # ZnDraw application data (zarr files)
  upload-temp:   # Temporary upload storage
```

### Managing Volumes

```bash
# List all volumes
docker volume ls

# Inspect volume details
docker volume inspect docker_zndraw-data

# Backup volume data
docker run --rm -v docker_zndraw-data:/data -v $(pwd):/backup \
  alpine tar czf /backup/zndraw-data-backup.tar.gz -C /data .

# Restore volume data
docker run --rm -v docker_zndraw-data:/data -v $(pwd):/backup \
  alpine sh -c "cd /data && tar xzf /backup/zndraw-data-backup.tar.gz"

# Access volume data (for debugging)
docker run --rm -it -v docker_zndraw-data:/data alpine sh
```

### Volume Locations

On the host system, Docker stores volumes in:
- Linux: `/var/lib/docker/volumes/`
- Mac: `~/Library/Containers/com.docker.docker/Data/vms/0/`
- Windows: `C:\ProgramData\Docker\volumes\`

**Note:** Use Docker commands to manage volumes, not direct filesystem access.

## Further Documentation

- See root `/DEPLOYMENT.md` for detailed deployment guide
- See root `/README.md` for project documentation
- See `nginx/nginx.conf` for nginx configuration details
- See each Dockerfile for build process details
