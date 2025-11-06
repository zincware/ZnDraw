# ZnDraw Development Docker Setup

Simple single-instance Docker Compose configuration for development and testing.

## When to Use This

**Use this setup for:**
- Local development and testing
- Debugging ZnDraw application behavior
- Testing configuration changes
- Small-scale testing with direct port access

**DO NOT use this for production.** For production deployment, use `/docker/docker-compose.yaml` instead.

## Quick Start

```bash
cd docker/zndraw

# Build and start services
docker compose build
docker compose up -d

# View logs
docker compose logs -f

# Check status
docker compose ps

# Stop services
docker compose down
```

## Access

- **Application**: http://localhost:5000
- **Redis**: localhost:6380 (for monitoring)

## Architecture

```
localhost:5000 → ZnDraw App (Gunicorn)
                      ↓
                  Redis
                      ↓
              Celery Worker
```

**Key differences from production:**
- No nginx reverse proxy
- No load balancing
- Single ZnDraw instance (not scaled)
- Single Celery worker (lower concurrency)
- Direct port exposure
- Debug logging enabled

## Configuration

Edit environment variables in `docker-compose.yaml`:

```yaml
environment:
  - GUNICORN_WORKERS=1
  - GUNICORN_THREADS=50        # Adjust for your workload
  - GUNICORN_LOG_LEVEL=debug   # Use 'info' for less verbosity
  - FLASK_SECRET_KEY=development-secret-key
```

## Testing

This setup is ideal for:
1. Testing Gunicorn configuration
2. Verifying Redis integration
3. Testing Celery tasks
4. Debugging WebSocket connections
5. Testing file uploads

## Monitoring

```bash
# View specific service logs
docker compose logs -f zndraw
docker compose logs -f celery-worker
docker compose logs -f redis

# Check resource usage
docker stats

# Inspect container
docker inspect zndraw-dev-app
```

## Troubleshooting

### Port 5000 already in use

```bash
# Find what's using port 5000
lsof -i :5000

# Change port in docker-compose.yaml
ports:
  - "5001:5000"  # Use different host port
```

### Can't connect to Redis

```bash
# Check Redis is running and healthy
docker compose ps redis
docker compose logs redis

# Test Redis connection
docker compose exec zndraw redis-cli -h redis ping
```

### Frontend not loading

```bash
# Check if static files were copied during build
docker compose exec zndraw ls -la /app/src/zndraw/static/

# Rebuild if needed
docker compose build --no-cache zndraw
```

## Switching to Production

When ready for production:

1. Review your configuration changes
2. Apply them to `/docker/docker-compose.yaml`
3. Follow production deployment checklist in `/docker/README.md`

```bash
# Switch to production setup
cd ../
docker compose build
docker compose up -d
```

## Files

- `docker-compose.yaml` - Development Docker Compose configuration
- `Dockerfile` - ZnDraw application image (shared with production)
- `README.md` - This file

## Further Documentation

- Production setup: `/docker/README.md`
- Deployment guide: `/DEPLOYMENT.md`
- Project documentation: `/README.md`
