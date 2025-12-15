# ZnDraw Standalone Deployment

Simple single-instance deployment with Redis and Celery worker.

## Quick Start

```bash
docker compose up -d
```

Access at http://localhost:5000

## With Custom Templates

Mount a local templates directory:

```yaml
services:
  zndraw:
    volumes:
      - zndraw-data:/app/data
      - ./templates:/app/templates:ro
```

Then run:

```bash
docker compose up -d
```

## Configuration

Edit `docker-compose.yaml` to customize:

| Variable | Description | Default |
|----------|-------------|---------|
| `ZNDRAW_REDIS_URL` | Redis connection | `redis://redis:6379` |
| `ZNDRAW_SERVER_PORT` | Server port | `5000` |
| `FLASK_SECRET_KEY` | Session secret | Change in production! |
| `ZNDRAW_ADMIN_USERNAME` | Admin user | Disabled |
| `ZNDRAW_ADMIN_PASSWORD` | Admin password | Disabled |

## Commands

```bash
# Start
docker compose up -d

# View logs
docker compose logs -f

# Stop
docker compose down

# Stop and remove data
docker compose down -v
```
