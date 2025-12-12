# ZnDraw Docker Deployment

This directory contains Docker configurations for both development and production deployments.

## Directory Structure

```
docker/
├── README.md           # This file
├── dev/                # Development environment with hot reload
│   ├── docker-compose.yaml
│   ├── docker-compose.dev.yaml
│   ├── nginx/
│   └── README.md       # Development setup guide
└── prod/               # Production deployment
    ├── docker-compose.yaml
    ├── docker-compose.with-templates.yaml
    ├── nginx/
    ├── zndraw/
    └── README.md       # Production deployment guide
```

## Quick Start

### Development

For local development with hot reload:

```bash
cd docker/dev
docker compose -f docker-compose.dev.yaml up -d
```

Access at http://localhost:5000

See [dev/README.md](dev/README.md) for full development setup instructions.

### Production

For production deployment:

```bash
cd docker/prod
docker compose up -d
```

Access at http://localhost

See [prod/README.md](prod/README.md) for full production deployment guide.

## Architecture Overview

```
                    ┌─────────────────────────────────────┐
                    │              Nginx                  │
                    │         (Reverse Proxy)             │
                    └──────────────┬──────────────────────┘
                                   │
              ┌────────────────────┼────────────────────┐
              │                    │                    │
              ▼                    ▼                    ▼
        ┌──────────┐        ┌──────────┐        ┌──────────┐
        │  ZnDraw  │        │  ZnDraw  │        │  ZnDraw  │
        │ Replica 1│        │ Replica 2│        │ Replica 3│
        └────┬─────┘        └────┬─────┘        └────┬─────┘
             │                   │                   │
             └───────────────────┼───────────────────┘
                                 │
                    ┌────────────┴────────────┐
                    │                         │
                    ▼                         ▼
              ┌──────────┐             ┌──────────┐
              │  Redis   │             │  Celery  │
              │ (Cache)  │◄────────────│ Workers  │
              └──────────┘             └──────────┘
```

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `ZNDRAW_REDIS_URL` | Redis connection URL | `redis://redis:6379/0` |
| `ZNDRAW_SERVER_HOST` | Server bind address | `0.0.0.0` |
| `ZNDRAW_SERVER_PORT` | Server port | `5000` |
| `ZNDRAW_ADMIN_USERNAME` | Admin username | None |
| `ZNDRAW_ADMIN_PASSWORD` | Admin password | None |
| `ZNDRAW_STORAGE_PATH` | LMDB storage directory | `/data` |
| `ZNDRAW_CELERY_ENABLED` | Enable Celery workers | `true` |
| `ZNDRAW_FILE_BROWSER_ENABLED` | Enable file browser | `false` |
| `FLASK_SECRET_KEY` | Flask secret key | (generate one!) |

### Resource Requirements

**Minimum (Development):**
- 2 CPU cores
- 4 GB RAM
- 10 GB disk space

**Recommended (Production):**
- 8+ CPU cores
- 16+ GB RAM
- 50+ GB disk space (depends on data size)

## Comparison

| Feature | Development | Production |
|---------|-------------|------------|
| Frontend | Vite dev server (hot reload) | Pre-built static files |
| Backend Replicas | 1 | 3+ |
| Celery Workers | 1 | 2+ |
| Log Level | DEBUG | INFO |
| Source Mounts | Yes | No |
| Use Case | Local development | Server deployment |
