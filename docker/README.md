# ZnDraw Docker Deployment

## Deployment Options

### [Standalone](standalone/) - Simple single-instance deployment

```bash
cd standalone
docker compose up -d
```

Access at http://localhost:5000

Best for: Personal use, small teams, development.

### [Production](production/) - Horizontal scaling with load balancer

```bash
cd production
docker compose up -d
```

Access at http://localhost

Best for: Production, high load, multiple concurrent users.

## Comparison

| Feature | Standalone | Production |
|---------|------------|--------|
| ZnDraw instances | 1 | 3+ (configurable) |
| Celery workers | 1 | 2+ (configurable) |
| Load balancer | No | Nginx |
| Port | 5000 | 80 |
| Complexity | Simple | Moderate |
