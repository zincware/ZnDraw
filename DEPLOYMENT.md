# ZnDraw Production Deployment Guide

## Architecture Overview

The production deployment uses a scalable, multi-container architecture optimized for high-performance servers:

```
Internet → Nginx (port 80)
             ├─ Static files (served directly from memory)
             ├─ /api/* → Load balanced across ZnDraw app replicas
             └─ /socket.io/* → Sticky session routing to ZnDraw apps
                    ↓
            ZnDraw App × 3 replicas (300 concurrent connections)
                    ↓
                Redis (pub/sub + state + Celery broker)
                    ↓
            Celery Workers × 2 replicas (16 parallel tasks)
```

### Components

**Nginx (1 container)**
- Serves static frontend files directly (faster than app servers)
- Load balances API requests across app replicas
- Sticky sessions for WebSocket connections (critical!)
- Resource limits: 2 CPU cores, 2GB RAM

**ZnDraw App (3 replicas)**
- Each replica: 100 threads, handles 100 concurrent connections
- Total capacity: 300 concurrent connections
- Flask-SocketIO with Redis message queue for multi-replica coordination
- Resource limits per replica: 16 CPU cores, 32GB RAM

**Celery Workers (2 replicas)**
- Each replica: 8 concurrent tasks
- Total capacity: 16 parallel background tasks
- Shared storage with app replicas
- Resource limits per replica: 16 CPU cores, 32GB RAM

**Redis (1 container)**
- SocketIO pub/sub messaging (coordinates WebSocket across app replicas)
- Celery task queue
- Application state storage
- Persistent data with volume
- Resource limits: 4 CPU cores, 8GB RAM

### Total Resource Usage (Default Configuration)

- **ZnDraw App**: 3 replicas × 16 cores × 32GB = 48 cores, 96GB
- **Celery Workers**: 2 replicas × 16 cores × 32GB = 32 cores, 64GB
- **Redis**: 4 cores, 8GB
- **Nginx**: 2 cores, 2GB
- **Total**: ~86 cores, ~170GB (leaves headroom on 128 core, 500GB server)

## Quick Start

### 1. Build and Deploy

```bash
# Stop existing containers
docker compose down

# Build all images (nginx, app, workers)
docker compose build

# Start all services with scaling
docker compose up -d

# Check status
docker compose ps
```

### 2. View Logs

```bash
# All services
docker compose logs -f

# Specific service
docker compose logs -f nginx
docker compose logs -f zndraw
docker compose logs -f celery-worker
docker compose logs -f redis
```

### 3. Access the Application

```bash
# Frontend and API
http://localhost

# Health check
curl http://localhost/health
```

## Scaling Configuration

### Adjusting Replica Counts

Edit `docker-compose.yaml`:

```yaml
# For light load (1-10 users)
zndraw:
  deploy:
    replicas: 2  # 200 concurrent connections

celery-worker:
  deploy:
    replicas: 1  # 8 parallel tasks

# For heavy load (50+ users)
zndraw:
  deploy:
    replicas: 5  # 500 concurrent connections

celery-worker:
  deploy:
    replicas: 4  # 32 parallel tasks
```

After changing replicas:
```bash
docker compose up -d --scale zndraw=5 --scale celery-worker=4
```

### Adjusting Threads per Container

Edit environment variables in `docker-compose.yaml`:

```yaml
zndraw:
  environment:
    # More threads = more concurrent connections per container
    # But also more memory usage
    - GUNICORN_THREADS=200  # 200 connections per replica
```

### Adjusting Resource Limits

Edit `docker-compose.yaml` deploy section for each service:

```yaml
deploy:
  resources:
    limits:
      cpus: '32'      # Maximum CPU cores
      memory: 64G     # Maximum memory
    reservations:
      cpus: '16'      # Minimum guaranteed cores
      memory: 32G     # Minimum guaranteed memory
```

## Monitoring

### Container Health

```bash
# Check health status of all containers
docker compose ps

# View specific service health
docker inspect $(docker compose ps -q zndraw) --format='{{.State.Health.Status}}'
```

### Resource Usage

```bash
# Real-time resource monitoring
docker stats

# Per-service resource usage
docker stats $(docker compose ps -q zndraw)
docker stats $(docker compose ps -q celery-worker)
```

### Application Logs

```bash
# Nginx access logs (includes response times)
docker compose logs nginx | grep "HTTP"

# Application errors
docker compose logs zndraw | grep "ERROR"

# Celery task processing
docker compose logs celery-worker | grep "Task"

# Redis operations
docker compose logs redis
```

## Performance Tuning

### For Upload-Heavy Workloads

Increase Celery worker replicas and concurrency:

```yaml
celery-worker:
  command: ["celery", "-A", "zndraw.app.make_celery", "worker", "--loglevel=info", "--concurrency=16"]
  deploy:
    replicas: 4  # 64 parallel upload processing tasks
```

### For WebSocket-Heavy Workloads

Increase ZnDraw app replicas and threads:

```yaml
zndraw:
  environment:
    - GUNICORN_THREADS=200  # More concurrent WebSocket connections
  deploy:
    replicas: 5  # 1000 total concurrent connections
```

### For Memory-Constrained Environments

Reduce threads per container and use more replicas:

```yaml
zndraw:
  environment:
    - GUNICORN_THREADS=50  # Fewer threads = less memory per container
  deploy:
    replicas: 6  # Same total capacity (300 connections) but spread across more containers
    resources:
      limits:
        memory: 16G  # Less memory per container
```

## Troubleshooting

### Nginx Can't Connect to Backend

```bash
# Check if zndraw replicas are running
docker compose ps zndraw

# Check if zndraw is healthy
docker compose ps | grep zndraw | grep "healthy"

# View nginx error logs
docker compose logs nginx | grep "upstream"
```

### WebSocket Disconnects

WebSocket issues usually indicate sticky session problems:

```bash
# Verify ip_hash is configured in nginx.conf
grep "ip_hash" nginx.conf

# Check if client IP is preserved
docker compose logs nginx | grep "X-Real-IP"
```

### High Memory Usage

```bash
# Identify which service is using memory
docker stats --no-stream

# Reduce threads or replicas for that service
# Edit docker-compose.yaml and restart:
docker compose up -d
```

### Slow Upload Processing

```bash
# Check Celery worker queue size
docker compose exec -it $(docker compose ps -q celery-worker | head -1) \
  celery -A zndraw.app.make_celery inspect active

# Increase worker replicas if queue is backed up
docker compose up -d --scale celery-worker=4
```

## Production Checklist

Before deploying to production:

- [ ] Change `FLASK_SECRET_KEY` in docker-compose.yaml
- [ ] Set `ZNDRAW_ADMIN_USERNAME` and `ZNDRAW_ADMIN_PASSWORD` if needed
- [ ] Configure SSL/TLS termination (use nginx with certbot or external load balancer)
- [ ] Set up log aggregation (ELK stack, Grafana, etc.)
- [ ] Configure monitoring and alerting (Prometheus + Grafana)
- [ ] Set up automated backups for `./data` and `redis-data` volumes
- [ ] Review and adjust resource limits based on server capacity
- [ ] Test failover by stopping individual containers
- [ ] Load test the application under expected peak load
- [ ] Document your specific scaling configuration

## Backup and Recovery

### Backup Data

```bash
# Stop services to ensure data consistency
docker compose stop

# Backup data directory
tar -czf backup-$(date +%Y%m%d).tar.gz ./data ./uploads

# Backup Redis data
docker compose exec redis redis-cli BGSAVE
docker cp $(docker compose ps -q redis):/data ./redis-backup

# Restart services
docker compose start
```

### Restore Data

```bash
# Stop services
docker compose stop

# Restore data
tar -xzf backup-YYYYMMDD.tar.gz

# Restore Redis
docker compose start redis
docker cp ./redis-backup $(docker compose ps -q redis):/data
docker compose restart redis

# Start all services
docker compose start
```

## Advanced: Manual Scaling

You can override the replicas defined in docker-compose.yaml:

```bash
# Start with custom replica counts
docker compose up -d \
  --scale zndraw=5 \
  --scale celery-worker=3

# Scale up during high load
docker compose up -d --scale zndraw=8 --no-recreate

# Scale down during low load
docker compose up -d --scale zndraw=2 --no-recreate
```

## SSL/TLS Configuration

For HTTPS, you have two options:

### Option 1: External Load Balancer (Recommended)

Use a cloud load balancer (AWS ALB, GCP Load Balancer, etc.) that handles SSL termination and forwards HTTP traffic to nginx.

### Option 2: Nginx SSL Termination

1. Obtain SSL certificates (Let's Encrypt recommended)
2. Mount certificates into nginx container
3. Update nginx.conf to listen on port 443

```yaml
nginx:
  ports:
    - "80:80"
    - "443:443"
  volumes:
    - ./ssl/cert.pem:/etc/nginx/cert.pem:ro
    - ./ssl/key.pem:/etc/nginx/key.pem:ro
```

Add to nginx.conf:
```nginx
server {
    listen 443 ssl http2;
    ssl_certificate /etc/nginx/cert.pem;
    ssl_certificate_key /etc/nginx/key.pem;
    # ... rest of config
}
```
