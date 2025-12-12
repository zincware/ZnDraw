# ZnDraw Development Docker Setup

This directory contains a development-optimized Docker Compose configuration with **hot reload** support for both frontend and backend development.

## Quick Start

```bash
# Start development environment
docker compose -f docker-compose.dev.yaml up -d

# Watch logs
docker compose -f docker-compose.dev.yaml logs -f

# Stop everything
docker compose -f docker-compose.dev.yaml down
```

## Architecture

- **Frontend Dev Server** (Vite) - Port 5173 with hot reload
- **Backend** (ZnDraw Flask app) - Port 5001 (also via 5000 through nginx)
- **Nginx** - Port 5000 (proxies to frontend & backend)
- **Redis** - Port 6380
- **Celery Worker** - Background tasks

## Development Workflow

### Frontend Changes (Hot Reload - No Restart Needed!)

1. Edit files in `../app/src/`
2. Vite automatically detects changes and hot-reloads in the browser
3. See changes instantly at http://localhost:5000

**Example:**
```bash
# Edit a React component
vim ../app/src/components/MyComponent.tsx

# Changes appear immediately in browser!
```

### Backend Python Changes (Restart Required)

1. Edit files in `../src/zndraw/`
2. Restart the backend container

```bash
# After editing Python code
docker compose -f docker-compose.dev.yaml restart zndraw

# Watch logs to see restart
docker compose -f docker-compose.dev.yaml logs -f zndraw
```

**Example:**
```bash
# Edit a Python file
vim ../src/zndraw/routes.py

# Restart backend
docker compose -f docker-compose.dev.yaml restart zndraw
```

### Dependencies Changes

**Frontend (package.json):**
```bash
# Rebuild frontend container
docker compose -f docker-compose.dev.yaml build frontend-dev
docker compose -f docker-compose.dev.yaml up -d frontend-dev
```

**Backend (pyproject.toml):**
```bash
# Rebuild backend containers
docker compose -f docker-compose.dev.yaml build zndraw celery-worker
docker compose -f docker-compose.dev.yaml up -d
```

## Access Points

- **Main Application**: http://localhost:5000 (via Nginx proxy)
- **Backend Direct**: http://localhost:5001 (bypass Nginx for debugging)
- **Frontend Direct**: http://localhost:5173 (Vite dev server with HMR)
- **Redis**: localhost:6380

## Admin Access

- **Username**: `admin`
- **Password**: `admin123`

Admin mode is enabled in development for testing admin features.

## Configuration

### Environment Variables

The development setup uses these key settings:

- `FLASK_DEBUG=1` - Flask debug mode enabled
- `ZNDRAW_LOG_LEVEL=DEBUG` - Verbose logging
- `NODE_ENV=development` - Frontend development mode
- `ZNDRAW_ADMIN_USERNAME=admin` - Admin username
- `ZNDRAW_ADMIN_PASSWORD=admin123` - Admin password

### Volume Mounts

**Frontend Source:**
- `../app/` → `/build/` (full frontend source for hot reload)

**Backend Source:**
- `../src/` → `/app/src/` (Python source code)
- `../pyproject.toml` → `/app/pyproject.toml/` (dependency manifest)

**Data:**
- `zndraw-data` - Persistent database storage
- `upload-temp` - Temporary upload storage
- `./template-files/` - Template files (s22.xyz)

## Comparison: Dev vs Production

| Feature | Development | Production |
|---------|------------|------------|
| **Frontend** | Vite dev server (hot reload) | Pre-built static files |
| **Backend Replicas** | 1 | 3 |
| **Celery Workers** | 1 | 2 |
| **Log Level** | DEBUG | INFO |
| **Flask Debug** | Enabled | Disabled |
| **Source Mounts** | Yes (hot reload) | No (baked in image) |
| **Rebuild on Change** | No (just restart) | Yes (full rebuild) |

## Debugging

### View Logs

```bash
# All services
docker compose -f docker-compose.dev.yaml logs -f

# Specific service
docker compose -f docker-compose.dev.yaml logs -f zndraw
docker compose -f docker-compose.dev.yaml logs -f frontend-dev
docker compose -f docker-compose.dev.yaml logs -f celery-worker
```

### Check Service Health

```bash
docker compose -f docker-compose.dev.yaml ps
```

### Access Container Shell

```bash
# Backend container
docker exec -it zndraw-backend-dev bash

# Frontend container
docker exec -it zndraw-frontend-dev sh
```

### Inspect Network

```bash
# Check if services can reach each other
docker exec zndraw-backend-dev ping frontend-dev
docker exec zndraw-backend-dev curl http://frontend-dev:5173
```

## Troubleshooting

### Frontend not hot-reloading

1. Check Vite dev server is running:
   ```bash
   docker compose -f docker-compose.dev.yaml logs frontend-dev
   ```

2. Verify port 5173 is accessible:
   ```bash
   curl http://localhost:5173
   ```

3. Check browser console for WebSocket connection errors

### Backend changes not appearing

1. Did you restart the container?
   ```bash
   docker compose -f docker-compose.dev.yaml restart zndraw
   ```

2. Check if Python files are mounted correctly:
   ```bash
   docker exec zndraw-backend-dev ls -la /app/src/zndraw/
   ```

3. Clear Python cache:
   ```bash
   docker exec zndraw-backend-dev find /app -name "*.pyc" -delete
   docker exec zndraw-backend-dev find /app -name "__pycache__" -type d -exec rm -rf {} +
   docker compose -f docker-compose.dev.yaml restart zndraw
   ```

### Port conflicts

If ports 5000, 5001, 5173, or 6380 are in use:

```bash
# Find what's using the port
lsof -i :5000

# Stop production containers if running
cd ../docker
docker compose down
```

## Tips

1. **Use production config for final testing**
   ```bash
   docker compose -f docker-compose.yaml up -d
   ```

2. **Keep dev and prod separate**
   - Dev uses different container names (`-dev` suffix)
   - Dev uses separate volumes
   - Can't run both simultaneously (port conflicts)

3. **Commit often during development**
   - Hot reload makes it easy to test quickly
   - Git commit tested changes frequently

4. **Monitor logs during development**
   ```bash
   docker compose -f docker-compose.dev.yaml logs -f
   ```

5. **Clean up regularly**
   ```bash
   # Stop and remove everything
   docker compose -f docker-compose.dev.yaml down --volumes

   # Remove dev images
   docker compose -f docker-compose.dev.yaml down --rmi all
   ```

## Performance Notes

- **First startup**: Takes 1-2 minutes (building frontend, installing deps)
- **Subsequent startups**: ~10 seconds (using cached layers)
- **Frontend hot reload**: Instant (< 1 second)
- **Backend restart**: ~5-10 seconds
- **Full rebuild**: 2-3 minutes (only needed for dependency changes)

## Next Steps

1. Make changes to frontend or backend code
2. See changes reflected immediately (frontend) or after restart (backend)
3. Test thoroughly in dev environment
4. Test in production config before deploying
5. Deploy to production

Happy coding! 🚀
