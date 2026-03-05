# Technology Stack

**Analysis Date:** 2026-03-05

## Languages

**Primary:**
- Python 3.11+ - Backend API, CLI, extensions, Python client (`src/zndraw/`)
- TypeScript 5.9+ - Frontend SPA (`frontend/src/`)

**Secondary:**
- JavaScript (ESNext) - One remark plugin (`frontend/src/utils/remark-frame-link.js`)

## Runtime

**Backend:**
- Python >=3.11 (`.python-version` = 3.11)
- Uvicorn ASGI server (`uvicorn>=0.40.0`)
- Production: Gunicorn + UvicornWorker (multi-process), see `src/zndraw/cli.py`

**Frontend:**
- Bun (preferred) or Node.js runtime
- Vite 7.x dev server with SWC-based React plugin (`@vitejs/plugin-react-swc`)

**Package Managers:**
- Backend: `uv` (lockfile: `uv.lock` present)
- Frontend: `bun` (lockfile: `bun.lock` present, config: `frontend/bunfig.toml`)

## Frameworks

**Core:**
- FastAPI >=0.128.0 - REST API framework (`src/zndraw/app.py`)
- python-socketio >=5.16.0 - Real-time WebSocket layer (`src/zndraw/socketio.py`)
- React 19.2 - Frontend UI (`frontend/src/frontend.tsx`)
- React Router 7.x - Client-side routing (`frontend/src/App.tsx`)

**State Management:**
- Zustand 5.x - Frontend global state (`frontend/src/store.tsx`, `frontend/src/stores/`)
- TanStack React Query 5.x - Server state / caching (`frontend/src/App.tsx`)
- Immer 10.x - Immutable state updates in Zustand slices

**3D Rendering:**
- Three.js 0.180 - WebGL 3D engine (`frontend/src/components/three/`)
- @react-three/fiber 9.x - React Three.js bindings
- @react-three/drei 10.x - Three.js helpers (controls, text, etc.)
- @react-three/gpu-pathtracer 0.3 - Path-traced rendering

**UI Components:**
- MUI (Material UI) 7.x - Component library (`@mui/material`, `@mui/icons-material`)
- MUI X Data Grid 8.x - Table components
- MUI X Tree View 8.x - Tree components
- JSON Forms 3.6 - Dynamic form generation from JSON Schema (`@jsonforms/react`, `@jsonforms/material-renderers`)

**Data Visualization:**
- Plotly.js 3.x - Interactive charts (`frontend/src/components/FigureWindow.tsx`)
- Recharts 3.x - React charting library
- Konva 10.x + react-konva 19.x - 2D canvas graphics

**Database ORM:**
- SQLModel >=0.0.31 - Pydantic-integrated ORM (`src/zndraw/models.py`)
- SQLAlchemy[asyncio] >=2.0.46 - Async database engine

**Task Queue:**
- TaskIQ with `taskiq-redis` (ListQueueBroker) - Background job execution (`src/zndraw/broker.py`)
- zndraw-joblib - Job management library (custom package)

**Testing:**
- pytest >=9.0.2 with pytest-asyncio >=1.3.0 - Backend tests
- Playwright >=1.58.0 - Frontend E2E tests (`frontend/e2e/`)
- pytest-cov >=7.0.0 - Coverage reporting

**Build/Dev:**
- Hatchling + hatch-vcs - Python build system (`pyproject.toml`)
- Custom build hook: `hatch_build.py` - Builds frontend during `uv build`
- Vite 7.x - Frontend bundler (`frontend/vite.config.ts`)
- Ruff >=0.14.14 - Python linting and formatting
- Pyright >=1.1.408 - Python type checking
- Biome - Frontend lint/format (referenced in `package.json` scripts, no config file yet)
- Prettier 3.6 - Frontend formatting (dev dependency)

## Key Dependencies

**Critical (Backend):**
- `ase` >=3.27.0 - Atomic Simulation Environment (atomic structure I/O)
- `asebytes[h5md,lmdb,mongodb]` >=0.3.0a3 - Binary frame storage backends
- `redis` >=7.1.0 - Async Redis client (ephemeral state, pub/sub, locks)
- `fakeredis` >=2.33.0 - In-process Redis replacement for development
- `zndraw-auth` >=0.2.3 - Authentication package (fastapi-users based)
- `zndraw-socketio` >=0.1.5 - Socket.IO wrapper with DI support
- `zndraw-joblib` >=0.1.4 - Job/task system with provider support
- `pydantic-settings` >=2.0.0 - Environment-based configuration
- `pyjwt` >=2.10.1 - JWT token handling
- `msgpack` + `msgpack-numpy` - Binary serialization for frame data

**Critical (Frontend):**
- `axios` 1.x - HTTP client (`frontend/src/myapi/client.ts`)
- `socket.io-client` 4.x - WebSocket client (`frontend/src/socket.ts`)
- `@msgpack/msgpack` 3.x - Binary message decoding
- `three` 0.180 - 3D rendering engine
- `zustand` 5.x - State management
- `@tanstack/react-query` 5.x - Server state caching

**Scientific (Backend):**
- `vesin` >=0.4.2 - Neighbor list computation
- `splines` >=0.3.3 - Spline interpolation for curves
- `networkx` >=3.6.1 - Graph algorithms (connectivity)
- `scipy` >=1.17.0 - Scientific computing
- `pandas` >=3.0.0 - Data analysis
- `plotly` >=6.5.2 - Server-side figure generation

**Molecule/Chemistry (Frontend):**
- `ketcher-react` + `ketcher-core` + `ketcher-standalone` 3.9 - SMILES molecule editor (lazy-loaded)

**Infrastructure (Backend):**
- `uvicorn` >=0.40.0 - ASGI server
- `httpx` >=0.28.0 - Async HTTP client (internal executor, tests)
- `typer` >=0.21.1 - CLI framework (`src/zndraw/cli.py`)
- `argon2-cffi` >=25.1.0 - Password hashing
- `asyncpg` >=0.31.0 - PostgreSQL async driver (production)
- `aiosqlite` (via sqlalchemy) - SQLite async driver (development)
- `lmdb` >=1.5.0 - LMDB storage backend

## Configuration

**Backend Environment (ZNDRAW_ prefix):**
- All settings via `pydantic_settings.BaseSettings` in `src/zndraw/config.py`
- `ZNDRAW_REDIS_URL` - Redis connection (None = auto-start fakeredis)
- `ZNDRAW_DATABASE_URL` - SQL database (default: `sqlite+aiosqlite://`)
- `ZNDRAW_STORAGE` - Frame storage URI (`memory://`, `*.lmdb`, `mongodb://host/db`)
- `ZNDRAW_HOST` / `ZNDRAW_PORT` - Server bind address (default: `0.0.0.0:8000`)
- `ZNDRAW_GUEST_PASSWORD` / `ZNDRAW_WORKER_PASSWORD` - Auth passwords
- `ZNDRAW_WORKER_ENABLED` - In-process TaskIQ worker toggle
- `ZNDRAW_INIT_DB_ON_STARTUP` - Auto-create tables (False for multi-worker)
- `ZNDRAW_SIMGEN_ENABLED` - Feature flag for SiMGen integration
- Auth settings: `ZNDRAW_AUTH_` prefix (via `zndraw-auth`)

**Frontend Build:**
- `VITE_APP_VERSION` - Injected at build time from Python package version
- Dev proxy: `/v1` and `/api` to `http://localhost:8000`, `/socket.io` to `ws://localhost:8000`

**Build Configuration Files:**
- `pyproject.toml` - Python project config, ruff settings, build system
- `frontend/vite.config.ts` - Vite build config with manual chunks (three, plotly, mui, vendor)
- `frontend/tsconfig.json` - TypeScript config (ESNext target, strict, `@/*` path alias)
- `frontend/package.json` - Frontend dependencies and scripts
- `hatch_build.py` - Custom hook to build frontend during Python package build

## Platform Requirements

**Development:**
- Python 3.11+
- Bun (preferred) or Node.js for frontend
- Redis (optional - fakeredis used by default)
- `uv` package manager

**Production (Docker):**
- Base image: `ghcr.io/astral-sh/uv:python3.12-bookworm-slim`
- Redis 7 (required, via Docker service)
- PostgreSQL 17 (recommended, async via `asyncpg`)
- MongoDB 7 (recommended for frame storage)
- Caddy 2 (reverse proxy, load balancer)
- Multi-replica deployment: N zndraw instances + dedicated TaskIQ workers

**CLI Entry Points:**
- `zndraw` - Main CLI (`src/zndraw/cli.py:app`)
- `zndraw-cli` - Agent CLI (`src/zndraw/cli_agent/`)
- `zndraw-db` - Database init CLI (`src/zndraw/cli.py:db_app`)

---

*Stack analysis: 2026-03-05*
