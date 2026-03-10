# External Integrations

**Analysis Date:** 2026-03-05

## APIs & External Services

**Internal REST API (Backend -> Frontend):**
- FastAPI app at `src/zndraw/app.py` serving all `/v1/*` endpoints
- Frontend API client: `frontend/src/myapi/client.ts` (axios-based)
- API prefix: `/v1/` for all REST endpoints
- Error format: RFC 9457 Problem JSON (`application/problem+json`)

**API Route Modules (all at `src/zndraw/routes/`):**
- `auth.py` - `/v1/auth/*` - Authentication (guest, JWT login, register, user management)
- `rooms.py` - `/v1/rooms/*` - Room CRUD, listing, frame selection
- `frames.py` - `/v1/rooms/{room_id}/frames` - Frame read/write (msgpack binary)
- `trajectory.py` - `/v1/rooms/{room_id}/trajectory` - File upload/download
- `geometries.py` - `/v1/rooms/{room_id}/geometries/*` - 3D geometry CRUD
- `bookmarks.py` - `/v1/rooms/{room_id}/bookmarks/*` - Frame bookmarks
- `chat.py` - `/v1/rooms/{room_id}/chat/*` - Chat messages
- `figures.py` - `/v1/rooms/{room_id}/figures/*` - Plotly figures
- `selection_groups.py` - `/v1/rooms/{room_id}/selection-groups/*` - Named selections
- `presets.py` - `/v1/rooms/{room_id}/presets/*` - Visual presets
- `screenshots.py` - `/v1/rooms/{room_id}/screenshots/*` - Screenshot capture
- `edit_lock.py` - `/v1/rooms/{room_id}/edit-lock` - Distributed edit locking
- `step.py` - `/v1/rooms/{room_id}/step` - Current frame step
- `progress.py` - `/v1/rooms/{room_id}/progress/*` - Progress trackers
- `tools.py` - `/v1/tools/*` - Utility tools (RDKit image conversion)
- `admin.py` - `/v1/admin/*` - Admin operations (shutdown, user management)
- `server_settings.py` - `/v1/server-settings/*` - Server-wide config (default room)
- `problems.py` - Error response utilities
- `utility.py` - `/v1/version`, `/v1/config/*` - Version and global settings
- zndraw-joblib routes: `/v1/joblib/rooms/{room_id}/jobs/*`, `/v1/joblib/rooms/{room_id}/tasks/*`, `/v1/joblib/rooms/{room_id}/providers/*`

## Real-Time Communication

**Socket.IO (WebSocket):**
- Server: `src/zndraw/socketio.py` - `socketio.AsyncServer(async_mode="asgi")`
- Client (frontend): `frontend/src/socket.ts` - `socket.io-client`
- Client (Python): `src/zndraw/client.py` - `python-socketio` sync client
- Transport: WebSocket with ASGI mount at `/socket.io`
- Authentication: JWT token passed via `auth.token` on connect
- Multi-process pub/sub: `socketio.AsyncRedisManager` for cross-worker events

**Socket.IO Events (defined in `src/zndraw/socket_events.py`):**

*Client -> Server (request):*
- `room_join` - Join a room (returns session_id, frame_count, step, camera_key)
- `room_leave` - Leave current room
- `user_get` - Get authenticated user info
- `typing_start` / `typing_stop` - Typing indicators

*Server -> Client (broadcast):*
- `session_joined` / `session_left` - User presence notifications
- `geometry_invalidate` - Geometry created/updated/deleted (triggers refetch)
- `frames_invalidate` - Frame data changed (triggers refetch)
- `lock_update` - Edit lock acquired/released/refreshed
- `typing` - Typing indicator broadcast
- `progress_start` / `progress_update` / `progress_complete` - Task progress
- `screenshot_request` - Server asks frontend to capture screenshot
- `message_new` / `message_edited` - Chat message events

**Socket.IO DI Pattern:**
- `tsio.app = app` set in lifespan enables `Depends()` in socket handlers
- Same dependency types as FastAPI: `RedisDep`, `SessionDep`, `StorageDep`
- See `src/zndraw/socketio.py` for handler registration

## Data Storage

**SQL Databases:**
- Development: SQLite via `aiosqlite` (in-memory default: `sqlite+aiosqlite://`)
- Production: PostgreSQL 17 via `asyncpg` (`postgresql+asyncpg://...`)
- ORM: SQLModel (`src/zndraw/models.py`)
- Connection: `ZNDRAW_DATABASE_URL` env var
- Engine creation: `zndraw_auth.db.create_engine_for_url()`
- Session factory: `async_sessionmaker` stored on `app.state.session_maker`
- SQLite locking: Automatic `asyncio.Lock` wrapper applied in `src/zndraw/database.py`

**SQL Models (`src/zndraw/models.py`):**
- `Room` - Room metadata (id, description, locked, step, default_camera)
- `Message` - Chat messages
- `RoomMembership` - Room access control with roles (member/moderator/owner)
- `RoomGeometry` - Persistent geometry instances (type + JSON config)
- `RoomBookmark` - Frame bookmarks
- `SelectionGroup` - Named atom selection groups
- `RoomFigure` - Plotly figure data
- `Screenshot` - Screenshot metadata (files stored on disk)
- `RoomPreset` - Visual presets with rules
- `ServerSettings` - Singleton server-wide config
- From zndraw-joblib: `Job`, `Task`, `Worker`, `WorkerJobLink`
- From zndraw-auth: `User` (via fastapi-users)

**Frame Storage (Binary Data):**
- Abstraction: `AsebytesStorage` (`src/zndraw/storage/asebytes_backend.py`)
- Router: `StorageRouter` (`src/zndraw/storage/router.py`) - wraps storage with provider mount support
- Backends (via `asebytes` library, selected by `ZNDRAW_STORAGE` URI):
  - `memory://` - In-memory (default, development)
  - `*.lmdb` - LMDB file on disk
  - `mongodb://host:port/db` - MongoDB (production)
- Format: `dict[bytes, bytes]` per frame (keys like `b"arrays.positions"`, values are msgpack-numpy)

**Redis (Ephemeral State):**
- Client: `redis.asyncio` stored on `app.state.redis` (`decode_responses=True`)
- Development: `fakeredis.TcpFakeServer` auto-started when `ZNDRAW_REDIS_URL` is None
- Production: Redis 7 (via Docker)
- Key patterns defined in `src/zndraw/redis.py` (`RedisKey` class):
  - `room:{room_id}:cameras` - Hash of session cameras (presence tracking)
  - `room:{room_id}:active-cameras` - Hash: session_id -> camera_key
  - `room:{room_id}:edit-lock` - JSON value with TTL (10s default)
  - `room:{room_id}:progress` - Hash of active progress trackers
  - `room:{room_id}:provider_frame_count` - Provider-backed frame count
  - `download-token:{token}` - Temporary download tokens
  - `provider-result:*` / `provider-inflight:*` - Provider caching

**File Storage (Screenshots):**
- Location: `ZNDRAW_MEDIA_PATH` (default: `zndraw-media/`)
- Screenshots stored as PNG files on disk
- Metadata in SQL `Screenshot` table

**Caching:**
- Provider result cache: `CompositeResultBackend` (`src/zndraw/result_backends.py`)
  - Frame data -> `StorageResultBackend` (same storage as frame storage)
  - Other data -> `RedisResultBackend` (Redis with TTL)
  - Pub/sub notifications for long-polling via Redis channels

## Authentication & Identity

**Auth Provider:**
- `zndraw-auth` package (wrapper around `fastapi-users`)
- JWT-based authentication
- Auth routes: `src/zndraw/routes/auth.py` at `/v1/auth/*`

**Authentication Flows:**
- Guest login: `POST /v1/auth/guest` - Creates anonymous user, returns JWT
- Password login: `POST /v1/auth/jwt/login` - OAuth2 form data
- Registration: `POST /v1/auth/register` - Creates user (auto-login after)
- User info: `GET /v1/auth/users/me` - Current user details
- CLI login: `/v1/auth/cli-login/*` - Browser-approved CLI authentication

**Frontend Auth (`frontend/src/utils/auth.ts`):**
- Token stored in `localStorage` as `zndraw_jwt_token`
- `acquireToken()` - Coalesced token acquisition (validates existing or creates guest)
- 401 interceptor on axios: auto-retry with fresh guest token
- Socket.IO auth: JWT passed via `socket.auth = { token }` on connect

**Internal Users:**
- Default admin: Created via `zndraw-auth` settings (`ZNDRAW_AUTH_DEFAULT_ADMIN_*`)
- Internal worker: `worker@internal.user` superuser (created in `src/zndraw/database.py`)

**Dependencies (`src/zndraw/dependencies.py`):**
- `CurrentUserDep` - Requires authenticated user
- `AdminUserDep` - Requires superuser
- `OptionalUserDep` - Optional user (for public endpoints)
- `WritableRoomDep` - Checks room existence + admin lock + edit lock
- `WritableGeometryDep` - Checks room + lock + geometry ownership

## Task Queue / Background Jobs

**TaskIQ + Redis:**
- Broker: `ListQueueBroker` from `taskiq-redis` (`src/zndraw/broker.py`)
- In-process worker: Started during lifespan when `ZNDRAW_WORKER_ENABLED=true`
- External workers: `taskiq worker zndraw.broker:broker`
- Executor: `InternalExtensionExecutor` (`src/zndraw/executor.py`)
  - Connects back to FastAPI as internal worker user
  - Runs extensions in `asyncio.to_thread()` (sync extensions in thread pool)

**Built-in Extensions (`src/zndraw/extensions/`):**
- Modifiers: Center, ChangeType, Delete, Duplicate, Empty, Replicate, Wrap, etc.
- Selections: All, ConnectedParticles, IdenticalSpecies, Invert, Neighbour, etc.
- Analysis: Distance, DihedralAngle, Properties1D, Properties2D
- Molecule Building: AddFromSMILES, PackBox (requires `molify`, `rdkit`)
- Filesystem: File browser provider

**Job Lifecycle:**
- Frontend submits task via `POST /v1/joblib/rooms/{room_id}/tasks/{job_name}`
- TaskIQ dispatches to worker (in-process or external)
- Worker creates `ZnDraw` client, runs extension, reports status
- Progress broadcast via Socket.IO events
- Sweeper task cleans up stale workers (`run_sweeper` in `src/zndraw/database.py`)

## Monitoring & Observability

**Error Tracking:**
- No external error tracking service (Sentry, etc.)
- RFC 9457 Problem JSON for all API errors (`src/zndraw/exceptions.py`)

**Logs:**
- Python `logging` module throughout backend
- Uvicorn access/error logs
- Docker healthcheck on port 8000

**Health Check:**
- Docker: `python -c "import urllib.request; urllib.request.urlopen('http://localhost:8000').read()"`
- Interval: 30s, timeout: 10s, start_period: 40-60s

## CI/CD & Deployment

**Docker:**
- `Dockerfile` - Multi-stage build (builder + runtime)
- Builder: `ghcr.io/astral-sh/uv:python3.12-bookworm-slim` + bun
- Frontend built during wheel build via `hatch_build.py`
- Static assets embedded in Python package at `src/zndraw/static/`

**Docker Compose Profiles (`docker/`):**
- `docker/standalone/docker-compose.yaml` - Single-node (zndraw + redis)
- `docker/production/docker-compose.yaml` - Multi-replica with:
  - Redis 7 (data persistence)
  - MongoDB 7 (frame storage)
  - PostgreSQL 17 (SQL database)
  - Caddy 2 (reverse proxy / load balancer)
  - N zndraw replicas (horizontally scaled)
  - Dedicated TaskIQ worker container
  - Template loader init container

**Image Registry:**
- `pythonf/zndraw:latest` (Docker Hub)

## Frontend-Backend Communication Summary

**HTTP (axios via `frontend/src/myapi/client.ts`):**
- Base URL: Same origin (Vite proxy in dev, same server in prod)
- Auth: `Authorization: Bearer {jwt}` header (auto-injected via interceptor)
- Session tracking: `X-Session-ID` header
- Lock tracking: `Lock-Token` header for mutations
- Binary data: `application/msgpack` for frame data (request + response)
- File uploads: `multipart/form-data`

**WebSocket (socket.io-client via `frontend/src/socket.ts`):**
- Connection: `io(undefined, { autoConnect: false })` - same origin
- Auth: `socket.auth = { token }` before connect
- Room lifecycle: `room_join` -> events -> `room_leave`
- Real-time events drive Zustand store updates (`frontend/src/hooks/useSocketManager.ts`)

**Frontend State Flow:**
1. Page load -> `acquireToken()` -> JWT from localStorage or guest login
2. Socket connect with JWT -> `room_join` event -> receive session_id, frame_count
3. REST API calls for data (frames, geometries, chat) with JWT + session headers
4. Socket.IO events trigger cache invalidation -> React Query refetches
5. Zustand store holds UI state (selections, playback, lock status)

## Environment Configuration

**Required env vars (production):**
- `ZNDRAW_REDIS_URL` - Redis connection string
- `ZNDRAW_DATABASE_URL` - PostgreSQL connection string
- `ZNDRAW_STORAGE` - Frame storage URI (mongodb:// for production)
- `ZNDRAW_AUTH_SECRET_KEY` - JWT signing secret (via zndraw-auth)
- `ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL` / `ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD`

**Optional env vars:**
- `ZNDRAW_HOST` / `ZNDRAW_PORT` - Server bind (default: 0.0.0.0:8000)
- `ZNDRAW_GUEST_PASSWORD` / `ZNDRAW_WORKER_PASSWORD` - Auth passwords
- `ZNDRAW_WORKER_ENABLED` - In-process worker (default: true, false in Docker)
- `ZNDRAW_SERVER_URL` - For external TaskIQ workers to reach FastAPI
- `ZNDRAW_INIT_DB_ON_STARTUP` - Auto-create tables (default: true)
- `ZNDRAW_SIMGEN_ENABLED` - SiMGen feature flag
- `ZNDRAW_EDIT_LOCK_TTL` - Edit lock TTL in seconds (default: 10)
- `ZNDRAW_MEDIA_PATH` - Screenshot storage path

**Secrets location:**
- `.env` file (referenced in docker-compose, gitignored)
- Never committed to repository

## Webhooks & Callbacks

**Incoming:**
- None (no external webhook endpoints)

**Outgoing:**
- None (no outbound webhook calls)

**Internal Callbacks:**
- TaskIQ executor connects back to FastAPI server via HTTP as internal worker user
- Socket.IO pub/sub via Redis for cross-worker event broadcasting
- Provider system: workers register via Socket.IO, receive dispatch via events

---

*Integration audit: 2026-03-05*
