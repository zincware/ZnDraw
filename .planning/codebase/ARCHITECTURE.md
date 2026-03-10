# Architecture

**Analysis Date:** 2026-03-05

## Pattern Overview

**Overall:** Full-stack real-time application with REST + Socket.IO dual communication.

**Key Characteristics:**
- FastAPI backend with SQLModel/SQLAlchemy for persistent state, Redis for ephemeral state
- React + Three.js frontend for 3D molecular visualization in the browser
- Socket.IO for real-time invalidation-based synchronization (server broadcasts events, clients refetch via REST)
- Pydantic models as the schema contract between backend and frontend (JSON Schema for geometry/extension forms)
- TaskIQ-based job system for running compute extensions (modifiers, analysis, selections)
- Multi-client architecture: browser frontend, Python `ZnDraw` client, and CLI agent all use the same REST + Socket.IO API

## Backend Layers

**Application Entry:**
- Purpose: FastAPI app assembly, router registration, Socket.IO mounting
- Location: `src/zndraw/app.py`
- Contains: Router includes, exception handler registration, static file serving, SPA catch-all
- Entry point: `zndraw.app:socket_app` (ASGI app wrapping FastAPI + Socket.IO)
- Dev runner: `main.py` runs uvicorn with reload

**Configuration:**
- Purpose: Environment-based settings with `ZNDRAW_` prefix
- Location: `src/zndraw/config.py`
- Contains: `Settings` (Pydantic BaseSettings), `SettingsDep` (FastAPI dependency)
- Key settings: `redis_url`, `storage`, `database_url`, `host`, `port`, `worker_enabled`, `simgen_enabled`

**Lifespan / Resource Management:**
- Purpose: Initializes and tears down all application resources in correct order
- Location: `src/zndraw/database.py`
- Contains: Database engine, session maker, Redis client, frame storage, Socket.IO manager, TaskIQ broker/worker, job sweeper
- Pattern: `@contextlib.asynccontextmanager` lifespan stored in `app.state.*`
- Auto-starts `fakeredis.TcpFakeServer` when no `REDIS_URL` is configured (dev mode)
- Creates internal worker user on startup for built-in extension execution

**REST API Routes:**
- Purpose: CRUD endpoints for all room resources
- Location: `src/zndraw/routes/` (one file per resource)
- Key files:
  - `src/zndraw/routes/rooms.py` - Room CRUD, listing, duplication
  - `src/zndraw/routes/frames.py` - Frame read/write with msgpack binary encoding
  - `src/zndraw/routes/geometries.py` - Geometry CRUD with Pydantic model validation
  - `src/zndraw/routes/chat.py` - Chat messages with pagination
  - `src/zndraw/routes/bookmarks.py` - Frame bookmarks
  - `src/zndraw/routes/figures.py` - Plotly figure management
  - `src/zndraw/routes/step.py` - Current frame step (with broadcast)
  - `src/zndraw/routes/edit_lock.py` - Distributed edit locks via Redis TTL
  - `src/zndraw/routes/screenshots.py` - Screenshot capture/upload
  - `src/zndraw/routes/presets.py` - Visual preset management
  - `src/zndraw/routes/selection_groups.py` - Named selection groups
  - `src/zndraw/routes/trajectory.py` - File upload/download (ASE formats)
  - `src/zndraw/routes/auth.py` - Guest login, JWT login, registration (delegates to zndraw-auth)
  - `src/zndraw/routes/admin.py` - User management, server shutdown
  - `src/zndraw/routes/server_settings.py` - Default room, global config
  - `src/zndraw/routes/tools.py` - RDKit molecule image generation
  - `src/zndraw/routes/utility.py` - Version endpoint
  - `src/zndraw/routes/problems.py` - RFC 9457 problem type documentation endpoint
  - `src/zndraw/routes/progress.py` - tqdm-style progress trackers via Redis
- API prefix: `/v1/` for all routes
- External package routes: `zndraw_joblib.router` mounted at `/v1/joblib/`

**Socket.IO Events:**
- Purpose: Real-time room management and event broadcasting
- Location: `src/zndraw/socketio.py`
- Contains: Connection lifecycle (connect/disconnect), room join/leave, typing indicators
- Pattern: `zndraw-socketio` wrapper with DI support (`Depends()` same as FastAPI)
- Event models: `src/zndraw/socket_events.py` (Pydantic models, snake_case event names derived from class names)
- Key events (client -> server): `room_join`, `room_leave`, `user_get`, `typing_start`, `typing_stop`
- Key broadcasts (server -> clients): `geometry_invalidate`, `frames_invalidate`, `frame_update`, `lock_update`, `session_joined`, `session_left`, `message_new`, `figure_invalidate`, `room_update`, `progress_start/update/complete`

**Dependencies (DI):**
- Purpose: FastAPI dependency injection for all shared resources
- Location: `src/zndraw/dependencies.py`
- Key dependencies:
  - `RedisDep` - async Redis client from `app.state.redis`
  - `StorageDep` - frame storage backend from `app.state.frame_storage`
  - `SioDep` - Socket.IO server wrapper from `app.state.tsio`
  - `CurrentUserDep` / `AdminUserDep` / `OptionalUserDep` - auth from zndraw-auth
  - `WritableRoomDep` - validates room exists + checks admin lock + edit lock
  - `WritableGeometryDep` - validates geometry write access (owner + lock checks)
  - `SessionMakerDep` / `CurrentUserFactoryDep` - scoped sessions for long-polling
  - `ResultBackendDep` - composite result backend for job system
- Override: `joblib_verify_writable_room` overridden with `get_writable_room_id` to enforce room locks

**Data Models (SQL):**
- Purpose: SQLModel tables for persistent room state
- Location: `src/zndraw/models.py`
- Tables: `Room`, `Message`, `RoomMembership`, `RoomGeometry`, `RoomBookmark`, `SelectionGroup`, `RoomFigure`, `Screenshot`, `RoomPreset`, `ServerSettings`
- External tables (imported at module level for `metadata.create_all`): `Job`, `Worker`, `Task`, `WorkerJobLink` from zndraw-joblib, `User` from zndraw-auth

**API Schemas:**
- Purpose: Pydantic models for request/response serialization
- Location: `src/zndraw/schemas.py`
- Pattern: Separate `*Request` and `*Response` models for each resource
- Generic envelopes: `CollectionResponse[T]`, `OffsetPage[T]`

**Error Handling:**
- Purpose: RFC 9457 Problem Details for all 4xx/5xx responses
- Location: `src/zndraw/exceptions.py`
- Pattern: `ProblemType` base class with `title`, `status` class vars, `.exception()` factory, `.raise_for_client()` for Python client
- Registry: `PROBLEM_TYPES` dict mapping kebab-case IDs to classes
- Problem URI: `/v1/problems/{kebab-case-name}` (documented at `src/zndraw/routes/problems.py`)

**Frame Storage:**
- Purpose: Room-scoped frame storage with pluggable backends
- Location: `src/zndraw/storage/`
- `AsebytesStorage` (`src/zndraw/storage/asebytes_backend.py`): wraps `asebytes.AsyncBlobIO` for memory, LMDB, MongoDB, or H5MD backends
- `StorageRouter` (`src/zndraw/storage/router.py`): delegates to `AsebytesStorage` + Redis for provider-backed frame counts, rejects writes to provider-mounted rooms
- Frame format: `dict[bytes, bytes]` (keys like `b"arrays.positions"`, values are msgpack-numpy encoded)

**Redis Ephemeral State:**
- Purpose: Session cameras, edit locks, progress trackers, provider caches
- Location: `src/zndraw/redis.py`
- Pattern: `RedisKey` class with static methods returning key patterns
- Key namespaces: `room:{id}:cameras`, `room:{id}:active-cameras`, `room:{id}:edit-lock`, `room:{id}:progress`, `provider-result:*`, `provider-inflight:*`

**Geometry System:**
- Purpose: Pydantic models for 3D scene objects (spheres, bonds, cameras, lights, etc.)
- Location: `src/zndraw/geometries/`
- Base: `BaseGeometry` (`src/zndraw/geometries/base.py`) with `owner`, `active`, `position`, `color`, `material`
- Registry: `geometries` dict in `src/zndraw/geometries/__init__.py` mapping type names to classes
- Types: Sphere, Arrow, Bond, Curve, CircleCurve, Cell, Floor, Camera, Box, Plane, Shape, DirectionalLight, AmbientLight, HemisphereLight, Fog, PathTracing, PropertyInspector
- JSON Schema generation: geometry models produce JSON Schema consumed by frontend JSONForms

**Extension System:**
- Purpose: Pluggable compute operations (modifiers, analysis, selections)
- Location: `src/zndraw/extensions/`
- Base: `Extension` (`src/zndraw/extensions/abc.py`) with `category` ClassVar and `run(vis)` method
- Categories: `modifiers`, `selections`, `analysis`
- Extension files: `src/zndraw/extensions/modifiers.py`, `src/zndraw/extensions/selections.py`, `src/zndraw/extensions/analysis.py`, `src/zndraw/extensions/molecule_building.py`, `src/zndraw/extensions/filesystem.py`
- Execution: `InternalExtensionExecutor` (`src/zndraw/executor.py`) runs extensions in threads via `asyncio.to_thread`, creates a `ZnDraw` client instance per task

**Job System (zndraw-joblib):**
- Purpose: Distributed task execution for extensions
- External package: `zndraw-joblib` (provides `Job`, `Worker`, `Task` models, `JobManager`, `ClaimedTask`, REST router)
- Broker: `taskiq_redis.ListQueueBroker` for task queuing
- Internal tasks: registered at startup from `_collect_extensions()`
- Worker: In-process TaskIQ receiver (disabled in Docker -- dedicated containers)
- Result backend: `CompositeResultBackend` routing frame data to storage, everything else to Redis

**Python Client:**
- Purpose: Synchronous client for programmatic access
- Location: `src/zndraw/client.py`
- Class: `ZnDraw` implements `MutableSequence[ase.Atoms]` for frame operations
- Accessor classes in `src/zndraw/accessors.py`: `Selections`, `Bookmarks`, `Geometries`, `Figures`, `ChatMessages`, `Extensions`, `Presets`, `Screenshots`, `Sessions`, `SelectionGroups`, `RoomMetadata`
- Uses `httpx` for REST, `python-socketio` for real-time events

**CLI:**
- Purpose: Command-line interface for server management and data upload
- Location: `src/zndraw/cli.py` (main CLI), `src/zndraw/cli_agent/` (agent CLI)
- Entry points: `zndraw` (main), `zndraw-cli` (agent), `zndraw-db` (database init)
- Agent CLI: structured JSON output for LLM agents, one sub-app per resource

## Frontend Architecture

**Entry Point:**
- Location: `frontend/src/frontend.tsx`
- Pattern: React 19 StrictMode, `createRoot`, HMR support via `import.meta.hot`

**App Shell:**
- Location: `frontend/src/App.tsx`
- Providers: `ThemeProvider` (MUI), `QueryClientProvider` (TanStack React Query), `RouterProvider` (react-router-dom)
- Theme: MUI `createTheme` with dark mode support (`colorSchemes: { dark: true }`)
- React Query: 30s stale time, 5min GC time

**Routing:**
- Framework: react-router-dom v7 (`createBrowserRouter`)
- Routes:
  - `/` -> `TemplateSelectionPage` (`frontend/src/pages/templateSelection.tsx`)
  - `/auth/cli` -> `CliLoginApprovePage` (`frontend/src/pages/cliLoginApprove.tsx`)
  - `/rooms` -> `RoomListPage` (`frontend/src/pages/roomList.tsx`)
  - `/rooms/:roomId/files` -> `FilesystemBrowserPage` (`frontend/src/pages/filesystemBrowser.tsx`)
  - `/rooms/:roomId` -> `MainPage` (`frontend/src/pages/landingPage.tsx`)
  - `/room/:roomId` -> `MainPage` (legacy alias)
- SPA fallback: backend serves `index.html` for unmatched routes (production)

**Main Page Layout (the room view):**
- Location: `frontend/src/pages/landingPage.tsx`
- Layout: Full-viewport flex column
  1. `AppBar` - toolbar with drawing/editing mode toggles, chat, screenshot, upload, room menu, user profile
  2. Content row: `SideBar` (left) + main area with `Canvas` (3D scene) + `WindowManager` (floating plot windows)
  3. `FrameProgressBar` - bottom playback slider
- Modals: `ChatWindow` (lazy loaded), `LoginDialog`, `RegisterDialog`, `AdminPanel`, `UserProfileDialog`, `ConnectionDialog`, `SiMGenTutorialDialog`
- Hooks: `useSocketManager` (real-time sync), `useKeyboardShortcuts`, `useDragAndDrop` (file upload)

**State Management:**
- Primary store: Zustand with slice pattern (`frontend/src/store.tsx`)
- Type: `AppState = ConnectionSlice & PlaybackSlice & SceneSlice & LockSlice & UISlice`
- Slices:
  - `ConnectionSlice` (`frontend/src/stores/slices/connectionSlice.ts`): roomId, user, sessionId, cameraKey, isConnected, serverVersion, globalSettings
  - `PlaybackSlice` (`frontend/src/stores/slices/playbackSlice.ts`): currentFrame, frameCount, playing, bookmarks, fps, frame_selection, synchronizedMode
  - `SceneSlice` (`frontend/src/stores/slices/sceneSlice.ts`): geometries, selections, selectionGroups, mode (view/drawing/editing), transformMode, activeCurveForDrawing, hoveredGeometryInstance, pendingFrameEdits
  - `LockSlice` (`frontend/src/stores/slices/lockSlice.ts`): superuserLock, userLock, lockToken, lock renewal timers
  - `UISlice` (`frontend/src/stores/slices/uiSlice.ts`): chatOpen, snackbar, progressTrackers, screenshotCapture
- Secondary stores:
  - `useGeometryStore` (`frontend/src/stores/geometryStore.ts`): geometry panel mode/selection (uses immer middleware)
  - `useWindowManagerStore` (`frontend/src/stores/windowManagerStore.ts`): floating plot window positions/sizes
  - `useFormStore` (`frontend/src/formStore.ts`): extension panel category/extension selection (uses immer middleware)
  - `useRoomsStore` (`frontend/src/roomsStore.tsx`): room list for overview page
- Selectors: `selectIsRoomReadOnly`, `selectActiveSelectionGroup` in `store.tsx`

**Backend Communication:**

1. **REST API Client:**
   - Location: `frontend/src/myapi/client.ts`
   - HTTP client: axios with request interceptor (JWT Bearer token, session ID, lock token)
   - Response interceptor: 401 -> acquire fresh guest token and retry once
   - Functions: one exported function per API endpoint (e.g., `listGeometries`, `getFrames`, `createRoom`, `submitTask`)
   - Frame data: msgpack binary encoding via `frontend/src/utils/msgpack-numpy.ts`

2. **Socket.IO Client:**
   - Location: `frontend/src/socket.ts`
   - Singleton: `socket = io(undefined, { autoConnect: false })`
   - Auth: `connectWithAuth()` acquires JWT, sets `socket.auth = { token }`, connects

3. **Real-time Sync (useSocketManager):**
   - Location: `frontend/src/hooks/useSocketManager.ts`
   - Pattern: Invalidation-based. Server broadcasts `*_invalidate` events, hook refetches via REST API and updates Zustand store
   - Room lifecycle: `connect` -> fetch version -> `room_join` -> fetch geometries (blocking) -> fetch secondary data (non-blocking)
   - Handles 404 on join by auto-creating room via REST, then retrying join
   - Event handlers for: frame_update, geometry_invalidate, frames_invalidate, lock_update, figure_invalidate, message_new, typing, progress_start/update/complete, room_update, selection_invalidate, bookmarks_invalidate

4. **React Query:**
   - Used for: frame data caching (`['frame', roomId, frameIndex, key]`), figure data, chat messages (infinite query), job schemas, task lists, default camera, metadata
   - Invalidation driven by Socket.IO events (targeted by query key predicates)
   - Custom hooks: `useFrameBatch`, `useFramePrefetcher`, `useGeometries`, `useJobs`, `useChat`, `useFigures`

**3D Rendering:**
- Framework: @react-three/fiber (React Three Fiber) + @react-three/drei (helpers)
- Canvas component: `frontend/src/components/Canvas.tsx`
- Geometry components in `frontend/src/components/three/`:
  - `Particles.tsx` (Sphere/atom rendering with instanced meshes)
  - `Bonds.tsx`, `Arrow.tsx`, `Box.tsx`, `Plane.tsx`, `Shape.tsx` (instanced geometry renderers)
  - `Camera.tsx` (camera viewports), `Lights.tsx` (DirectionalLight, AmbientLight, HemisphereLight)
  - `Curve.tsx`, `CircleCurve.tsx` (spline curves for drawing mode)
  - `Cell.tsx` (unit cell wireframe), `Floor.tsx`, `Fog.tsx`
  - `DrawingIndicator.tsx`, `EditingIndicator.tsx` (mode indicators)
  - `MultiGeometryTransformControls.tsx` (editing mode transform handles)
  - `HoverInfoBox.tsx`, `StaticInfoBox.tsx` (atom/geometry info overlays)
  - `ScreenshotProvider.tsx`, `PathtracingScreenshotCapture.tsx` (capture utilities)
  - `VirtualCanvas.tsx` (offscreen rendering for screenshots)
  - `materials.tsx` (Three.js material definitions)
- Path tracing: `frontend/src/components/PathTracingRenderer.tsx` using @react-three/gpu-pathtracer
- Camera controls: OrbitControls + custom `CameraManager.js` + `useCameraControls` hook

**UI Component Hierarchy (Room View):**
```
MainPage (landingPage.tsx)
  AppBar
    SiMGenButtons, DrawingModeToggle, EditingModeToggle
    ColorModeToggle, ConnectionDialog, ChatToggle, AddPlotButton
    ScreenshotButton, UploadButton, RoomManagementMenu, UserProfile
  SideBar
    PrimaryDrawer (icon navigation: selections, modifiers, analysis, geometries)
    SecondaryPanel (extension form with JSONForms, task submission)
    SelectionsPanel / GeometryPanel (context-specific panels)
  Canvas (R3F Canvas)
    Sphere, Bonds, Arrow, Box, Plane, Shape (instanced geometries)
    Camera, Lights, Floor, Fog, Cell (scene objects)
    Curve, CircleCurve, DrawingIndicator (drawing mode)
    MultiGeometryTransformControls, EditingIndicator (editing mode)
    CameraManager, OrbitControls (navigation)
    ScreenshotProvider, PathTracingRenderer (utilities)
  WindowManager (floating Plotly plot windows)
  FrameProgressBar (bottom timeline with bookmarks)
  ChatWindow (lazy-loaded overlay)
  ProgressNotifications (toast notifications)
```

**JSONForms Integration:**
- Purpose: Auto-generated forms from Pydantic JSON Schema (geometry editing, extension parameters)
- Renderers: `frontend/src/components/jsonforms-renderers/` (custom renderers for color pickers, SMILES editors, vec3 inputs, material selectors, etc.)
- Pattern: Backend geometry/extension Pydantic models -> JSON Schema -> JSONForms renders form -> user edits -> PUT/PATCH to backend

## Data Flow

**Frame Rendering Pipeline:**

1. User navigates to frame N (slider, playback, or socket event)
2. `useFrameBatch` hook requests `GET /v1/rooms/{roomId}/frames?indices=N&keys=arrays.positions,arrays.numbers,...`
3. Backend reads from `StorageRouter` -> `AsebytesStorage` (memory/LMDB/MongoDB)
4. If frame missing and provider mounted: dispatches `ProviderRequest` via Socket.IO, long-polls until result cached
5. Response: msgpack binary (`application/x-msgpack`)
6. Frontend decodes via `unpackBinary()` -> TypedArrays (Float32Array for positions, etc.)
7. React Query caches decoded data
8. `Particles` component reads positions/colors/numbers, renders instanced spheres via Three.js

**Geometry CRUD Flow:**

1. User edits geometry in JSONForms sidebar panel
2. Frontend calls `PUT /v1/rooms/{roomId}/geometries/{key}` with `{ type, data }`
3. Backend validates via Pydantic geometry model, stores in `RoomGeometry` SQL table
4. Backend broadcasts `GeometryInvalidate(operation="set", key=...)` via Socket.IO
5. All connected frontends receive event, refetch `GET /v1/rooms/{roomId}/geometries/{key}`
6. Updated geometry renders in Three.js scene

**Extension Execution Flow:**

1. User fills extension form (JSONForms from schema), clicks "Run"
2. Frontend calls `POST /v1/joblib/rooms/{roomId}/tasks/{jobName}` with payload
3. Backend creates `Task` record, dispatches to TaskIQ broker
4. TaskIQ worker picks up task, runs `InternalExtensionExecutor`
5. Executor creates `ZnDraw` client, instantiates extension class with payload, calls `extension.run(vis)`
6. Extension modifies frames/geometries/selections via `ZnDraw` client (REST API calls)
7. Each modification triggers Socket.IO broadcasts to all room members
8. Frontend receives invalidation events, refetches and re-renders

**Authentication Flow:**

1. Page loads -> `acquireToken()` checks localStorage for JWT
2. If token exists: validate via `GET /v1/auth/users/me`
3. If no token or invalid: `POST /v1/auth/guest` creates anonymous guest user
4. JWT stored in localStorage, sent as Bearer token on all REST requests
5. Socket.IO: `connectWithAuth()` passes JWT in `socket.auth = { token }`
6. Server validates JWT in `on_connect` handler, stores `user_id` in Socket.IO session

**Edit Lock Flow:**

1. User enters drawing or editing mode -> `acquireLock(msg)` called
2. Frontend: `PUT /v1/rooms/{roomId}/edit-lock` with optional Lock-Token header
3. Backend: stores lock in Redis with TTL (10s), returns `lock_token`
4. Frontend: starts 5s renewal interval, includes Lock-Token in all mutation requests
5. Backend: validates Lock-Token on all write endpoints (via `WritableRoomDep`)
6. Lock broadcast: `LockUpdate` event to all room members (UI shows lock indicator)
7. On mode exit: `DELETE /v1/rooms/{roomId}/edit-lock` releases lock
8. On disconnect: server auto-releases lock in `on_disconnect` handler

**State Management:**
- Persistent state (SQL): rooms, geometries, messages, bookmarks, figures, screenshots, presets, selection groups, memberships, server settings, jobs/tasks/workers
- Ephemeral state (Redis): session cameras, edit locks, progress trackers, provider frame counts, provider results/inflight
- Client state (Zustand): current frame, playing, mode, selections (optimistic), UI panels, connection status
- Cache (React Query): frame data, figure data, chat messages, job schemas, metadata

## Key Abstractions

**Room:**
- Purpose: Top-level container for a visualization session
- Backend: `Room` SQL model + Redis ephemeral state + storage frames
- Frontend: roomId in URL, all data scoped to room
- Pattern: Room ID is a UUID string, created on first visit if missing

**Geometry:**
- Purpose: 3D scene object with Pydantic-defined properties
- Backend: `BaseGeometry` subclasses in `src/zndraw/geometries/`, stored as JSON in `RoomGeometry.config`
- Frontend: rendered by matching Three.js component in `frontend/src/components/three/`
- Pattern: Type discriminator string maps to Pydantic model (backend) and React component (frontend)
- Special: Camera geometries stored in Redis (ephemeral session cameras) vs SQL (persistent shared cameras)

**Extension:**
- Purpose: Compute operation on frames/selections
- Backend: `Extension` subclass with `category` and `run(vis)` method
- Frontend: form rendered from JSON Schema, submitted as job task
- Dispatch: zndraw-joblib broker -> TaskIQ worker -> `InternalExtensionExecutor`

**Frame:**
- Purpose: Single snapshot of atomic system (positions, numbers, cell, properties)
- Format: `dict[bytes, bytes]` with msgpack-numpy encoded values
- Storage: `AsebytesStorage` (pluggable: memory/LMDB/MongoDB/H5MD)
- Wire format: msgpack binary response (`application/x-msgpack`)
- Client interface: `ZnDraw` implements `MutableSequence[ase.Atoms]`

## Entry Points

**ASGI Application:**
- Location: `src/zndraw/app.py` -> `socket_app`
- Triggers: uvicorn (via `main.py` or `zndraw` CLI)
- Responsibilities: Serves REST API, Socket.IO, and static frontend

**CLI (zndraw):**
- Location: `src/zndraw/cli.py`
- Entry point: `zndraw` console script
- Commands: start server, upload files, check status, connect, shutdown

**Agent CLI (zndraw-cli):**
- Location: `src/zndraw/cli_agent/__init__.py`
- Entry point: `zndraw-cli` console script
- Commands: structured JSON interface for LLM agents (rooms, frames, geometries, etc.)

**Database CLI (zndraw-db):**
- Location: `src/zndraw/cli.py` (`db_app`)
- Entry point: `zndraw-db` console script
- Commands: `db-init` for production multi-worker table initialization

**Python Client:**
- Location: `src/zndraw/client.py`
- Class: `ZnDraw(url, room, user, password)`
- Usage: `from zndraw import ZnDraw; vis = ZnDraw(url="http://localhost:8000")`

**Frontend Dev Server:**
- Location: `frontend/vite.config.ts`
- Command: `bun run dev` (Vite dev server with proxy to FastAPI at :8000)

## Error Handling

**Strategy:** RFC 9457 Problem Details for all HTTP/Socket.IO errors

**Patterns:**
- Backend: `ProblemType.exception(detail=...)` raises `ProblemException`
- Exception handler converts to `application/problem+json` response
- Socket.IO: `@tsio.exception_handler(ProblemException)` returns problem dict
- Frontend: axios interceptor handles 401 (token refresh), components handle specific status codes
- Python client: `ProblemType.raise_for_client()` maps to Python exceptions (`KeyError`, `IndexError`, `PermissionError`, etc.)

## Cross-Cutting Concerns

**Logging:** Python `logging` module, no structured logging framework

**Validation:** Pydantic models for all request bodies, `json_schema_extra` annotations for frontend form hints

**Authentication:** JWT tokens via fastapi-users (zndraw-auth package). Guest users auto-created. Superuser flag for admin operations.

**Authorization:** Room-level (membership, admin lock, edit lock), geometry-level (owner), task-level (worker identity)

**Serialization:** msgpack for frame data (binary efficiency), JSON for everything else, base64 for frame key encoding in REST

---

*Architecture analysis: 2026-03-05*
