# Coding Conventions

**Analysis Date:** 2026-03-05

## Backend (Python)

### Naming Patterns

**Files:**
- snake_case for all Python files: `socket_events.py`, `frame_source.py`
- Router files named by resource: `src/zndraw/routes/bookmarks.py`, `src/zndraw/routes/frames.py`
- Test files mirror source structure with `test_` prefix: `tests/test_routes_bookmarks.py`
- CLI agent tests in subdirectory: `tests/test_cli_agent/test_frames.py`

**Functions:**
- snake_case for all functions: `verify_room()`, `get_writable_room()`
- Prefix private/internal names with underscore: `_check_locks()`, `_camel_to_kebab()`
- Async functions keep same naming — no `async_` prefix
- FastAPI endpoint handlers named as verbs: `list_bookmarks()`, `set_bookmark()`, `delete_bookmark()`

**Variables:**
- snake_case for all variables and parameters
- Type aliases use PascalCase: `CurrentUserDep`, `RedisDep`, `WritableRoomDep`
- Constants are UPPER_SNAKE_CASE: `EDIT_LOCK_REFRESH = 5`, `PROBLEM_TYPES`

**Types/Classes:**
- PascalCase for all classes: `BaseGeometry`, `RoomBookmark`, `ProblemType`
- Pydantic models: PascalCase with descriptive suffix — `RoomCreate`, `BookmarkResponse`, `EditLockRequest`
- SQLModel tables: PascalCase singular nouns — `Room`, `Message`, `RoomMembership`, `Screenshot`
- Enums: PascalCase class, UPPER_SNAKE values: `MemberRole.OWNER`

### Code Style

**Formatting:**
- Tool: `ruff format` (via `uv run ruff format .`)
- Line length: 88 characters
- Quote style: double quotes
- Indent style: spaces (4 spaces)

**Linting:**
- Tool: `ruff` with extensive rule selection (see `pyproject.toml` `[tool.ruff.lint]`)
- Key enabled rule sets: E, W, F, I, B, C4, UP, ASYNC, S, SIM, RUF, PT, PERF, TCH, FA, FAST, N
- isort configured: `ruff check --select I --fix .`
- Fix import order: `uv run ruff check --select I --fix .`

**Type Checking:**
- Tool: `pyright` (via `uv run pyright .`)
- Known false positives: Redis async stubs for `smembers`/`srem`/`sadd`/`scard`
- Known false positives: `AsyncServerWrapper` dynamic attrs (`manager`, `manager_initialized`)
- Use `# type: ignore[...]` with specific error codes, never bare `# type: ignore`

### Import Organization

**Order (enforced by ruff isort):**
1. Standard library: `from collections.abc import AsyncIterator`
2. Third-party: `from fastapi import APIRouter`
3. First-party (`zndraw`): `from zndraw.models import Room`

**Path Aliases:**
- None — use relative or absolute `zndraw.*` imports
- `known-first-party = ["zndraw"]` in ruff config

**`from __future__ import annotations`:**
- Used in some modules (e.g., `src/zndraw/schemas.py`, `src/zndraw/exceptions.py`)
- Use when forward references are needed

**TYPE_CHECKING pattern:**
- Use `if TYPE_CHECKING:` for imports only needed for type hints: `src/zndraw/schemas.py` line 42

### Error Handling

**RFC 9457 Problem Details (backend standard):**
All HTTP errors use RFC 9457 `application/problem+json` format.

Define problem types as classes inheriting from `ProblemType` in `src/zndraw/exceptions.py`:

```python
class BookmarkNotFound(ProblemType):
    """The requested bookmark does not exist."""
    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)
```

Raise errors using the `.exception()` class method:
```python
raise BookmarkNotFound.exception(f"Bookmark '{index}' not found")
```

Document possible errors on endpoints with `problem_responses()`:
```python
@router.get("", responses=problem_responses(NotAuthenticated, RoomNotFound))
```

All problem types are registered in `PROBLEM_TYPES` dict at bottom of `src/zndraw/exceptions.py`.

**FastAPI validation errors** are also converted to RFC 9457 format in `src/zndraw/app.py` `_validation_exception_handler()`.

### Logging

**Framework:** Python `logging` (standard library)
- No custom logging framework — uses module-level loggers
- Frontend uses `console.error()`, `console.warn()`, `console.debug()`

### Comments

**When to Comment:**
- Use docstrings on all public functions, classes, and modules
- Inline comments for non-obvious logic or workarounds
- Section separators in long files: `# ====== Section Name ======`

**Docstrings:**
- Numpy-style docstrings for public APIs:
```python
def verify_room(session: AsyncSession, room_id: str) -> Room:
    """Verify room exists and return it, or raise RoomNotFound."""
```
- Concise single-line docstrings for simple functions
- Multi-line numpy-style for complex functions with parameters/returns

### Function Design

**Parameters:**
- Use FastAPI `Depends()` for dependency injection — annotated types like `CurrentUserDep`, `RedisDep`
- Max 10 args allowed for FastAPI endpoints (configured in ruff)
- Use `Path()` for path parameters with validation

**Return Values:**
- Always return Pydantic response models from endpoints, never raw dicts
- Use typed response models: `BookmarkResponse`, `StatusResponse`, `BookmarksResponse`

### Module Design

**Exports:**
- `__init__.py` files re-export key public API items
- Per-file ignores: `F401` in `__init__.py` for re-exports

**Barrel Files:**
- Used for geometry models: `src/zndraw/geometries/__init__.py`
- Used for extensions: `src/zndraw/extensions/__init__.py`

### Dependency Injection Pattern

**Define dependencies as functions + Annotated type aliases** in `src/zndraw/dependencies.py`:

```python
def get_redis(request: Request) -> AsyncRedis:
    return request.app.state.redis

RedisDep = Annotated[AsyncRedis, Depends(get_redis)]
```

Use in endpoints:
```python
async def list_bookmarks(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> BookmarksResponse:
```

### Pydantic Model Conventions

**All defaults in the Pydantic model** (single source of truth — CLAUDE.md rule):
```python
class Settings(BaseSettings):
    edit_lock_ttl: int = 10  # seconds
    storage: str = "memory://"
```

**Use `model_config = {"from_attributes": True}`** for ORM-backed response models.

**Use `Field()` for validation:**
```python
class FrameCreateRequest(BaseModel):
    frames: list[dict[str, Any]] = Field(min_length=1, max_length=1000)
```

### Router/Endpoint Conventions

**Router pattern** — each resource domain gets its own file in `src/zndraw/routes/`:
```python
router = APIRouter(prefix="/v1/rooms/{room_id}/bookmarks", tags=["bookmarks"])
```

**URL structure:** `/v1/rooms/{room_id}/{resource}` for room-scoped resources.

**Socket.IO events** use Pydantic models defined in `src/zndraw/socket_events.py`. Event names derived from class name in snake_case: `BookmarksInvalidate` -> `bookmarks_invalidate`.

**Broadcast pattern** after mutations:
```python
await sio.emit(
    BookmarksInvalidate(room_id=room_id, index=index, operation="set"),
    room=room_channel(room_id),
)
```

---

## Frontend (TypeScript/React)

### Naming Patterns

**Files:**
- Components: PascalCase `.tsx`: `ChatWindow.tsx`, `SelectionsPanel.tsx`, `PathTracingRenderer.tsx`
- Hooks: camelCase with `use` prefix: `useChat.ts`, `useGeometries.ts`, `useSocketManager.ts`
- Stores: camelCase with `Store` suffix: `geometryStore.ts`, `formStore.ts`, `windowManagerStore.ts`
- Slices: camelCase with `Slice` suffix: `sceneSlice.ts`, `connectionSlice.ts`, `uiSlice.ts`
- Utils: camelCase: `colorUtils.ts`, `cameraUtils.ts`, `msgpack-numpy.ts`
- Types: camelCase: `chat.ts`, `jobs.ts`, `property-inspector.ts`
- Constants: camelCase: `layout.ts`, `fileTypes.ts`
- Pages: camelCase: `landingPage.tsx`, `roomList.tsx`, `templateSelection.tsx`

**Functions/Variables:**
- camelCase for functions and variables: `handleSendMessage`, `getInitialPosition`
- Constants: UPPER_SNAKE_CASE or `as const` objects: `LAYOUT_CONSTANTS`, `CHAT_WIDTH`
- React hooks: `use` prefix: `useAppStore`, `useChatMessages`
- Event handlers: `handle` prefix for React handlers, `on` prefix for socket handlers

**Types/Interfaces:**
- PascalCase: `ChatMessage`, `GeometryData`, `AppState`
- Interface for component props: `ChatWindowProps`, `SocketManagerOptions`
- Type for union/utility types: `GeometryMode`, `ProgressColor`

### Code Style

**Formatting:**
- Tool: Biome (`biome format --write .`) and Prettier (both available)
- Scripts: `bun run format` / `bun run format:check`
- Indent: tabs (Biome default)
- Semicolons: yes

**Linting:**
- Tool: Biome (`biome lint .` / `biome check --fix .`)
- Scripts: `bun run lint` / `bun run lint:fix`

**TypeScript:**
- Strict mode enabled (`"strict": true` in `tsconfig.json`)
- Target: ESNext
- Module resolution: bundler
- JSX: react-jsx
- Path alias: `@/*` -> `./src/*`

### Import Organization

**Order (observed pattern):**
1. React / React ecosystem (`react`, `react-dom`, `react-router-dom`)
2. Third-party libraries (`@mui/*`, `@tanstack/*`, `zustand`, `socket.io-client`)
3. Local imports — stores, hooks, utils, types
4. CSS imports last

**Path Aliases:**
- `@/*` maps to `./src/*` (configured in `tsconfig.json`)
- Usage observed but not ubiquitous — relative imports also used

### Component Patterns

**Functional components only** — no class components:
```typescript
const ChatWindow = ({ open, onClose }: ChatWindowProps) => {
    // ...
};
export default memo(ChatWindow);
```

**Performance optimization:**
- Use `memo()` for components that re-render frequently: `export default memo(ChatWindow);`
- Use `useCallback` for stable event handler references
- Use individual Zustand selectors to prevent unnecessary re-renders:
  ```typescript
  const setCurrentFrame = useAppStore((state) => state.setCurrentFrame);
  const resetChatUnread = useAppStore((state) => state.resetChatUnread);
  ```

**Page components** in `frontend/src/pages/` — one per route:
- `landingPage.tsx` -> `/rooms/:roomId`
- `roomList.tsx` -> `/rooms`
- `templateSelection.tsx` -> `/`

### State Management

**Zustand** is the primary state management library (`frontend/src/store.tsx`).

**Slice pattern** — large store split into domain slices:
```
frontend/src/stores/slices/
  connectionSlice.ts  # roomId, user, session, server version
  playbackSlice.ts    # frame navigation, playback
  sceneSlice.ts       # geometries, selections, editing modes
  lockSlice.ts        # edit lock management
  uiSlice.ts          # chat, snackbar, progress, screenshots
```

Composed into a single store:
```typescript
export const useAppStore = create<AppState>((...a) => ({
    ...createConnectionSlice(...a),
    ...createPlaybackSlice(...a),
    ...createSceneSlice(...a),
    ...createLockSlice(...a),
    ...createUISlice(...a),
}));
```

**Slice definition pattern:**
```typescript
export interface ConnectionSlice {
    // State
    roomId: string | null;
    // Actions
    setRoomId: (roomId: string) => void;
}

export const createConnectionSlice: StateCreator<AppState, [], [], ConnectionSlice> = (set) => ({
    roomId: null,
    setRoomId: (roomId) => set({ roomId }),
});
```

**Domain-specific stores** use Zustand + Immer middleware:
```typescript
export const useGeometryStore = create<GeometryState & GeometryActions>()(
    immer((set) => ({
        mode: "list",
        setMode: (mode) => { set((state) => { state.mode = mode; }); },
    })),
);
```

**Derived selectors** defined as standalone functions:
```typescript
export const selectIsRoomReadOnly = (state: AppState): boolean => { ... };
```

### Data Fetching

**TanStack React Query** for server state:
```typescript
const queryClient = new QueryClient({
    defaultOptions: {
        queries: { staleTime: 30000, gcTime: 5 * 60 * 1000 },
    },
});
```

**Custom hooks** wrap React Query for each resource domain in `frontend/src/hooks/`:
```typescript
export const useGeometriesList = (roomId: string | null) => {
    return useQuery<GeometryListResponse>({
        queryKey: ["geometries", roomId, "list"],
        queryFn: () => listGeometries(roomId!),
        enabled: !!roomId,
        staleTime: 30000,
    });
};
```

**Query key conventions:**
- `["geometries", roomId, "list"]` — list queries
- `["geometries", roomId, "detail", key]` — detail queries
- `["chat", roomId]` — domain-scoped queries
- `["frame", roomId, frameIndex, key]` — frame data queries

**Mutations with optimistic updates:**
```typescript
export const useCreateGeometry = () => {
    const queryClient = useQueryClient();
    return useMutation({
        mutationFn: (...) => createGeometry(...),
        onMutate: async (variables) => {
            // Cancel, snapshot, optimistically update
        },
        onError: (_err, variables, context) => {
            // Rollback
        },
        onSettled: (_data, _err, variables) => {
            // Invalidate
        },
    });
};
```

### API Client

**Single API client** in `frontend/src/myapi/client.ts` using Axios:
- Interceptors for JWT auth, session ID, lock token headers
- Automatic 401 retry with token refresh
- All API functions exported as named exports: `listGeometries()`, `createRoom()`
- Response types defined alongside functions

### Socket.IO Integration

**Single socket instance** in `frontend/src/socket.ts`:
```typescript
export const socket = io(undefined, { autoConnect: false });
```

**Central socket manager hook** `frontend/src/hooks/useSocketManager.ts`:
- Registers all socket event handlers
- Handles connect/disconnect/reconnect lifecycle
- Updates Zustand store and React Query cache from socket events
- Uses invalidation pattern — socket events trigger REST re-fetches

### Styling

**MUI (Material UI)** is the component library:
- Components: `@mui/material`, `@mui/icons-material`
- Data display: `@mui/x-data-grid`, `@mui/x-tree-view`
- Theme: created with `createTheme()` supporting light/dark mode

**Inline `sx` prop** for component-specific styles (no CSS modules):
```typescript
<Box sx={{ p: 2, display: "flex", gap: 1, alignItems: "center" }}>
```

**CSS file** `frontend/src/index.css` for global styles only (body, root, animations).

**Layout constants** centralized in `frontend/src/constants/layout.ts`:
```typescript
export const LAYOUT_CONSTANTS = {
    APPBAR_HEIGHT: 50,
    PROGRESSBAR_HEIGHT: 45,
    PRIMARY_DRAWER_WIDTH: 50,
    SECONDARY_DRAWER_WIDTH: 600,
} as const;
```

### 3D Rendering

**React Three Fiber** (`@react-three/fiber`) + **Drei** (`@react-three/drei`):
- Three.js components in `frontend/src/components/three/`
- Each geometry type gets its own component: `Particles.tsx`, `Bonds.tsx`, `Camera.tsx`, `Curve.tsx`
- Path tracing via `@react-three/gpu-pathtracer`

### Data Flow: Backend to Frontend

1. Frontend connects via Socket.IO with JWT auth
2. Socket `room_join` event returns initial state (step, frame count, etc.)
3. Geometries and selections fetched via REST (`GET /v1/rooms/{roomId}/geometries`)
4. Frame binary data fetched via REST with MessagePack encoding
5. Server mutations broadcast invalidation events via Socket.IO
6. Frontend re-fetches data via REST on invalidation (not via socket payload)
7. React Query caches all REST responses, socket events invalidate cache keys

---

*Convention analysis: 2026-03-05*
