# Phase 3: Socket Handler Modules - Research

**Researched:** 2026-03-06
**Domain:** React hook decomposition / Socket.IO event handler extraction
**Confidence:** HIGH

## Summary

Phase 3 is a pure structural refactor of `frontend/src/hooks/useSocketManager.ts` (949 lines) into 7 domain-grouped handler modules plus a slim orchestrator (~150 lines). The file contains 24 `socket.on()` registrations and 23 `socket.off()` calls inside a single `useEffect` with a 35-field dependency array. All handler functions are closures that capture store setters, `queryClient`, and `roomId` from the hook scope.

The refactoring pattern is well-established: extract handler functions into factory modules that receive a flat context object, then have the orchestrator build context, call factories, and register/deregister all handlers in one `useEffect`. Phase 2's barrel pattern (`stores/slices/scene/index.ts`) provides the exact template for directory structure and re-exports.

**Primary recommendation:** Create a `socketHandlers/` subdirectory under `hooks/`, define a single `HandlerContext` interface, implement 7 `createXxxHandlers(ctx)` factories, a shared `createInvalidateHandler` utility, and reduce `useSocketManager.ts` to context-building + handler registration.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- 7 handler modules (not 6): connectionHandlers, frameHandlers, geometryHandlers, chatHandlers, sceneInvalidationHandlers, figureHandlers, roomHandlers
- Selections, bookmarks, and cameras grouped with geometries (cameras are a geometry type per Phase 2)
- Figures get their own module (different store dependency: windowManagerStore)
- Room join flow extracted as separate `handleRoomJoin` function within the connection module (onConnect is ~230 lines)
- Retry state (retryDelay, backoff constants) owned by the connection module, not the orchestrator
- `cancelled` flag stays in the orchestrator -- passed to handlers via HandlerContext as `isCancelled: () => boolean`
- Single flat `HandlerContext` type with all dependencies (~30 fields: roomId, isOverview, isCancelled, queryClient, store setters)
- Each module exports a `createXxxHandlers(ctx: HandlerContext)` factory that returns an object of named handler functions
- No sub-types or per-module param types -- one flat context, modules pick what they need
- Chat module owns its mutable state (typingTimeouts Map) and provides a cleanup function
- `createInvalidateHandler` extracted to shared `socketHandlers/utils.ts`
- Each handler module defines and exports typed interfaces for its socket events (replaces `any` params -- fulfills SOCK-08)
- Existing chat types (MessageNewEvent, MessageEditedEvent) stay in `types/chat.ts` -- shared with useChat.ts
- New `frontend/src/hooks/socketHandlers/` subdirectory
- Barrel `index.ts` re-exports all create functions + HandlerContext type
- `types.ts` for HandlerContext definition
- `utils.ts` for createInvalidateHandler factory
- 7 handler module files
- `useSocketManager.ts` stays in `hooks/` and imports from `socketHandlers/`

### Claude's Discretion
- Exact field list for HandlerContext (which setters go in)
- Internal structure of each handler module (helper functions, const organization)
- Whether any handler needs additional store access via `useAppStore.getState()` beyond what's in context
- Barrel index.ts organization and re-export grouping

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| SOCK-01 | Connection/lifecycle handlers extracted to separate module | connectionHandlers.ts -- onConnect (with handleRoomJoin), onDisconnect, onConnectError; owns retry state |
| SOCK-02 | Frame handlers extracted to separate module | frameHandlers.ts -- onFrameUpdate, onFramesInvalidate, onFrameSelectionUpdate |
| SOCK-03 | Geometry handlers extracted to separate module | geometryHandlers.ts -- onGeometriesInvalidate, onSelectionsInvalidate, onSelectionGroupsInvalidate, onBookmarksInvalidate, onDefaultCameraInvalidate, onActiveCameraUpdate |
| SOCK-04 | Chat handlers extracted to separate module | chatHandlers.ts -- onChatMessageNew, onChatMessageUpdated, onTyping; owns typingTimeouts Map; provides cleanup function |
| SOCK-05 | Scene invalidation handlers extracted to separate module | sceneInvalidationHandlers.ts -- onInvalidate, onSchemaInvalidate |
| SOCK-06 | Room/lock/progress handlers extracted to separate module | roomHandlers.ts -- onRoomUpdate, onRoomDelete, onLockUpdate, onProgressStarted, onProgressUpdate, onProgressComplete |
| SOCK-07 | Orchestrator hook reduced to ~150 lines | useSocketManager.ts becomes: imports, HandlerContext construction, createXxxHandlers() calls, socket.on/off registration, connectWithAuth/cleanup |
| SOCK-08 | Handler parameters use typed interfaces instead of `any` | Each module exports typed interfaces for its event payloads; existing chat types reused from types/chat.ts |
| SOCK-09 | 13 E2E Playwright specs pass unchanged | Pure structural refactor -- no functional changes; verified by `tsc` + E2E suite |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| React | 19.x | Hook lifecycle (useEffect) | Project standard |
| socket.io-client | 4.x | Socket event binding | Project standard |
| zustand | 5.x | State management (useAppStore, useWindowManagerStore, useRoomsStore) | Project standard |
| @tanstack/react-query | 5.x | Query cache invalidation (queryClient) | Project standard |
| TypeScript | ESNext target | Type safety, interfaces | Project tsconfig.json |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| vite | 6.x | Build tool (bundles the refactored modules) | Build verification |
| biome | (scripts only) | Formatting and linting | Code style enforcement |

### Alternatives Considered
None -- all decisions are locked. This is a pure internal refactoring with no new dependencies.

## Architecture Patterns

### Recommended Project Structure
```
frontend/src/hooks/
├── useSocketManager.ts           # Slim orchestrator (~150 lines)
└── socketHandlers/
    ├── index.ts                  # Barrel re-exports
    ├── types.ts                  # HandlerContext interface
    ├── utils.ts                  # createInvalidateHandler factory
    ├── connectionHandlers.ts     # SOCK-01
    ├── frameHandlers.ts          # SOCK-02
    ├── geometryHandlers.ts       # SOCK-03
    ├── chatHandlers.ts           # SOCK-04
    ├── sceneInvalidationHandlers.ts  # SOCK-05
    ├── figureHandlers.ts         # (new -- split from scene)
    └── roomHandlers.ts           # SOCK-06
```

### Pattern 1: Handler Context (Dependency Injection)

**What:** A single flat interface carrying all dependencies that handler factories need.
**When to use:** Every handler factory receives this as its sole argument.
**Why flat:** The orchestrator already extracts ~30 selectors at the top of the hook. These become HandlerContext fields directly.

**HandlerContext fields (derived from current code analysis):**

```typescript
// Source: analysis of useSocketManager.ts lines 34-73
import type { QueryClient } from "@tanstack/react-query";
import type { GeometryData, GlobalSettings } from "../../myapi/client";
import type { UserInfo } from "../../utils/auth";
import type { InitializationError, Progress } from "../../store";

export interface HandlerContext {
  // Identity / routing
  roomId: string | undefined;
  appStoreRoomId: string | null;
  isOverview: boolean;
  isCancelled: () => boolean;

  // Query client
  queryClient: QueryClient;

  // Connection setters
  setConnected: (status: boolean) => void;
  setInitializationError: (error: InitializationError | null) => void;
  setUser: (user: UserInfo) => void;
  setSessionId: (sessionId: string | null) => void;
  setCameraKey: (cameraKey: string | null) => void;
  setServerVersion: (version: string | null) => void;
  setGlobalSettings: (settings: GlobalSettings | null) => void;

  // Playback setters
  setFrameCount: (count: number) => void;
  setCurrentFrame: (frame: number) => void;
  setFrameSelection: (selection: number[] | null) => void;
  setBookmarks: (bookmarks: Record<number, string> | null) => void;

  // Scene/geometry setters
  setSelections: (selections: Record<string, number[]>) => void;
  setSelectionGroups: (groups: Record<string, Record<string, number[]>>) => void;
  setGeometries: (geometries: Record<string, any>) => void;
  setGeometrySchemas: (schemas: Record<string, any>) => void;
  setGeometryDefaults: (defaults: Record<string, any>) => void;
  updateGeometry: (key: string, geometry: any) => void;
  removeGeometry: (key: string) => void;
  setActiveCurveForDrawing: (key: string | null) => void;

  // Lock setters
  setSuperuserLock: (locked: boolean) => void;
  setUserLock: (email: string | null, message?: string | null) => void;

  // Progress setters
  setProgressTrackers: (trackers: Record<string, Progress>) => void;
  addProgressTracker: (tracker: Progress) => void;
  updateProgressTracker: (update: Partial<Progress> & { progress_id: string }) => void;
  removeProgressTracker: (progressId: string) => void;

  // Window manager
  openWindow: (figureKey: string) => void;
}
```

**Note:** Some handlers also use `useAppStore.getState()` for imperative reads (e.g., `playing`, `chatOpen`, `sessionId`, `lockToken`, `activeCurveForDrawing`). These are NOT in the context -- they are read via `useAppStore.getState()` inside the handler at call time, which is the correct pattern for values that change during the effect lifetime.

### Pattern 2: Handler Factory

**What:** Each module exports a `createXxxHandlers(ctx: HandlerContext)` function returning a record of named handler functions.
**When to use:** Called once in the orchestrator's `useEffect` to produce the handler functions that get registered.

```typescript
// Source: project convention from CONTEXT.md decisions
interface ConnectionHandlers {
  onConnect: () => Promise<void>;
  onDisconnect: () => void;
  onConnectError: (err: Error) => Promise<void>;
}

export function createConnectionHandlers(ctx: HandlerContext): ConnectionHandlers {
  // Retry state owned by this module (not in HandlerContext)
  const INITIAL_RETRY_DELAY = 1000;
  const MAX_RETRY_DELAY = 30000;
  let retryDelay = INITIAL_RETRY_DELAY;

  async function handleRoomJoin(response: RoomJoinResponse): Promise<void> {
    // ... extracted from onConnect lines 161-287
  }

  return {
    onConnect: async () => { /* ... */ },
    onDisconnect: () => { /* ... */ },
    onConnectError: async (err: Error) => { /* ... */ },
  };
}
```

### Pattern 3: Chat Module with Cleanup

**What:** Chat handlers own mutable state (typingTimeouts Map) and return a cleanup function alongside handlers.
**When to use:** Only for chat module -- it's the only module with mutable state that needs explicit cleanup.

```typescript
interface ChatHandlers {
  onChatMessageNew: (data: MessageNewEvent) => void;
  onChatMessageUpdated: (data: MessageEditedEvent) => void;
  onTyping: (data: TypingEvent) => void;
}

interface ChatHandlersResult {
  handlers: ChatHandlers;
  cleanup: () => void;  // Clears typingTimeouts
}

export function createChatHandlers(ctx: HandlerContext): ChatHandlersResult {
  const typingTimeouts = new Map<string, ReturnType<typeof setTimeout>>();

  return {
    handlers: { /* ... */ },
    cleanup: () => {
      for (const timeout of typingTimeouts.values()) clearTimeout(timeout);
      typingTimeouts.clear();
    },
  };
}
```

### Pattern 4: Shared Invalidate Handler Factory

**What:** Generic factory for the common "fetch then update store" pattern.
**When to use:** Used by geometryHandlers (selections, bookmarks, selection groups) and potentially others.

```typescript
// Source: useSocketManager.ts lines 85-99
export function createInvalidateHandler<T>(
  fetchFn: (roomId: string) => Promise<T>,
  updateStoreFn: (data: T) => void,
  eventName: string,
  getRoomId: () => string | undefined,
): (data: unknown) => Promise<void> {
  return async () => {
    const roomId = getRoomId();
    if (!roomId) return;
    try {
      const response = await fetchFn(roomId);
      updateStoreFn(response);
    } catch (error) {
      console.error(`Error fetching ${eventName}:`, error);
    }
  };
}
```

**Key change from current code:** The current `createInvalidateHandler` closes over `roomId` from the `useEffect` closure. The extracted version takes a `getRoomId` function (or reads from `ctx.roomId`) to preserve the same behavior.

### Pattern 5: Orchestrator Registration

**What:** The slim orchestrator builds HandlerContext, calls all factories, and does bulk socket.on/off registration.
**When to use:** The core of the refactored `useSocketManager.ts`.

```typescript
// Source: project convention -- single useEffect pattern
useEffect(() => {
  let cancelled = false;
  const ctx: HandlerContext = {
    roomId, isOverview, isCancelled: () => cancelled,
    queryClient, setConnected, setFrameCount, /* ... all ~30 fields ... */
  };

  const connection = createConnectionHandlers(ctx);
  const frames = createFrameHandlers(ctx);
  const geometries = createGeometryHandlers(ctx);
  const { handlers: chat, cleanup: chatCleanup } = createChatHandlers(ctx);
  const sceneInvalidation = createSceneInvalidationHandlers(ctx);
  const figures = createFigureHandlers(ctx);
  const room = createRoomHandlers(ctx);

  // Register ALL handlers
  socket.on("connect", connection.onConnect);
  socket.on("disconnect", connection.onDisconnect);
  socket.on("connect_error", connection.onConnectError);
  // ... all 24 registrations ...

  // Auth/connect sequence (unchanged)
  if (socket.connected) {
    connection.onConnect();
  } else {
    connectWithAuth().then(/* ... */);
  }

  return () => {
    cancelled = true;
    // Cleanup logic (room_leave, sessionId clear, etc.)
    chatCleanup();
    socket.off("connect", connection.onConnect);
    // ... all socket.off() calls ...
  };
}, [/* same dependency array */]);
```

### Anti-Patterns to Avoid

- **Splitting socket.on/off across multiple useEffects:** Would cause React StrictMode double-mount issues and handler ordering bugs. All registration MUST stay in one useEffect (this is an explicit Out of Scope constraint).
- **Per-module context sub-types:** Adds complexity without benefit. One flat HandlerContext is simpler and matches the current code structure.
- **Putting mutable state in HandlerContext:** The `retryDelay` and `typingTimeouts` are module-internal. HandlerContext is for dependencies, not mutable state.
- **Importing useAppStore inside handler modules at module level:** Handler modules should NOT call hooks. They receive dependencies via context. However, `useAppStore.getState()` (imperative, non-hook access) IS allowed inside handler functions for reading current state at event-fire time.
- **Exporting handler functions directly (not via factory):** Would lose closure over mutable module state and the context dependency injection.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Query invalidation patterns | Custom cache busting | `queryClient.invalidateQueries()` / `queryClient.setQueryData()` | Already in use, battle-tested |
| Socket event typing | Runtime type checking | TypeScript interfaces on handler params | Compile-time safety, zero runtime cost |
| Module re-exports | Manual import paths everywhere | Barrel `index.ts` | Phase 2 established this pattern |
| Factory for fetch-then-update | Inline try/catch in each handler | `createInvalidateHandler` from utils.ts | DRY -- used by 3+ handlers |

## Common Pitfalls

### Pitfall 1: Stale Closure Over roomId
**What goes wrong:** Handler functions capture `roomId` at creation time. If the factory is called once but roomId changes, handlers use stale values.
**Why it happens:** The factory pattern creates closures in the useEffect body.
**How to avoid:** This is actually fine -- the useEffect dependency array includes `roomId`, so the entire effect re-runs when roomId changes, creating fresh handlers with the new roomId. This is the existing behavior and MUST be preserved.
**Warning signs:** If someone tries to "optimize" by memoizing handler factories outside the useEffect.

### Pitfall 2: Missing socket.off for Every socket.on
**What goes wrong:** Event handler leak -- old handler runs alongside new handler after re-render.
**Why it happens:** Forgetting to add a corresponding `socket.off()` in the cleanup function.
**How to avoid:** The orchestrator should register and deregister using the same handler reference. The factory return object provides stable references.
**Warning signs:** `socket.listeners("event_name").length > 1` in dev tools.

### Pitfall 3: Chat Cleanup Race
**What goes wrong:** Typing timeout fires after component unmount, calling `removeTypingUser` on unmounted store.
**Why it happens:** `setTimeout` callbacks outlive the useEffect cleanup if not explicitly cleared.
**How to avoid:** `chatCleanup()` MUST be called in the useEffect cleanup, which clears all pending timeouts.
**Warning signs:** Console warnings about state updates on unmounted components.

### Pitfall 4: Connection Handler Accessing Wrong `cancelled`
**What goes wrong:** onConnectError has early returns checking `cancelled`, but if `cancelled` is captured by value instead of by reference, it's always false.
**Why it happens:** JavaScript primitive capture semantics.
**How to avoid:** Pass `isCancelled: () => boolean` function in HandlerContext (a getter, not a value). This is already the decided approach.
**Warning signs:** Stale connection attempts continuing after component unmount.

### Pitfall 5: useAppStore.getState() vs Context Fields
**What goes wrong:** Reading a value from context that was stale (captured at effect creation time) when the handler needs the current value.
**Why it happens:** Some values change during the effect lifetime (e.g., `playing`, `chatOpen`, `sessionId`, `lockToken`).
**How to avoid:** Context fields are for dependencies that trigger effect re-creation (the useEffect dep array). For values that change independently, use `useAppStore.getState()` inside the handler at call time. The current code already does this correctly -- preserve the pattern.
**Warning signs:** `useAppStore.getState().playing` in onFrameUpdate, `useAppStore.getState().chatOpen` in onChatMessageNew -- these are correct and should NOT be moved to context.

### Pitfall 6: Breaking the useEffect Dependency Array
**What goes wrong:** TypeScript or React lint errors if the dependency array doesn't include all referenced values.
**Why it happens:** Extracting handler creation to factories means the factory calls still reference context fields that must be in deps.
**How to avoid:** The orchestrator still creates the HandlerContext inline in the useEffect, referencing all the same hook values. The dependency array stays the same (~35 entries).
**Warning signs:** React exhaustive-deps lint warnings.

### Pitfall 7: Import Cycles
**What goes wrong:** Handler module imports from store.tsx which imports from slices which might import from handlers.
**Why it happens:** Circular import chains.
**How to avoid:** Handler modules should only import types from store (not the store instance). They receive store setters via HandlerContext. The only exception is `useAppStore.getState()` for imperative reads -- this is a function call, not a module-level import of state.
**Warning signs:** Runtime "undefined" errors or bundler warnings about circular dependencies.

## Code Examples

### Handler Event Type Definitions (SOCK-08)

Each handler module defines typed interfaces for its socket event payloads:

```typescript
// connectionHandlers.ts
// Source: analysis of useSocketManager.ts lines 108-340

/** Room join socket response -- success case (no `status` field). */
export interface RoomJoinResponse {
  session_id: string;
  step: number;
  frame_count: number;
  camera_key: string | null;
  locked: boolean;
  progress_trackers?: Record<string, import("../../store").Progress>;
}

/** Room join socket response -- error case (RFC 9457 problem JSON). */
export interface RoomJoinError {
  status: number;
  detail?: string;
}
```

```typescript
// frameHandlers.ts
// Source: analysis of useSocketManager.ts lines 348-456

export interface FrameUpdateEvent {
  frame: number;
}

export interface FramesInvalidateEvent {
  room_id: string;
  action: "add" | "delete" | "modify" | "clear";
  indices?: number[];
  count?: number | null;
  reason?: string | null;
}

export interface FrameSelectionUpdateEvent {
  indices: number[] | null;
}
```

```typescript
// geometryHandlers.ts
// Source: analysis of useSocketManager.ts lines 527-718

export interface GeometryInvalidateEvent {
  operation?: "set" | "delete";
  key?: string;
}

export interface DefaultCameraInvalidateEvent {
  room_id: string;
  default_camera: string | null;
}

export interface ActiveCameraUpdateEvent {
  active_camera: string;
}
```

```typescript
// chatHandlers.ts -- reuses types from types/chat.ts
// Source: analysis of useSocketManager.ts lines 458-525

export interface TypingEvent {
  user_id: string;
  email: string;
  is_typing: boolean;
}
```

```typescript
// sceneInvalidationHandlers.ts
// Source: analysis of useSocketManager.ts lines 362-385

export interface InvalidateEvent {
  roomId: string;
  userName: string;
  category: string;
  extension: string;
  sessionId?: string;
}

export interface SchemaInvalidateEvent {
  category: string;
}
```

```typescript
// figureHandlers.ts
// Source: analysis of useSocketManager.ts lines 639-696

export interface FigureInvalidateEvent {
  key: string;
  operation?: "set" | "delete";
}
```

```typescript
// roomHandlers.ts
// Source: analysis of useSocketManager.ts lines 720-796

export interface RoomUpdateEvent {
  id: string;
  frame_count?: number | null;
  locked?: boolean | null;
  [key: string]: unknown;  // Room snapshot may contain other fields
}

export interface RoomDeleteEvent {
  room_id: string;
}

export interface LockUpdateEvent {
  action: "acquired" | "refreshed" | "released";
  user_id?: string | null;
  sid?: string;
  msg?: string | null;
  ttl?: number;
}

export interface ProgressStartedEvent {
  progress_id: string;
  description: string;
  unit?: string;
}

export interface ProgressUpdateEvent {
  progress_id: string;
  n?: number;
  total?: number | null;
  elapsed?: number;
}

export interface ProgressCompleteEvent {
  progress_id: string;
}
```

### Handler-to-Event Registration Map

Complete mapping of socket events to handler modules (24 events total):

| Event Name | Handler Function | Module |
|-----------|-----------------|--------|
| `connect` | `onConnect` | connectionHandlers |
| `disconnect` | `onDisconnect` | connectionHandlers |
| `connect_error` | `onConnectError` | connectionHandlers |
| `frame_update` | `onFrameUpdate` | frameHandlers |
| `frames_invalidate` | `onFramesInvalidate` | frameHandlers |
| `frame_selection_update` | `onFrameSelectionUpdate` | frameHandlers |
| `geometry_invalidate` | `onGeometriesInvalidate` | geometryHandlers |
| `selection_invalidate` | `onSelectionsInvalidate` | geometryHandlers |
| `selection_groups_invalidate` | `onSelectionGroupsInvalidate` | geometryHandlers |
| `bookmarks_invalidate` | `onBookmarksInvalidate` | geometryHandlers |
| `default_camera_invalidate` | `onDefaultCameraInvalidate` | geometryHandlers |
| `active_camera_update` | `onActiveCameraUpdate` | geometryHandlers |
| `message_new` | `onChatMessageNew` | chatHandlers |
| `message_edited` | `onChatMessageUpdated` | chatHandlers |
| `typing` | `onTyping` | chatHandlers |
| `invalidate` | `onInvalidate` | sceneInvalidationHandlers |
| `schema_invalidate` | `onSchemaInvalidate` | sceneInvalidationHandlers |
| `figure_invalidate` | `onFiguresInvalidate` | figureHandlers |
| `room_update` | `onRoomUpdate` | roomHandlers |
| `room_delete` | `onRoomDelete` | roomHandlers |
| `lock_update` | `onLockUpdate` | roomHandlers |
| `progress_start` | `onProgressStarted` | roomHandlers |
| `progress_update` | `onProgressUpdate` | roomHandlers |
| `progress_complete` | `onProgressComplete` | roomHandlers |

### Module Size Estimates

| Module | Handlers | Approx Lines | Notes |
|--------|----------|-------------|-------|
| connectionHandlers.ts | 3 | ~270 | Largest -- onConnect is ~230 lines with handleRoomJoin |
| frameHandlers.ts | 3 | ~110 | onFramesInvalidate is complex (query predicate logic) |
| geometryHandlers.ts | 6 | ~150 | Includes createInvalidateHandler usage + onGeometriesInvalidate (~100 lines) |
| chatHandlers.ts | 3 | ~80 | Includes typingTimeouts Map + cleanup |
| sceneInvalidationHandlers.ts | 2 | ~30 | Smallest -- two simple query invalidations |
| figureHandlers.ts | 1 | ~60 | Self-contained -- only uses windowManagerStore |
| roomHandlers.ts | 6 | ~90 | Mostly simple delegations to store |
| types.ts | - | ~50 | HandlerContext interface |
| utils.ts | - | ~20 | createInvalidateHandler generic |
| index.ts | - | ~15 | Barrel re-exports |
| **Total new files** | - | **~875** | |
| **Orchestrator** | - | **~150** | Context build + registration |

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Monolithic hook with all handlers inline | Domain-grouped factory modules | This phase | Readability, maintainability, testability |
| `any` typed event params | Typed interfaces per event | This phase | Compile-time safety, IDE autocomplete |
| `createInvalidateHandler` defined inside useEffect | Shared utility in utils.ts | This phase | DRY, reusable across modules |

**Not changing:**
- Socket.IO client library version -- stays at 4.x
- Single useEffect pattern -- required by React lifecycle
- `useAppStore.getState()` for imperative reads -- correct pattern for event handlers
- Dependency array entries -- same hook values referenced

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Playwright 1.58.x + tsc (TypeScript compiler) |
| Config file | `frontend/playwright.config.ts` |
| Quick run command | `cd frontend && npx tsc --noEmit` |
| Full suite command | `cd frontend && npx playwright test` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| SOCK-01 | Connection handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-02 | Frame handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-03 | Geometry handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-04 | Chat handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-05 | Scene invalidation handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-06 | Room/lock/progress handlers extracted | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-07 | Orchestrator reduced to ~150 lines | manual-only | `wc -l frontend/src/hooks/useSocketManager.ts` | N/A |
| SOCK-08 | Typed interfaces (no `any` params) | compilation | `cd frontend && npx tsc --noEmit` | N/A (structural) |
| SOCK-09 | 12 E2E Playwright specs pass | e2e | `cd frontend && npx playwright test` | 12 spec files exist |

**Note on SOCK-09:** The requirement states "13 Playwright E2E specs" but the repository contains 12 spec files. This may refer to test blocks rather than files, or the count may have changed. The verification criterion is: all existing E2E specs pass unchanged.

### Sampling Rate
- **Per task commit:** `cd frontend && npx tsc --noEmit` (type check ~10s)
- **Per wave merge:** `cd frontend && npx tsc --noEmit && npx vite build` (full build ~30s)
- **Phase gate:** E2E suite requires running server -- verify manually or skip per Phase 2 precedent (02-02 decision)

### Wave 0 Gaps
None -- existing test infrastructure covers all phase requirements. TypeScript compilation is the primary verification tool for a structural refactor. No new test files needed.

## Open Questions

1. **E2E spec count discrepancy**
   - What we know: Repository has 12 `.spec.ts` files in `frontend/e2e/`. Requirements say "13 Playwright E2E specs".
   - What's unclear: Whether the count refers to files or test blocks, or if a spec was added/removed.
   - Recommendation: Treat as "all existing E2E specs must pass" -- the exact count doesn't affect planning.

2. **`appStoreRoomId` vs `roomId` in HandlerContext**
   - What we know: The orchestrator has both `roomId` (from options or appStore) and `appStoreRoomId` (raw from store). `onSchemaInvalidate` uses `appStoreRoomId` while most others use `roomId`.
   - What's unclear: Whether this distinction is intentional or a bug.
   - Recommendation: Include both in HandlerContext to preserve exact behavior. The sceneInvalidationHandlers module needs `appStoreRoomId`.

3. **Connection handler imports**
   - What we know: `onConnect` imports from `../myapi/client` (createRoom, listGeometries, etc.), `../socket` (socket, connectWithAuth), and `../utils/*`.
   - What's unclear: Whether these should be passed via HandlerContext or imported directly by the module.
   - Recommendation: API functions and utils are stateless -- import them directly in the handler module. Only stateful/hook-derived values go in HandlerContext. This matches the principle that HandlerContext carries React state, not utility functions.

## Sources

### Primary (HIGH confidence)
- Direct code analysis of `frontend/src/hooks/useSocketManager.ts` (949 lines, 24 socket events)
- Direct code analysis of `frontend/src/stores/slices/*.ts` (ConnectionSlice, LockSlice, PlaybackSlice, UISlice, SceneSlice)
- Direct code analysis of `frontend/src/stores/windowManagerStore.ts`
- Direct code analysis of `frontend/src/types/chat.ts`
- Direct code analysis of `frontend/src/socket.ts`
- Direct code analysis of `frontend/src/roomsStore.tsx`
- Direct code analysis of `frontend/src/store.tsx`
- Phase 2 barrel pattern in `frontend/src/stores/slices/scene/index.ts`
- `frontend/tsconfig.json` -- ESNext target, strict mode, bundler resolution
- `frontend/package.json` scripts -- `tsc && vite build`
- `frontend/playwright.config.ts` -- 12 spec files, Chromium, no retries

### Secondary (MEDIUM confidence)
- Phase 3 CONTEXT.md decisions (user-locked, comprehensive)
- REQUIREMENTS.md SOCK-01 through SOCK-09 definitions
- STATE.md project history and Phase 2 precedent (skipping E2E in favor of tsc)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- no new libraries, pure refactoring of existing code
- Architecture: HIGH -- factory pattern is well-understood, CONTEXT.md decisions are comprehensive and specific
- Pitfalls: HIGH -- identified from direct code analysis of closures, cleanup, and React lifecycle
- Event types: HIGH -- derived directly from current handler parameter usage in source code

**Research date:** 2026-03-06
**Valid until:** Indefinite -- this is a structural refactoring of existing code with no external dependency concerns
