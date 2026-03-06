# Phase 3: Socket Handler Modules - Context

**Gathered:** 2026-03-06
**Status:** Ready for planning

<domain>
## Phase Boundary

Decompose `useSocketManager.ts` (949 lines) into domain-grouped handler modules with typed parameters, leaving a slim orchestrator (~150 lines) that registers and cleans up handlers in a single `useEffect`. Pure structural refactor — no functional changes. All 13 Playwright E2E specs must pass unchanged.

</domain>

<decisions>
## Implementation Decisions

### Domain grouping
- 7 handler modules (not the 6 originally listed in roadmap):
  1. **connectionHandlers** — onConnect (with extracted handleRoomJoin), onDisconnect, onConnectError
  2. **frameHandlers** — onFrameUpdate, onFramesInvalidate, onFrameSelectionUpdate
  3. **geometryHandlers** — onGeometriesInvalidate, onSelectionsInvalidate, onSelectionGroupsInvalidate, onBookmarksInvalidate, onDefaultCameraInvalidate, onActiveCameraUpdate
  4. **chatHandlers** — onChatMessageNew, onChatMessageUpdated, onTyping (owns typingTimeouts Map + cleanup)
  5. **sceneInvalidationHandlers** — onInvalidate, onSchemaInvalidate (small module, kept separate for clarity)
  6. **figureHandlers** — onFiguresInvalidate (own module — interacts with windowManagerStore, self-contained)
  7. **roomHandlers** — onRoomUpdate, onRoomDelete, onLockUpdate, onProgressStarted, onProgressUpdate, onProgressComplete
- Selections, bookmarks, and cameras grouped with geometries (cameras are a geometry type per Phase 2)
- Figures get their own module (different store dependency: windowManagerStore)

### Connection handler design
- Room join flow extracted as separate `handleRoomJoin` function within the connection module (onConnect is ~230 lines)
- Retry state (retryDelay, backoff constants) owned by the connection module, not the orchestrator
- `cancelled` flag stays in the orchestrator — passed to handlers via HandlerContext as `isCancelled: () => boolean`

### Dependency passing
- Single flat `HandlerContext` type with all dependencies (~30 fields: roomId, isOverview, isCancelled, queryClient, store setters)
- Each module exports a `createXxxHandlers(ctx: HandlerContext)` factory that returns an object of named handler functions
- No sub-types or per-module param types — one flat context, modules pick what they need
- Chat module owns its mutable state (typingTimeouts Map) and provides a cleanup function

### Factory pattern
- `createInvalidateHandler` extracted to shared `socketHandlers/utils.ts`
- Used by geometryHandlers (selections, bookmarks) and potentially others

### Event types
- Each handler module defines and exports typed interfaces for its socket events (replaces `any` params — fulfills SOCK-08)
- Existing chat types (MessageNewEvent, MessageEditedEvent) stay in `types/chat.ts` — shared with useChat.ts

### Directory structure
- New `frontend/src/hooks/socketHandlers/` subdirectory (follows Phase 2's scene/ barrel pattern)
- Barrel `index.ts` re-exports all create functions + HandlerContext type
- `types.ts` for HandlerContext definition
- `utils.ts` for createInvalidateHandler factory
- 7 handler module files (connectionHandlers.ts, frameHandlers.ts, etc.)
- `useSocketManager.ts` stays in `hooks/` and imports from `socketHandlers/`

### Claude's Discretion
- Exact field list for HandlerContext (which setters go in)
- Internal structure of each handler module (helper functions, const organization)
- Whether any handler needs additional store access via `useAppStore.getState()` beyond what's in context
- Barrel index.ts organization and re-export grouping

</decisions>

<code_context>
## Existing Code Insights

### Reusable Assets
- `createInvalidateHandler` factory (lines 85-99): Generic fetch-then-update pattern — will move to shared utils
- `types/chat.ts`: MessageNewEvent, MessageEditedEvent already typed — chat handlers import these
- `socket.ts`: Exports `connectWithAuth` and `socket` singleton — connection module imports these
- `windowManagerStore.ts`: Used only by figure handlers — clean separation point

### Established Patterns
- Phase 2 barrel pattern: `stores/slices/scene/index.ts` re-exports sub-slice creators and types
- `StateCreator` / factory pattern: Sub-slices use `(set, get)` factory — handler modules will use `(ctx: HandlerContext)` factory
- `useAppStore.getState()` for imperative access inside callbacks (already used in onConnect for showSnackbar, resetChatUnread)
- Single useEffect with socket.on/socket.off registration — MUST be preserved (REQUIREMENTS out-of-scope)

### Integration Points
- `useSocketManager.ts` (orchestrator): Builds HandlerContext, calls create functions, registers all handlers in one useEffect
- `useAppStore` selectors (lines 34-67): ~30 selectors extracted at hook top — these become HandlerContext fields
- `useQueryClient()` (line 68): Passed via HandlerContext
- `useWindowManagerStore` (line 69): Only needed by figure handlers — passed via HandlerContext or direct import
- `useRoomsStore` (line 15): Used by onRoomUpdate/onRoomDelete — room handlers import directly or receive via context

</code_context>

<specifics>
## Specific Ideas

No specific requirements — open to standard approaches

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 03-socket-handler-modules*
*Context gathered: 2026-03-06*
