---
phase: 03-socket-handler-modules
plan: 02
subsystem: ui
tags: [typescript, socket.io, react-query, zustand, handler-pattern]

# Dependency graph
requires:
  - phase: 03-01
    provides: HandlerContext type and createInvalidateHandler utility
provides:
  - connectionHandlers factory (onConnect, onDisconnect, onConnectError)
  - frameHandlers factory (onFrameUpdate, onFramesInvalidate, onFrameSelectionUpdate)
  - geometryHandlers factory (6 handlers for geometry/selection/bookmark/camera invalidation)
  - roomHandlers factory (6 handlers for room/lock/progress state)
  - Typed event interfaces for all socket events
affects: [03-03]

# Tech tracking
tech-stack:
  added: []
  patterns: [handler-factory-with-typed-events, imperative-store-reads-via-getState]

key-files:
  created:
    - frontend/src/hooks/socketHandlers/connectionHandlers.ts
    - frontend/src/hooks/socketHandlers/frameHandlers.ts
    - frontend/src/hooks/socketHandlers/geometryHandlers.ts
    - frontend/src/hooks/socketHandlers/roomHandlers.ts
  modified: []

key-decisions:
  - "RoomJoinResponse union type with RoomJoinError replaces `any` for socket callback typing"
  - "RoomUpdateEvent uses index signature + as Room cast for setRoom since server sends full snapshots"

patterns-established:
  - "Handler factory: createXxxHandlers(ctx: HandlerContext) returns named handler functions"
  - "Imperative store reads: useAppStore.getState() for values that change during effect lifetime"
  - "createInvalidateHandler from utils for simple fetch-and-update handlers with getRoomId getter"

requirements-completed: [SOCK-01, SOCK-02, SOCK-03, SOCK-06, SOCK-08]

# Metrics
duration: 4min
completed: 2026-03-06
---

# Phase 03 Plan 02: Handler Modules Summary

**4 handler factories extracted with typed event interfaces: connection lifecycle (270+ lines with handleRoomJoin), frame invalidation with query predicates, geometry CRUD with cache management, and room/lock/progress state**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-06T09:08:07Z
- **Completed:** 2026-03-06T09:12:39Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Extracted connectionHandlers with handleRoomJoin as separate internal function and retry state owned by module
- Extracted frameHandlers preserving complex query predicate invalidation logic exactly
- Extracted geometryHandlers using createInvalidateHandler from utils for 3 handlers (selections, selectionGroups, bookmarks)
- Extracted roomHandlers correctly interacting with useRoomsStore and useAppStore via getState()
- All event parameters typed with exported interfaces (no `any` on event parameters)

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract connectionHandlers and frameHandlers** - `fdb7ea2` (feat)
2. **Task 2: Extract geometryHandlers and roomHandlers** - `639c505` (feat)

## Files Created/Modified
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts` - Connection lifecycle factory with handleRoomJoin, version checking, retry backoff
- `frontend/src/hooks/socketHandlers/frameHandlers.ts` - Frame update/invalidation factory with query predicate logic
- `frontend/src/hooks/socketHandlers/geometryHandlers.ts` - Geometry, selection, bookmark, camera invalidation factory
- `frontend/src/hooks/socketHandlers/roomHandlers.ts` - Room update/delete, lock, progress tracking factory

## Decisions Made
- Used `RoomJoinResponse | RoomJoinError` discriminated union (via `"status" in response`) to replace `any` for socket callback typing in connectionHandlers
- Used `as Room` cast in roomHandlers.setRoom call since RoomUpdateEvent has optional fields but server always sends full Room snapshots

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed RoomUpdateEvent to Room type mismatch in setRoom call**
- **Found during:** Task 2 (roomHandlers extraction)
- **Issue:** TypeScript error TS2345 -- RoomUpdateEvent has optional frame_count/locked but Room type requires them
- **Fix:** Added `as Room` cast with comment explaining server sends full snapshots; imported Room type
- **Files modified:** frontend/src/hooks/socketHandlers/roomHandlers.ts
- **Verification:** TypeScript compilation passes
- **Committed in:** 639c505 (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 type error from typed extraction)
**Impact on plan:** Necessary fix when replacing `any` with typed interfaces. No scope creep.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All 7 handler modules now exist under socketHandlers/ (types.ts, utils.ts, chatHandlers.ts, sceneInvalidationHandlers.ts, figureHandlers.ts + 4 new from this plan)
- Ready for Plan 03 (orchestrator wiring) to integrate these modules into useSocketManager.ts

## Self-Check: PASSED

- All 4 created files verified present on disk
- Both task commits (fdb7ea2, 639c505) verified in git log
- TypeScript compilation passes (only pre-existing errors in cliLoginApprove.tsx)

---
*Phase: 03-socket-handler-modules*
*Completed: 2026-03-06*
