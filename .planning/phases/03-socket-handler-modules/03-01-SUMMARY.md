---
phase: 03-socket-handler-modules
plan: 01
subsystem: ui
tags: [typescript, zustand, socket.io, react-query, handler-factory]

# Dependency graph
requires:
  - phase: 02-scene-sub-slices
    provides: SceneSlice sub-slice types used in HandlerContext setter signatures
provides:
  - HandlerContext interface with ~30 typed dependency fields
  - createInvalidateHandler generic factory with getRoomId parameter
  - chatHandlers module with typed events and cleanup function
  - sceneInvalidationHandlers module with appStoreRoomId for schema queries
  - figureHandlers module with windowManagerStore imperative access
affects: [03-02, 03-03]

# Tech tracking
tech-stack:
  added: []
  patterns: [handler-factory-pattern, typed-socket-events, handler-context-DI]

key-files:
  created:
    - frontend/src/hooks/socketHandlers/types.ts
    - frontend/src/hooks/socketHandlers/utils.ts
    - frontend/src/hooks/socketHandlers/chatHandlers.ts
    - frontend/src/hooks/socketHandlers/sceneInvalidationHandlers.ts
    - frontend/src/hooks/socketHandlers/figureHandlers.ts
  modified: []

key-decisions:
  - "HandlerContext uses exact Zustand slice setter signatures (including optional source param on updateGeometry)"
  - "Pre-existing tsc errors in cliLoginApprove.tsx documented as out-of-scope (2 errors, unrelated to refactoring)"

patterns-established:
  - "Handler factory: createXxxHandlers(ctx: HandlerContext) returns object of named handler functions"
  - "Chat cleanup pattern: factory returns { handlers, cleanup } for modules with mutable state"
  - "Typed socket events: each module exports interfaces for its event payloads (no any on data params)"

requirements-completed: [SOCK-04, SOCK-05, SOCK-08]

# Metrics
duration: 2min
completed: 2026-03-06
---

# Phase 3 Plan 1: Socket Handler Foundation Summary

**HandlerContext type with ~30 fields, createInvalidateHandler utility, and 3 handler modules (chat, scene invalidation, figures) with typed event interfaces**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-06T09:07:50Z
- **Completed:** 2026-03-06T09:09:57Z
- **Tasks:** 2
- **Files created:** 5

## Accomplishments
- HandlerContext interface with all ~30 dependency fields matching exact Zustand slice signatures
- createInvalidateHandler extracted to utils.ts with getRoomId getter parameter (closure-free)
- Three handler modules extracted: chatHandlers, sceneInvalidationHandlers, figureHandlers
- All socket event parameters typed with exported interfaces (SOCK-08 partially fulfilled)
- Chat module owns typingTimeouts Map and exposes cleanup function for useEffect teardown

## Task Commits

Each task was committed atomically:

1. **Task 1: Create HandlerContext type and createInvalidateHandler utility** - `6c37413` (feat)
2. **Task 2: Extract chatHandlers, sceneInvalidationHandlers, and figureHandlers** - `0c9fcf2` (feat)

## Files Created/Modified
- `frontend/src/hooks/socketHandlers/types.ts` - HandlerContext interface with ~30 typed fields
- `frontend/src/hooks/socketHandlers/utils.ts` - createInvalidateHandler generic factory
- `frontend/src/hooks/socketHandlers/chatHandlers.ts` - Chat handlers with typingTimeouts cleanup
- `frontend/src/hooks/socketHandlers/sceneInvalidationHandlers.ts` - Scene invalidation handlers using appStoreRoomId
- `frontend/src/hooks/socketHandlers/figureHandlers.ts` - Figure handlers with windowManagerStore.getState()

## Decisions Made
- HandlerContext uses exact Zustand slice setter signatures (e.g., `updateGeometry` includes optional `source` param)
- Pre-existing tsc errors in `cliLoginApprove.tsx` (2 errors, `string | null` to `string | number | boolean`) documented as out-of-scope

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

- Pre-existing TypeScript errors in `src/pages/cliLoginApprove.tsx` (2 errors). Confirmed pre-existing by checking compilation on base commit. Not related to this refactoring.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Foundation files (types.ts, utils.ts) ready for use by remaining handler modules (03-02, 03-03)
- Handler factory pattern validated and ready for larger modules (connectionHandlers, frameHandlers, geometryHandlers, roomHandlers)
- No blockers for next plan

## Self-Check: PASSED

- All 5 created files verified present on disk
- Both task commits verified in git log (6c37413, 0c9fcf2)

---
*Phase: 03-socket-handler-modules*
*Completed: 2026-03-06*
