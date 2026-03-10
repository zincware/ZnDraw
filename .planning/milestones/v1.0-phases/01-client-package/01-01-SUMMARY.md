---
phase: 01-client-package
plan: 01
subsystem: api
tags: [python, packaging, exceptions, serialization, ase, dataclass]

# Dependency graph
requires: []
provides:
  - "src/zndraw/client/ package scaffold with __init__.py"
  - "Shared ZnDrawError and RoomLockedError in src/zndraw/exceptions.py"
  - "NotConnectedError in src/zndraw/client/exceptions.py"
  - "Frame serialization helpers in src/zndraw/client/serialization.py"
  - "ZnDrawLock context manager in src/zndraw/client/lock.py"
affects: [01-02, 01-03]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "client/ package structure with single-responsibility modules"
    - "Shared exceptions in zndraw.exceptions, client-specific in client/exceptions"
    - "TYPE_CHECKING guard for forward references to avoid circular imports"

key-files:
  created:
    - src/zndraw/client/__init__.py
    - src/zndraw/client/exceptions.py
    - src/zndraw/client/serialization.py
    - src/zndraw/client/lock.py
  modified:
    - src/zndraw/exceptions.py
    - src/zndraw/__init__.py
    - src/zndraw/_client_legacy.py

key-decisions:
  - "Renamed client.py to _client_legacy.py to coexist with client/ directory during incremental split"
  - "ZnDrawError and RoomLockedError placed before ProblemType definitions in exceptions.py"

patterns-established:
  - "Rename-then-create approach: rename monolith to _legacy, create package directory, update imports temporarily"
  - "TYPE_CHECKING guard for cross-module type hints (lock.py references APIManager)"

requirements-completed: [CLNT-01, CLNT-02, CLNT-03, CLNT-04]

# Metrics
duration: 2min
completed: 2026-03-05
---

# Phase 1 Plan 01: Extract Foundations Summary

**Shared exceptions moved to exceptions.py, client/ package created with serialization helpers, NotConnectedError, and ZnDrawLock as independent leaf modules**

## Performance

- **Duration:** 2 min
- **Started:** 2026-03-05T22:07:38Z
- **Completed:** 2026-03-05T22:10:18Z
- **Tasks:** 2
- **Files modified:** 7

## Accomplishments
- Eliminated lazy import cycle in RoomLocked.raise_for_client() by moving ZnDrawError and RoomLockedError to src/zndraw/exceptions.py
- Created src/zndraw/client/ package with four modules: __init__.py, exceptions.py, serialization.py, lock.py
- All new modules importable independently without circular import errors
- Legacy client still works via _client_legacy.py rename approach

## Task Commits

Each task was committed atomically:

1. **Task 1: Move shared exceptions to exceptions.py and create client package scaffold** - `544c59c` (feat)
2. **Task 2: Extract serialization helpers and ZnDrawLock to client submodules** - `c5bc508` (feat)

## Files Created/Modified
- `src/zndraw/client/__init__.py` - Package marker (empty at this stage)
- `src/zndraw/client/exceptions.py` - NotConnectedError inheriting from ZnDrawError
- `src/zndraw/client/serialization.py` - atoms_to_json_dict, json_dict_to_atoms, raw_frame_to_atoms, _estimate_frame_size, chunk constants
- `src/zndraw/client/lock.py` - ZnDrawLock dataclass context manager with auto-refresh
- `src/zndraw/exceptions.py` - Added ZnDrawError and RoomLockedError, removed lazy import cycle
- `src/zndraw/__init__.py` - Temporarily imports from _client_legacy
- `src/zndraw/_client_legacy.py` - Renamed from client.py, imports shared exceptions from zndraw.exceptions

## Decisions Made
- Renamed client.py to _client_legacy.py (instead of deleting) to allow the client/ directory to be created while keeping all existing imports working during the incremental split
- Placed ZnDrawError and RoomLockedError before ProblemType definitions in exceptions.py since they are plain exceptions, not problem types

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All leaf modules are in place for Plan 02 to extract APIManager, SocketManager, and ZnDraw core
- The _client_legacy.py file will be consumed by Plan 02 (extract remaining classes, wire re-exports, delete legacy)
- Plan 03 can then add unit tests for the extracted modules

## Self-Check: PASSED

All 6 created/modified files verified on disk. Both task commits (544c59c, c5bc508) verified in git log.

---
*Phase: 01-client-package*
*Completed: 2026-03-05*
