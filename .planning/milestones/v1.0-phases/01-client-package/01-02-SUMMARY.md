---
phase: 01-client-package
plan: 02
subsystem: api
tags: [python, packaging, refactoring, dataclass, mutablesequence, socketio]

# Dependency graph
requires:
  - phase: 01-01
    provides: "client/ package scaffold with leaf modules (serialization.py, lock.py, exceptions.py)"
provides:
  - "APIManager class in src/zndraw/client/api.py"
  - "SocketManager class in src/zndraw/client/socket.py"
  - "ZnDraw main class in src/zndraw/client/core.py"
  - "Complete __init__.py re-export surface (14 names) for backwards compatibility"
  - "Clean import chain: from zndraw import ZnDraw -> from zndraw.client -> from zndraw.client.core"
affects: [01-03]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "TYPE_CHECKING guard in socket.py for ZnDraw forward reference (avoids circular import)"
    - "Lazy imports for ProviderTimeoutError in api.py get_frame/get_frames methods"
    - "__init__.py re-export layer preserves all previously-public import paths"

key-files:
  created:
    - src/zndraw/client/api.py
    - src/zndraw/client/socket.py
    - src/zndraw/client/core.py
  modified:
    - src/zndraw/client/__init__.py
    - src/zndraw/__init__.py

key-decisions:
  - "ProviderTimeoutError imported lazily inside get_frame/get_frames methods rather than at module top level, matching the codebase lazy-import pattern for worker dependencies"
  - "Kept _decode_raw_frame as a staticmethod on ZnDraw (not extracted to serialization.py) because it uses msgpack_numpy which is specific to the .get() API"

patterns-established:
  - "Module-per-class: api.py, socket.py, core.py each own exactly one @dataclass"
  - "Re-export __init__.py with __all__ ensures from zndraw.client import X backwards compatibility"

requirements-completed: [CLNT-05, CLNT-06, CLNT-07, CLNT-08, CLNT-09]

# Metrics
duration: 17min
completed: 2026-03-05
---

# Phase 1 Plan 02: Extract Core Classes Summary

**APIManager, SocketManager, and ZnDraw extracted to single-responsibility modules with full backwards-compatible re-export __init__.py and legacy file deletion**

## Performance

- **Duration:** 17 min
- **Started:** 2026-03-05T22:13:42Z
- **Completed:** 2026-03-05T22:30:57Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- Extracted APIManager (~880 lines) to src/zndraw/client/api.py with all REST API methods
- Extracted SocketManager (~130 lines) to src/zndraw/client/socket.py with TYPE_CHECKING import for ZnDraw to avoid circular import
- Extracted ZnDraw (~750 lines) to src/zndraw/client/core.py implementing MutableSequence[ase.Atoms]
- Wired client/__init__.py with complete re-export surface (14 names) preserving all existing import paths
- Updated src/zndraw/__init__.py to import from zndraw.client instead of _client_legacy
- Deleted src/zndraw/_client_legacy.py (2289 lines removed)
- All 941 tests pass (1 pre-existing failure in test_storage_asebytes unrelated to changes)

## Task Commits

Each task was committed atomically:

1. **Task 1: Extract APIManager and SocketManager to separate modules** - `2937bd9` (feat)
2. **Task 2: Extract ZnDraw core, wire __init__.py re-exports, delete legacy file** - `fc49f16` (feat)

## Files Created/Modified
- `src/zndraw/client/api.py` - APIManager dataclass with all REST API methods (auth, rooms, frames, geometries, bookmarks, figures, presets, sessions, chat, screenshots, progress, extensions, tasks)
- `src/zndraw/client/socket.py` - SocketManager dataclass with Socket.IO connection management, event handlers, and TYPE_CHECKING import for ZnDraw
- `src/zndraw/client/core.py` - ZnDraw main class implementing MutableSequence[ase.Atoms] with all properties, accessor patterns, mount/unmount, extension execution
- `src/zndraw/client/__init__.py` - Complete re-export surface: ZnDraw, APIManager, SocketManager, ZnDrawLock, NotConnectedError, ZnDrawError, RoomLockedError, Sessions, atoms_to_json_dict, json_dict_to_atoms, raw_frame_to_atoms, _estimate_frame_size, _TARGET_CHUNK_BYTES, _MAX_CHUNK_FRAMES
- `src/zndraw/__init__.py` - Changed import from _client_legacy to zndraw.client
- `src/zndraw/_client_legacy.py` - Deleted (was 2289 lines, now split across api.py + socket.py + core.py)

## Decisions Made
- ProviderTimeoutError imported lazily inside get_frame/get_frames methods rather than at module top level, matching the codebase lazy-import pattern for worker dependencies
- Kept _decode_raw_frame as a staticmethod on ZnDraw (not extracted to serialization.py) because it uses msgpack_numpy which is specific to the .get() API

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Pre-existing test failure in `tests/test_storage_asebytes.py::test_close_clears_rooms` (unrelated to refactoring, confirmed fails on prior commit too). Causes test pollution for `test_router_extend_delegates_to_default`. Documented in deferred-items.md.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All three core classes extracted and independently importable
- Plan 03 can add unit tests for the extracted modules
- No legacy file remnants -- clean package structure complete

## Self-Check: PASSED

All 5 created/modified files verified on disk. Legacy file deletion confirmed. Both task commits (2937bd9, fc49f16) verified in git log.

---
*Phase: 01-client-package*
*Completed: 2026-03-05*
