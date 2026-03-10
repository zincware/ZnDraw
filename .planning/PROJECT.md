# Split Monolithic Files

## What This Is

Refactor the three largest files in ZnDraw into logically grouped modules: `client.py` (backend), `useSocketManager.ts` (frontend), and `sceneSlice.ts` (frontend). Each monolith has been decomposed into single-responsibility modules.

## Core Value

Every extracted module has a single, clear responsibility -- files are grouped by concern, not by historical accident.

## Requirements

### Validated

- ✓ `client.py` split into `src/zndraw/client/` package with modules grouped by concern -- v1.0
- ✓ `useSocketManager.ts` split into grouped handler modules under `socketHandlers/` directory -- v1.0
- ✓ `sceneSlice.ts` split into dedicated sub-slices (geometry, selection, editing, drawing) -- v1.0
- ✓ Unit tests added for newly extracted modules (24 tests) -- v1.0
- ✓ All existing functionality works (499 backend tests pass) -- v1.0
- ✓ `ZnDraw` client accessible for programmatic use (import paths preserved) -- v1.0
- ✓ Socket.IO event handlers manage real-time sync in frontend (24 typed events) -- v1.0
- ✓ Zustand scene slice manages geometries, selections, editing, and drawing state -- v1.0

### Active

(None -- milestone complete. Use `/gsd:new-milestone` to define next work.)

### Out of Scope

- Functional changes to any module -- this was a pure structural refactor
- Other tech debt items from CONCERNS.md (type:ignore, except Exception, etc.)
- Frontend unit test infrastructure (Vitest) -- deferred to v2 (TEST-01, TEST-02, TEST-03)
- `APIManager` split by REST resource group -- deferred to v2 (DECO-01)
- Deprecation removal (`vis.log()`, `progress_bar`) -- deferred to v2 (DECO-02)
- Re-typing `geometries: Record<string, any>` -- separate concern

## Context

Shipped v1.0 with 12 feat commits across 29 files.
- Backend: FastAPI + python-socketio + SQLModel + Redis
- Frontend: React 19 + Zustand 5 (slice pattern) + Socket.IO client + Three.js
- Python client: `src/zndraw/client/` package (7 modules, 2364 lines total)
- Scene state: `frontend/src/stores/slices/scene/` (4 sub-slices, 67 members)
- Socket handlers: `frontend/src/hooks/socketHandlers/` (7 modules, 24 typed events)
- Tests: 499 backend (pytest) + 24 client unit tests, 13 E2E specs (Playwright)

## Constraints

- **Test suite**: All 499 existing tests must pass after each split
- **No functional changes**: Behavior must be identical before and after
- **Package manager**: `uv` for Python, `bun` for frontend

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Clean break on imports | Internal code, no external consumers | ✓ Good -- re-exports preserve all paths |
| Logical grouping over size limits | Cohesion matters more than line count | ✓ Good -- modules are cohesive |
| Add unit tests for extracted modules | Validate the split works | ✓ Good -- 24 tests catch regressions |
| Risk-ordered phases (backend→frontend) | Strongest test coverage first | ✓ Good -- backend tests caught issues early |
| Camera state in geometry sub-slice | Camera is geometry concern | ✓ Good -- clean separation |
| Handler factory pattern for socket modules | Context injection without globals | ✓ Good -- HandlerContext enables clean DI |
| Single useEffect for socket registration | Avoid React StrictMode double-mount bugs | ✓ Good -- 24 on/24 off symmetric cleanup |
| TYPE_CHECKING guard for ZnDraw forward ref | Avoid circular import in socket.py | ✓ Good -- no runtime circular imports |

---
*Last updated: 2026-03-06 after v1.0 milestone*
