# Split Monolithic Files

## What This Is

Refactor the three largest files in ZnDraw into logically grouped modules: `client.py` (2297 lines, backend), `useSocketManager.ts` (949 lines, frontend), and `sceneSlice.ts` (596 lines, frontend). Each file combines multiple unrelated concerns that should be separate modules.

## Core Value

Every extracted module has a single, clear responsibility — files are grouped by concern, not by historical accident.

## Requirements

### Validated

- All existing ZnDraw functionality works (499 tests pass)
- `ZnDraw` client accessible for programmatic use
- Socket.IO event handlers manage real-time sync in frontend
- Zustand scene slice manages geometries, selections, editing, and drawing state

### Active

- [ ] `client.py` split into `src/zndraw/client/` package with modules grouped by concern
- [ ] `useSocketManager.ts` split into grouped handler modules under a hooks directory
- [ ] `sceneSlice.ts` split into dedicated slices (geometry, selection, editing, drawing)
- [ ] Unit tests added for newly extracted modules

### Out of Scope

- Functional changes to any of the three files — this is a pure structural refactor
- Other tech debt items from CONCERNS.md (type:ignore, except Exception, etc.)
- Frontend unit test infrastructure setup beyond what's needed for the split modules
- Changing public API contracts (REST endpoints, socket events, Zustand selectors used by components)

## Context

- Brownfield project with full codebase map at `.planning/codebase/`
- Backend: FastAPI + python-socketio + SQLModel + Redis
- Frontend: React 19 + Zustand 5 (slice pattern) + Socket.IO client + Three.js
- Python client uses httpx for REST, python-socketio for real-time
- Zustand store composed as `AppState = ConnectionSlice & PlaybackSlice & SceneSlice & LockSlice & UISlice`
- Clean break allowed on imports — no backwards-compatibility shims needed
- Tests: 499 backend tests (pytest), 13 E2E specs (Playwright), zero frontend unit tests

## Constraints

- **Test suite**: All 499 existing tests must pass after each split
- **No functional changes**: Behavior must be identical before and after
- **Package manager**: `uv` for Python, `bun` for frontend

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| Clean break on imports | Internal code, no external consumers | -- Pending |
| Logical grouping over size limits | Cohesion matters more than line count | -- Pending |
| Add unit tests for extracted modules | Validate the split works and prevents regression | -- Pending |

---
*Last updated: 2026-03-05 after initialization*
