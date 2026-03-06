# Requirements: Split Monolithic Files

**Defined:** 2026-03-05
**Core Value:** Every extracted module has a single, clear responsibility — files grouped by concern, not by historical accident.

## v1 Requirements

Requirements for this milestone. Each maps to roadmap phases.

### Client Package

- [x] **CLNT-01**: `client.py` converted to `src/zndraw/client/` package with `__init__.py`
- [x] **CLNT-02**: Serialization helpers extracted to `client/serialization.py`
- [x] **CLNT-03**: Exception classes extracted to `client/exceptions.py`
- [x] **CLNT-04**: `ZnDrawLock` extracted to `client/lock.py`
- [x] **CLNT-05**: `APIManager` extracted to `client/api.py`
- [x] **CLNT-06**: `SocketManager` extracted to `client/socket.py`
- [x] **CLNT-07**: `ZnDraw` main class extracted to `client/core.py`
- [x] **CLNT-08**: `from zndraw import ZnDraw` continues to work
- [x] **CLNT-09**: All 499 existing tests pass unchanged
- [x] **CLNT-10**: Unit tests added for extracted modules

### Socket Handler Modules

- [ ] **SOCK-01**: Connection/lifecycle handlers extracted to separate module
- [ ] **SOCK-02**: Frame handlers extracted to separate module
- [ ] **SOCK-03**: Geometry handlers extracted to separate module
- [ ] **SOCK-04**: Chat handlers extracted to separate module
- [ ] **SOCK-05**: Scene invalidation handlers extracted to separate module
- [ ] **SOCK-06**: Room/lock/progress handlers extracted to separate module
- [ ] **SOCK-07**: Orchestrator hook reduced to ~150 lines of registration
- [ ] **SOCK-08**: Handler parameters use typed interfaces instead of `any`
- [ ] **SOCK-09**: 13 E2E Playwright specs pass unchanged

### Scene Sub-Slices

- [x] **SCEN-01**: Geometry sub-slice extracted (`geometrySlice.ts`)
- [x] **SCEN-02**: Selection sub-slice extracted (`selectionSlice.ts`)
- [x] **SCEN-03**: Editing sub-slice extracted (`editingSlice.ts`)
- [x] **SCEN-04**: Drawing sub-slice extracted (`drawingSlice.ts`)
- [x] **SCEN-05**: Camera concern handled (separate or merged into geometry)
- [x] **SCEN-06**: `SceneSlice` interface composed from sub-interfaces
- [x] **SCEN-07**: `createSceneSlice` composes sub-creators following existing pattern

## v2 Requirements

Deferred to future milestone. Tracked but not in current roadmap.

### Frontend Testing

- **TEST-01**: Vitest infrastructure setup for frontend unit tests
- **TEST-02**: Unit tests for extracted socket handler functions
- **TEST-03**: Unit tests for Zustand sub-slices

### Further Decomposition

- **DECO-01**: `APIManager` split by REST resource group (frames, geometries, rooms, etc.)
- **DECO-02**: Deprecation removal (`vis.log()`, `progress_bar`)

## Out of Scope

| Feature | Reason |
|---------|--------|
| Functional changes to any module | Pure structural refactor — behavior must be identical |
| Shared socket event type package (Python <-> TypeScript) | Requires codegen infrastructure, separate project |
| Splitting `useSocketManager` into multiple hooks | Would cause React StrictMode double-mount and cleanup ordering bugs |
| Converting Zustand slices to separate stores | Would break `AppState` composition and cross-slice `get()` |
| Moving `accessors.py` into client package | Already extracted, not part of this scope |
| Re-typing `geometries: Record<string, any>` | Separate concern, would make refactor diff unreviewable |

## Traceability

| Requirement | Phase | Status |
|-------------|-------|--------|
| CLNT-01 | Phase 1 | Complete |
| CLNT-02 | Phase 1 | Complete |
| CLNT-03 | Phase 1 | Complete |
| CLNT-04 | Phase 1 | Complete |
| CLNT-05 | Phase 1 | Complete |
| CLNT-06 | Phase 1 | Complete |
| CLNT-07 | Phase 1 | Complete |
| CLNT-08 | Phase 1 | Complete |
| CLNT-09 | Phase 1 | Complete |
| CLNT-10 | Phase 1 | Complete |
| SCEN-01 | Phase 2 | Complete |
| SCEN-02 | Phase 2 | Complete |
| SCEN-03 | Phase 2 | Complete |
| SCEN-04 | Phase 2 | Complete |
| SCEN-05 | Phase 2 | Complete |
| SCEN-06 | Phase 2 | Complete |
| SCEN-07 | Phase 2 | Complete |
| SOCK-01 | Phase 3 | Pending |
| SOCK-02 | Phase 3 | Pending |
| SOCK-03 | Phase 3 | Pending |
| SOCK-04 | Phase 3 | Pending |
| SOCK-05 | Phase 3 | Pending |
| SOCK-06 | Phase 3 | Pending |
| SOCK-07 | Phase 3 | Pending |
| SOCK-08 | Phase 3 | Pending |
| SOCK-09 | Phase 3 | Pending |

**Coverage:**
- v1 requirements: 26 total
- Mapped to phases: 26
- Unmapped: 0

---
*Requirements defined: 2026-03-05*
*Last updated: 2026-03-05 after roadmap creation*
