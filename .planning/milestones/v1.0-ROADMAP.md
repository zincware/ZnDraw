# Roadmap: Split Monolithic Files

## Overview

Three oversized files become three well-structured packages/directories: `client.py` (2297 lines) splits into a Python package, `sceneSlice.ts` (596 lines) decomposes into Zustand sub-slices, and `useSocketManager.ts` (949 lines) becomes domain-grouped handler modules. Each phase delivers one complete split, verified by existing tests before moving to the next. The order follows risk profile: strongest test coverage first (499 backend tests), then compiler-verified Zustand work, then the most nuanced React hook restructuring last.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Client Package** - Convert `client.py` into `src/zndraw/client/` package with single-responsibility modules and unit tests
- [x] **Phase 2: Scene Sub-Slices** - Decompose `sceneSlice.ts` into sub-slices following the existing Zustand slice composition pattern
- [ ] **Phase 3: Socket Handler Modules** - Extract `useSocketManager.ts` handlers into domain-grouped modules with typed interfaces

## Phase Details

### Phase 1: Client Package
**Goal**: The monolithic `client.py` is replaced by a well-organized `client/` package where each module has a single responsibility, all public imports are preserved, and extracted modules have unit test coverage
**Depends on**: Nothing (first phase)
**Requirements**: CLNT-01, CLNT-02, CLNT-03, CLNT-04, CLNT-05, CLNT-06, CLNT-07, CLNT-08, CLNT-09, CLNT-10
**Success Criteria** (what must be TRUE):
  1. `from zndraw import ZnDraw` and all other existing import paths work without changes to consuming code
  2. `src/zndraw/client/` is a package with separate modules for serialization, exceptions, lock, API management, socket management, and the main ZnDraw class
  3. All 499 existing backend tests pass without modification
  4. Unit tests exist for serialization helpers, exception classes, and other extracted modules
  5. The original monolithic `client.py` file no longer exists
**Plans:** 3 plans

Plans:
- [x] 01-01-PLAN.md -- Extract foundations: shared exceptions, serialization, lock, package scaffold
- [x] 01-02-PLAN.md -- Extract APIManager, SocketManager, ZnDraw core; wire re-exports; delete legacy file
- [x] 01-03-PLAN.md -- Add unit tests for serialization and exception hierarchy

### Phase 2: Scene Sub-Slices
**Goal**: The monolithic `sceneSlice.ts` is decomposed into focused sub-slices (geometry, selection, editing, drawing) that compose back into the existing `SceneSlice` interface without changing `store.tsx` or any consuming components
**Depends on**: Phase 1
**Requirements**: SCEN-01, SCEN-02, SCEN-03, SCEN-04, SCEN-05, SCEN-06, SCEN-07
**Success Criteria** (what must be TRUE):
  1. Four sub-slice files exist under a scene slices directory: geometry, selection, editing, and drawing
  2. `SceneSlice` type is composed from sub-slice interfaces and the composed type matches the original interface
  3. `createSceneSlice` spreads sub-slice creators following the same pattern used by `store.tsx` for top-level slices
  4. Camera state is handled explicitly (separate sub-slice or merged into geometry)
  5. All 12 Playwright E2E specs pass unchanged
**Plans:** 2 plans

Plans:
- [x] 02-01-PLAN.md -- Create four sub-slice files (geometry, selection, editing, drawing) plus barrel index.ts
- [x] 02-02-PLAN.md -- Wire store.tsx to scene barrel, delete old sceneSlice.ts, verify E2E

### Phase 3: Socket Handler Modules
**Goal**: The monolithic `useSocketManager.ts` is decomposed into domain-grouped handler modules with typed parameters, leaving a slim orchestrator that registers and cleans up handlers in a single `useEffect`
**Depends on**: Phase 2
**Requirements**: SOCK-01, SOCK-02, SOCK-03, SOCK-04, SOCK-05, SOCK-06, SOCK-07, SOCK-08, SOCK-09
**Success Criteria** (what must be TRUE):
  1. Handler modules exist for each domain: connection/lifecycle, frames, geometries, chat, scene invalidation, figures, and room/lock/progress
  2. The orchestrator file (`useSocketManager.ts`) is reduced to approximately 150 lines of handler registration and cleanup
  3. All handler parameters use typed interfaces instead of `any`
  4. All `socket.on()`/`socket.off()` calls remain in a single `useEffect` (no split registration)
  5. All Playwright E2E specs pass unchanged
**Plans:** 2/3 plans executed

Plans:
- [ ] 03-01-PLAN.md -- Create socketHandlers/ foundation (types, utils) and extract chat, scene invalidation, figure handlers
- [ ] 03-02-PLAN.md -- Extract connection, frame, geometry, and room handlers
- [ ] 03-03-PLAN.md -- Create barrel index.ts and rewrite useSocketManager.ts as slim orchestrator

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Client Package | 3/3 | Complete | - |
| 2. Scene Sub-Slices | 2/2 | Complete | 2026-03-06 |
| 3. Socket Handler Modules | 2/3 | In Progress|  |
