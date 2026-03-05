# Project Research Summary

**Project:** Split Monolithic Files (ZnDraw)
**Domain:** Structural refactoring -- monolithic file decomposition in a Python/TypeScript codebase
**Researched:** 2026-03-05
**Confidence:** HIGH

## Executive Summary

This project is a pure structural refactor of three oversized files in a FastAPI + React/Zustand application: `client.py` (2297 lines), `useSocketManager.ts` (949 lines), and `sceneSlice.ts` (596 lines). No new libraries, no behavior changes, no API modifications. The "stack" is Python's module-to-package conversion pattern and Zustand's existing slice composition pattern -- both already proven in this codebase. The research confirms that all three files have clear, natural module boundaries marked by class definitions, section headers, and domain concerns. The split is mechanical, not architectural.

The recommended approach is sequential: split `client.py` first (clearest boundaries, 499 tests for instant verification, zero frontend tooling), then `sceneSlice.ts` (follows the codebase's existing Zustand slice pattern), then `useSocketManager.ts` last (most nuanced due to React effect lifecycle and handler dependency injection). Each split should be done in multiple commits -- structural setup first, then move code one concern at a time -- to preserve git history and enable incremental verification.

The primary risks are all related to import/export contract preservation. On the Python side, 20+ files import specific names from `zndraw.client`; missing a re-export in `__init__.py` breaks everything. On the TypeScript side, Zustand sub-slices must use `StateCreator<AppState>` (not the narrowed sub-slice type) to preserve cross-slice `get()` access, and socket handler registration/cleanup must remain centralized in a single `useEffect`. All of these risks are detectable immediately through compilation or test failures -- no silent bugs if verification is run after each commit.

## Key Findings

### Recommended Stack

No new dependencies are needed. This refactor uses only built-in language and framework mechanisms already present in the codebase. See [STACK.md](./STACK.md) for full details.

**Core techniques:**
- **Python `__init__.py` re-exports with `__all__`:** Convert `client.py` to `client/` package. Re-export all public names. Use tautological `as` syntax (`from .x import Y as Y`) for pyright compatibility. Already used by `geometries/` and `extensions/` in this codebase.
- **Zustand `StateCreator<AppState>` sub-slice composition:** Split `sceneSlice` into sub-slices using the same spread pattern that `store.tsx` uses for top-level slices. Cross-slice access via `get()` is preserved because all sub-slices are parameterized on full `AppState`.
- **Handler factory functions for socket events:** Extract handler logic from `useSocketManager` into domain-grouped factory functions. Keep a single `useEffect` orchestrator for registration/cleanup. Optionally adopt `useEffectEvent` (React 19.2) later for cleaner dependency management.
- **`TYPE_CHECKING` guards for circular imports:** `SocketManager` references `ZnDraw` and vice versa. The codebase already uses `from __future__ import annotations`. Use `TYPE_CHECKING` to break the circular import when these move to separate files.

### Expected Features

All deliverables are structural -- the "features" are the module boundaries themselves. See [FEATURES.md](./FEATURES.md) for full details.

**Must have (table stakes):**
- `client.py` decomposes into 6 modules: `_serialization.py`, `_exceptions.py`, `_lock.py`, `_api.py`, `_socket.py`, `_zndraw.py` plus `__init__.py` re-exports
- `sceneSlice.ts` decomposes into 4 sub-slices: geometry, selection, mode/editing, frame editing
- `useSocketManager.ts` decomposes into 8-9 handler modules plus slim orchestrator
- All 499 backend tests pass unchanged
- All 13 Playwright E2E specs pass unchanged
- All import paths preserved via re-exports

**Should have (differentiators):**
- Unit tests for extracted Python modules (serialization round-trips, APIManager dispatch)
- Typed handler signatures replacing `data: any` in socket handlers
- Removal of deprecated methods (`vis.log()`, `progress_bar`) during the split

**Defer (v2+):**
- Further decomposition of `APIManager` (850 lines but uniform pattern -- not worth splitting)
- Frontend unit test infrastructure (vitest setup is orthogonal to the split)
- Shared socket event type package between Python and TypeScript

### Architecture Approach

Three monolithic files become three packages/directories, each with an entry point that re-exports the public API. Internal modules follow a strict top-down dependency flow: leaf modules (serialization, exceptions, utilities) at the bottom, facade/orchestrator modules at the top. No module imports from a sibling at the same level -- cross-concern access uses runtime parameters (Python) or Zustand's `get()` (TypeScript). See [ARCHITECTURE.md](./ARCHITECTURE.md) for full details.

**Major components after split:**

1. **`src/zndraw/client/`** (Python package) -- 6 internal modules. `_api.py` is the leaf (HTTP client). `_zndraw.py` is the root (composes everything). Dependency flow is strictly one-directional.
2. **`frontend/src/stores/slices/sceneSlices/`** (TypeScript sub-slices) -- 4 sub-slice creators spread into existing `createSceneSlice`. `store.tsx` does not change.
3. **`frontend/src/hooks/socketHandlers/`** (TypeScript handler modules) -- 8-9 handler files grouped by domain. Single `useSocketManager.ts` orchestrator imports and registers all handlers.

### Critical Pitfalls

The top 5 risks, all with concrete prevention strategies. See [PITFALLS.md](./PITFALLS.md) for full details.

1. **Broken import contracts (Python)** -- 20+ files import from `zndraw.client`. Grep all imports before starting; create comprehensive `__init__.py` re-exports; run tests after the very first structural commit.
2. **Cross-slice type narrowing (Zustand)** -- Sub-slices must use `StateCreator<AppState, [], [], SubSlice>`, not `StateCreator<SubSlice>`. The wrong type parameter silently breaks `get().acquireLock()`, which could allow lock-free editing and data corruption.
3. **Socket handler registration/cleanup mismatch** -- Keep all `socket.on()`/`socket.off()` calls in a single `useEffect`. Extract handler *functions*, not handler *hooks*. Multiple `useEffect` calls cause registration ordering bugs and StrictMode double-mount issues.
4. **Circular imports between ZnDraw and SocketManager** -- Use `TYPE_CHECKING` guards. The codebase already has `from __future__ import annotations` which defers annotation evaluation. No runtime `isinstance` checks cross these boundaries.
5. **Re-export chain breakage (TypeScript)** -- `store.tsx` re-exports `getActiveCurves` and `selectPreferredCurve` from `sceneSlice.ts`. When helpers move to sub-slice files, update the re-export source path.

## Implications for Roadmap

Based on research, the project naturally divides into 3 phases ordered by risk profile and verification infrastructure.

### Phase 1: Split `client.py` into `client/` Package

**Rationale:** Clearest module boundaries (each class is self-contained). Richest test coverage (499 tests). Python's module-to-package conversion is the most straightforward of the three splits. No frontend tooling needed. Establishes the pattern for the project.

**Delivers:** `src/zndraw/client/` package with 6 internal modules and `__init__.py` re-exports. All public import paths preserved.

**Addresses:** Table stakes features: serialization extraction, exception extraction, lock extraction, APIManager extraction, SocketManager extraction, ZnDraw extraction.

**Avoids:** Pitfall 1 (broken imports -- mitigated by comprehensive re-exports), Pitfall 4 (circular imports -- mitigated by TYPE_CHECKING), Pitfall 9 (SocketManager coupling -- accepted, not fixed).

**Commit strategy:** (1) Create `client/` directory with `__init__.py` that imports from original `client.py` -- verify tests pass. (2) Extract leaf modules one at a time: `_serialization.py`, `_exceptions.py`, `_lock.py`. (3) Extract `_api.py`. (4) Extract `_socket.py` and `_zndraw.py` together (circular reference pair). (5) Remove original `client.py`.

### Phase 2: Split `sceneSlice.ts` into Sub-Slices

**Rationale:** Follows the exact Zustand slice composition pattern already used by 5 top-level slices in `store.tsx`. The interface definition in `sceneSlice.ts` already groups properties by concern. Only `store.tsx` imports from `sceneSlice.ts`, so the blast radius is minimal. Depends on understanding the pattern, which Phase 1 establishes.

**Delivers:** `frontend/src/stores/slices/sceneSlices/` directory with 4 sub-slice creators. `createSceneSlice` becomes a composition function that spreads sub-slices. `store.tsx` unchanged.

**Addresses:** Table stakes features: geometry slice, selection slice, mode/editing slice, frame edit slice, composed SceneSlice type.

**Avoids:** Pitfall 2 (cross-slice type narrowing -- all sub-slices use `StateCreator<AppState>`), Pitfall 5 (re-export chain -- update `store.tsx` re-exports), Pitfall 6 (duplicate property keys -- atomic moves, remove from old and add to new in same commit).

**Commit strategy:** (1) Create sub-slice directory and `geometrySlice.ts` (most independent, no cross-slice writes). (2) Extract `selectionSlice.ts`. (3) Extract `modeSlice.ts` (most complex -- lock acquisition, mode transitions). (4) Extract `frameEditSlice.ts`. (5) Update parent `sceneSlice.ts` to compose sub-slices.

### Phase 3: Split `useSocketManager.ts` into Handler Modules

**Rationale:** Most nuanced split. Handlers form closures over many dependencies. The handler factory pattern requires designing a dependency injection interface. However, this is lowest structural risk because handlers are pure side-effect functions that can be tested in isolation. Saved for last because the existing 949-line hook, while ugly, is functional.

**Delivers:** `frontend/src/hooks/socketHandlers/` directory with 8-9 handler modules. `useSocketManager.ts` slimmed to ~150 lines: just the `useEffect` orchestrator that imports and registers handlers.

**Addresses:** Table stakes features: connection handlers, frame handlers, geometry handlers, chat handlers, room handlers, lock handlers, progress handlers, cache handlers, orchestrator.

**Avoids:** Pitfall 3 (registration/cleanup mismatch -- centralized registration in single useEffect), Pitfall 7 (dependency array divergence -- preserve current pattern of stable setters + roomId/isOverview).

**Commit strategy:** (1) Create handler directory and `utils.ts` with shared `createInvalidateHandler`. (2) Extract smallest handler groups first: `lockHandlers.ts`, `progressHandlers.ts`, `chatHandlers.ts`. (3) Extract `frameHandlers.ts`, `geometryHandlers.ts`. (4) Extract `connectionHandlers.ts` last (most complex -- room creation, auth, initial state fetch). (5) Slim down orchestrator.

### Phase Ordering Rationale

- **Dependency order:** Phase 1 (Python) is completely independent of Phases 2-3 (TypeScript). Phases 2 and 3 are independent of each other but share the same verification infrastructure (Playwright E2E). Doing Phase 1 first means backend stability is locked in before touching the frontend.
- **Risk gradient:** Phase 1 has the best safety net (499 fast backend tests). Phase 2 has good safety (TypeScript compiler + 13 E2E specs). Phase 3 has adequate safety (same E2E specs) but the most subtle failure modes (stale closures, handler duplication).
- **Pattern establishment:** Phase 1 establishes the "extract module, re-export, verify" cadence that Phases 2 and 3 follow. The mental model transfers even though the languages differ.
- **Pitfall avoidance:** The most severe pitfalls (Pitfalls 1, 4) are in Phase 1, which has the strongest verification. The most subtle pitfalls (Pitfalls 3, 7) are in Phase 3, which is why it comes last -- by then the team has practiced the extract-verify loop.

### Research Flags

Phases likely needing deeper research during planning:
- **Phase 3 (useSocketManager split):** The handler factory dependency injection pattern needs concrete design. How many dependency objects? One per handler group or one shared? How are React Query `queryClient` and Zustand store actions passed? The current code uses closures -- the factory pattern is a meaningful design shift. Recommend `/gsd:research-phase` before implementation.

Phases with standard, well-documented patterns (skip research):
- **Phase 1 (client.py split):** Python module-to-package conversion is textbook. The codebase already has examples (`geometries/`, `extensions/`). No further research needed.
- **Phase 2 (sceneSlice split):** Zustand sub-slice composition is documented in the official wiki and already used in this codebase. No further research needed.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | No new technologies. All techniques are standard Python/TypeScript patterns already used in this codebase. |
| Features | HIGH | Based on direct line-by-line analysis of the three source files. Module boundaries align with class/concern boundaries. |
| Architecture | HIGH | Proposed architecture mirrors existing patterns in the codebase (package re-exports, slice composition, handler extraction). |
| Pitfalls | HIGH | All pitfalls identified from concrete import graph analysis and dependency tracing. Prevention strategies are verified against existing codebase patterns. |

**Overall confidence:** HIGH

This is a well-scoped structural refactor in a codebase with strong test coverage and established patterns for the exact decomposition techniques being applied. The research is based entirely on direct codebase analysis -- no external API documentation or novel technology patterns were needed.

### Gaps to Address

- **`useEffectEvent` adoption:** STACK.md recommends considering React 19.2's `useEffectEvent` for the socket handler split but at MEDIUM confidence. This should be evaluated during Phase 3 planning -- it is an optimization, not a requirement. The handler factory pattern works without it.
- **Frontend unit test setup:** FEATURES.md identifies vitest-based handler tests as a differentiator. The vitest infrastructure does not exist yet. Decide during Phase 3 planning whether to set it up as part of the split or defer entirely.
- **`Sessions` class export:** The import analysis found `from zndraw.client import Sessions` in one test file. This class is not prominently documented in the research. Verify during Phase 1 that `Sessions` is properly re-exported or that the test import is updated.

## Sources

### Primary (HIGH confidence)
- Direct analysis of `src/zndraw/client.py` (2297 lines) -- class boundaries, import graph, dependency flow
- Direct analysis of `frontend/src/hooks/useSocketManager.ts` (949 lines) -- handler grouping, dependency array, effect lifecycle
- Direct analysis of `frontend/src/stores/slices/sceneSlice.ts` (596 lines) -- state grouping, cross-slice access patterns
- Existing codebase patterns: `geometries/`, `extensions/`, `connectionSlice.ts`, `lockSlice.ts`, `store.tsx` -- verified in-repo precedents
- Import graph analysis via grep across `src/` and `tests/` directories

### Secondary (MEDIUM confidence)
- [Zustand slices pattern (DeepWiki)](https://deepwiki.com/pmndrs/zustand/7.1-slices-pattern) -- sub-slice composition
- [Python `__init__.py` re-exports (RealPython)](https://realpython.com/python-init-py/) -- package re-export patterns
- [React 19.2 useEffectEvent](https://react.dev/reference/react/useEffectEvent) -- stable callback refs (considered, not required)

### Tertiary (LOW confidence)
- None. All findings are backed by codebase analysis or official documentation.

---
*Research completed: 2026-03-05*
*Ready for roadmap: yes*
