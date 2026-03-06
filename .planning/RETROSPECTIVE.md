# Project Retrospective

*A living document updated after each milestone. Lessons feed forward into future planning.*

## Milestone: v1.0 -- Split Monolithic Files

**Shipped:** 2026-03-06
**Phases:** 3 | **Plans:** 8

### What Was Built
- Python client package (`src/zndraw/client/`) with 7 modules replacing 2297-line monolith
- Zustand scene sub-slices (4 files) replacing 596-line monolith
- Socket handler modules (7 files + types + utils) replacing 949-line hook
- 24 unit tests for serialization and exception hierarchy

### What Worked
- Risk-ordered phase execution (strongest tests first) caught import chain issues early in Phase 1
- Pure structural refactoring with zero functional changes made verification straightforward
- Handler factory pattern (HandlerContext DI) cleanly separated handler logic from React lifecycle
- Zustand intersection type composition preserved type safety without interface duplication

### What Was Inefficient
- E2E Playwright verification deferred to "human needed" in Phases 2-3 (no running server in execution context)
- ROADMAP.md state tracking fell behind -- Phase 3 still showed unchecked after completion
- SUMMARY one_liner frontmatter consistently null -- tooling didn't populate this field

### Patterns Established
- **Barrel index pattern**: Sub-modules export from barrel `index.ts`; consumers import from barrel only
- **Handler factory pattern**: `create*Handlers(ctx: HandlerContext)` returns `{ handlers, cleanup }`
- **Sub-slice composition**: `Type = SubA & SubB & SubC`, `createSlice = { ...createSubA(), ...createSubB() }`
- **Re-export layer**: `__init__.py` preserves all public import paths during package extraction

### Key Lessons
1. Pure structural refactors benefit from compiler verification (TypeScript `tsc --noEmit`) as primary confidence signal
2. Keep socket.on/off registration in a single useEffect -- splitting causes React StrictMode cleanup ordering bugs
3. TYPE_CHECKING guards cleanly solve circular import issues in Python package splits

### Cost Observations
- Model mix: ~70% sonnet (execution, verification), ~30% opus (planning, audit, integration check)
- Sessions: ~6 (codebase map, 3 phase executions, audit, completion)
- Notable: Sub-30min total execution time for 8 plans -- mechanical refactoring is fast with clear plans

---

## Cross-Milestone Trends

### Process Evolution

| Milestone | Phases | Plans | Key Change |
|-----------|--------|-------|------------|
| v1.0 | 3 | 8 | First milestone -- established barrel index and handler factory patterns |

### Cumulative Quality

| Milestone | Backend Tests | Frontend Tests | E2E Specs |
|-----------|---------------|----------------|-----------|
| v1.0 | 499 + 24 new | 0 (deferred to v2) | 13 (unchanged) |

### Top Lessons (Verified Across Milestones)

1. Risk-order phases by test coverage strength -- strongest safety net first
2. Pure structural refactors verify via compiler, not runtime -- save E2E for functional changes
