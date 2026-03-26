# Repository Consolidation: auth + joblib into zndraw

**Date:** 2026-03-25
**Scope:** zndraw-fastapi only — copy source from zndraw-auth and zndraw-joblib, no changes to those repos
**Goal:** Consolidate zndraw-auth and zndraw-joblib into the zndraw-fastapi repo as colocated top-level packages, preserving import compatibility, test isolation, and separation of concerns.

## Problem

Four repos serve one product:

| Repo | Source LOC | Standalone? |
|------|-----------|-------------|
| zndraw-fastapi | ~20k | main app |
| zndraw-auth | ~870 | no — useless without zndraw |
| zndraw-joblib | ~3,600 | no — useless without zndraw |
| zndraw-socketio | 2,648 | **yes** — generic typed wrapper |

The auth and joblib packages cannot function independently, yet they live in separate repos with separate CI, separate releases, and separate version management. This creates friction:

1. **Cross-repo changes** require coordinated PRs and version bumps across 2-3 repos.
2. **CI duplication** — auth has no CI at all (relies on zndraw-fastapi to catch breakage); joblib has its own GitHub Actions.
3. **Dependency version drift** — each repo pins its own versions of shared deps (fastapi, sqlmodel, pydantic-settings).
4. **Development friction** — editable installs with relative paths break in worktrees; developers must manage 3-4 repos in sync.

zndraw-socketio is genuinely standalone (zero zndraw imports, Python >=3.10, generic library) and stays its own repo.

## Design

### Approach: Colocated Top-Level Packages

Move `zndraw_auth` and `zndraw_joblib` source and tests into the zndraw-fastapi repo as **top-level Python packages** alongside `zndraw`. All three ship in a single wheel under a single version.

This preserves:
- **All existing imports** — `from zndraw_auth import User` and `from zndraw_joblib import JobManager` work unchanged.
- **Test isolation** — each package's tests live in a separate directory with their own conftest and can run independently.
- **Separation of concerns** — auth and joblib remain self-contained packages; the repo boundary is replaced by a directory boundary.

### Repository Structure (After)

```
src/
├── zndraw/              # main app (unchanged)
├── zndraw_auth/         # from zndraw-auth/src/zndraw_auth/ (~8 modules, ~870 LOC)
└── zndraw_joblib/       # from zndraw-joblib/src/zndraw_joblib/ (~13 modules, ~3,600 LOC)
tests/
├── zndraw/              # current tests/ contents moved here
│   └── conftest.py      # full app fixtures (FastAPI + Redis + auth + joblib)
├── zndraw_auth/         # from zndraw-auth/tests/ (~6 test files)
│   └── conftest.py      # minimal FastAPI app with auth routers only
└── zndraw_joblib/       # from zndraw-joblib/tests/ (~28 files incl. conftest)
    └── conftest.py      # minimal FastAPI app with joblib router only
```

### pyproject.toml Changes

**Build targets** (preserve existing `artifacts` and `exclude` entries):
```toml
[tool.hatch.build.targets.wheel]
packages = ["src/zndraw", "src/zndraw_auth", "src/zndraw_joblib"]
```

**Dependencies — absorb from auth + joblib:**
- From zndraw-auth: `fastapi-users[sqlalchemy]>=14.0.0`, verify `aiosqlite>=0.19.0` is covered by existing `sqlalchemy[aiosqlite]` extra
- From zndraw-joblib: `taskiq>=0.12.1`, `taskiq-redis>=1.2.2` (already have httpx, pydantic, sqlmodel, zndraw-socketio)
- Remove: `zndraw-auth>=0.2.3`, `zndraw-joblib>=0.1.7`

**Ruff / isort** — add absorbed packages to first-party list:
```toml
[tool.ruff.lint.isort]
known-first-party = ["zndraw", "zndraw_auth", "zndraw_joblib"]
```

**Test configuration:**
```toml
[tool.pytest.ini_options]
testpaths = ["tests/zndraw", "tests/zndraw_auth", "tests/zndraw_joblib"]
asyncio_mode = "auto"
asyncio_default_fixture_loop_scope = "function"
```

### Testing Strategy

**All tests run against real services (Redis, etc.).** No mocks, no fakeredis in tests — ever.

**CI runs each suite sequentially to enforce isolation, plus a combined run:**
```yaml
- run: uv run pytest tests/zndraw_auth/
- run: uv run pytest tests/zndraw_joblib/
- run: uv run pytest tests/zndraw/
- run: uv run pytest tests/
```

The first three enforce isolation (each suite's conftest constructs a minimal FastAPI app with only its own routers). The final `pytest tests/` run covers everything together to catch cross-suite fixture collisions.

The zndraw_auth tests must not depend on zndraw or joblib fixtures. The zndraw_joblib tests must not depend on zndraw fixtures (joblib already depends on auth + socketio at the code level, which is fine).

### Migration Steps

1. **Copy source** — `zndraw-auth/src/zndraw_auth/` → `src/zndraw_auth/`, same for joblib.
2. **Copy tests** — `zndraw-auth/tests/` → `tests/zndraw_auth/`, same for joblib.
3. **Move existing tests** — `tests/*.py` → `tests/zndraw/`.
4. **Handle `_version.py`** — delete `_version.py` from auth and joblib. Update their `__init__.py` to re-export the version from zndraw: `from zndraw._version import __version__`. hatch-vcs only generates one version file (for zndraw); all three packages share a single version.
5. **Update pyproject.toml** — wheel packages, absorbed deps, test paths, `known-first-party`, `asyncio_default_fixture_loop_scope`.
6. **Lint compliance** — run `uvx prek --all-files` to fix formatting, linting, and import order.
7. **Update CI** — sequential per-suite runs + combined run.
8. **Fix test conftest paths** — adjust any relative imports or path assumptions in copied conftest files. Ensure all conftest files use real Redis.
9. **Ensure all tests pass** — `pytest tests/zndraw_auth/`, `pytest tests/zndraw_joblib/`, `pytest tests/zndraw/`, then `pytest tests/`.
10. **Remove editable-install references** — no more `../../zndraw-auth` paths in pyproject.toml.

### What Does NOT Change

- **All Python imports** — zero find-and-replace needed.
- **zndraw-socketio** — stays its own repo, its own PyPI package, its own CI.
- **Frontend** — completely unaffected.
- **Public API** — `zndraw`, `zndraw-cli`, `zndraw-db` entry points unchanged.

### Git History

No history preservation from auth/joblib repos. Clean copy with a single commit per migration step.

## Risks

| Risk | Mitigation |
|------|-----------|
| Accidental cross-package imports (auth importing from zndraw) | Per-suite CI enforcement; code review |
| conftest fixture collisions when running `pytest tests/` | Each suite uses isolated conftest; no shared fixtures at `tests/` root |
| Dependency conflicts from merging three dep lists | All three already share the same core deps (fastapi, sqlmodel, pydantic); conflicts unlikely |
| joblib's `zndraw_auth` import resolves to local copy, not installed package | This is the desired behavior — they're colocated, same wheel |
| PyPI shadowing — old `zndraw-auth`/`zndraw-joblib` packages remain on PyPI | Not a conflict: they are removed from zndraw's dependency list, so pip won't co-install them |
| Copied source fails stricter ruff config | Migration step 6 explicitly addresses this before tests run |

## Success Criteria

1. `pytest tests/zndraw_auth/` passes with 0 imports from zndraw or zndraw_joblib in auth source.
2. `pytest tests/zndraw_joblib/` passes with 0 imports from zndraw in joblib source.
3. `pytest tests/zndraw/` passes (all existing tests).
4. `uv build` produces a single wheel containing all three packages.
5. `pip install zndraw` (from built wheel) makes `import zndraw`, `import zndraw_auth`, `import zndraw_joblib` all work.
6. CI runs each suite sequentially and independently.
