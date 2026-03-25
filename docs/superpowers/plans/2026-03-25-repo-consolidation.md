# Repository Consolidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Consolidate zndraw-auth and zndraw-joblib into the zndraw-fastapi repo as colocated top-level packages with zero import changes.

**Architecture:** Copy source and tests from `../zndraw-auth` and `../zndraw-joblib` into `src/` and `tests/` respectively. Move existing zndraw tests into `tests/zndraw/`. Update `pyproject.toml` to build all three packages in a single wheel. Each test suite keeps its own conftest and runs independently.

**Tech Stack:** Python, hatchling, uv, pytest, GitHub Actions

**Spec:** `docs/superpowers/specs/2026-03-25-repo-consolidation-design.md`

---

## File Structure Map

**Create (copy from external repos):**
- `src/zndraw_auth/` — 7 files from `../zndraw-auth/src/zndraw_auth/` (excluding `_version.py`)
- `src/zndraw_joblib/` — 13 files from `../zndraw-joblib/src/zndraw_joblib/` (excluding `_version.py`)
- `tests/zndraw_auth/` — 5 test files + conftest from `../zndraw-auth/tests/`
- `tests/zndraw_joblib/` — 28 files + conftest from `../zndraw-joblib/tests/`

**Move (within this repo):**
- `tests/*.py`, `tests/test_cli_agent/`, `tests/test_client/`, `tests/worker/` → `tests/zndraw/`

**Modify:**
- `pyproject.toml` — wheel packages, deps, ruff isort, pytest config
- `.github/workflows/test.yaml` — sequential per-suite test runs

---

### Task 1: Copy zndraw-auth source

**Files:**
- Create: `src/zndraw_auth/__init__.py`, `admin.py`, `cli_login.py`, `db.py`, `schemas.py`, `settings.py`, `users.py`

- [ ] **Step 1: Copy auth source, excluding `_version.py` and `__pycache__`**

```bash
cp -r ../zndraw-auth/src/zndraw_auth/ src/zndraw_auth/
rm -f src/zndraw_auth/_version.py
rm -rf src/zndraw_auth/__pycache__
```

- [ ] **Step 2: Add version re-export to `__init__.py`**

Add this line to the top of `src/zndraw_auth/__init__.py` (after the docstring, before the first import):
```python
from zndraw._version import __version__ as __version__
```

- [ ] **Step 3: Verify the copy**

```bash
ls src/zndraw_auth/
```

Expected: `__init__.py  admin.py  cli_login.py  db.py  schemas.py  settings.py  users.py`

- [ ] **Step 4: Verify import works**

```bash
uv run python -c "from zndraw_auth import User; print('OK')"
```

Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_auth/
git commit -m "feat: copy zndraw-auth source into repo"
```

---

### Task 2: Copy zndraw-joblib source

**Files:**
- Create: `src/zndraw_joblib/__init__.py`, `client.py`, `dependencies.py`, `events.py`, `exceptions.py`, `models.py`, `provider.py`, `registry.py`, `router.py`, `schemas.py`, `settings.py`, `sweeper.py`

- [ ] **Step 1: Copy joblib source, excluding `_version.py` and `__pycache__`**

```bash
cp -r ../zndraw-joblib/src/zndraw_joblib/ src/zndraw_joblib/
rm -f src/zndraw_joblib/_version.py
rm -rf src/zndraw_joblib/__pycache__
```

- [ ] **Step 2: Add version re-export to `__init__.py`**

Add this line to the top of `src/zndraw_joblib/__init__.py` (after the docstring, before the first import):
```python
from zndraw._version import __version__ as __version__
```

- [ ] **Step 3: Verify the copy**

```bash
ls src/zndraw_joblib/
```

Expected: `__init__.py  client.py  dependencies.py  events.py  exceptions.py  models.py  provider.py  registry.py  router.py  schemas.py  settings.py  sweeper.py`

- [ ] **Step 4: Verify import works**

```bash
uv run python -c "from zndraw_joblib import JobManager; print('OK')"
```

Expected: `OK`

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_joblib/
git commit -m "feat: copy zndraw-joblib source into repo"
```

---

### Task 3: Move existing zndraw tests into `tests/zndraw/`

**Files:**
- Move: all `tests/*.py` files and `tests/test_cli_agent/`, `tests/test_client/`, `tests/worker/` → `tests/zndraw/`

- [ ] **Step 1: Create target directory and move all test files**

```bash
mkdir -p tests/zndraw
# Move all .py files (including conftest.py)
mv tests/*.py tests/zndraw/
# Move subdirectories
mv tests/test_cli_agent tests/zndraw/
mv tests/test_client tests/zndraw/
mv tests/worker tests/zndraw/
# Clean stale bytecode
rm -rf tests/__pycache__
```

- [ ] **Step 2: Verify the move**

```bash
ls tests/zndraw/conftest.py tests/zndraw/test_admin.py
ls tests/zndraw/test_cli_agent/ tests/zndraw/test_client/ tests/zndraw/worker/
```

Expected: all files present, `tests/` root has no `.py` files left.

- [ ] **Step 3: Verify tests still discover correctly**

```bash
uv run pytest tests/zndraw/ --collect-only -q 2>&1 | tail -5
```

Expected: shows collected test count (no import errors).

- [ ] **Step 4: Commit**

```bash
git add tests/
git commit -m "refactor: move zndraw tests into tests/zndraw/"
```

---

### Task 4: Copy zndraw-auth tests

**Files:**
- Create: `tests/zndraw_auth/conftest.py`, `test_admin_token.py`, `test_auth.py`, `test_cli_login.py`, `test_scoped_session.py`, `test_users_router.py`

- [ ] **Step 1: Copy auth tests, excluding `__pycache__`**

```bash
cp -r ../zndraw-auth/tests/ tests/zndraw_auth/
find tests/zndraw_auth -name __pycache__ -exec rm -rf {} + 2>/dev/null; true
```

- [ ] **Step 2: Verify the copy**

```bash
ls tests/zndraw_auth/
```

Expected: `conftest.py  test_admin_token.py  test_auth.py  test_cli_login.py  test_scoped_session.py  test_users_router.py`

- [ ] **Step 3: Verify test collection**

```bash
uv run pytest tests/zndraw_auth/ --collect-only -q 2>&1 | tail -5
```

Expected: shows collected tests, no import errors.

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw_auth/
git commit -m "feat: copy zndraw-auth tests into repo"
```

---

### Task 5: Copy zndraw-joblib tests

**Files:**
- Create: `tests/zndraw_joblib/` — all 28+ files from `../zndraw-joblib/tests/`

- [ ] **Step 1: Copy joblib tests, excluding `__pycache__`**

```bash
cp -r ../zndraw-joblib/tests/ tests/zndraw_joblib/
find tests/zndraw_joblib -name __pycache__ -exec rm -rf {} + 2>/dev/null; true
```

- [ ] **Step 2: Verify the copy**

```bash
ls tests/zndraw_joblib/ | wc -l
```

Expected: ~29 files (28 test/support files + conftest).

- [ ] **Step 3: Verify test collection**

```bash
uv run pytest tests/zndraw_joblib/ --collect-only -q 2>&1 | tail -5
```

Expected: shows collected tests, no import errors.

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw_joblib/
git commit -m "feat: copy zndraw-joblib tests into repo"
```

---

### Task 6: Update pyproject.toml

**Files:**
- Modify: `pyproject.toml`

This task makes 5 changes to `pyproject.toml`:

- [ ] **Step 1: Update wheel packages**

Change:
```toml
packages = ["src/zndraw"]
```
To:
```toml
packages = ["src/zndraw", "src/zndraw_auth", "src/zndraw_joblib"]
```

- [ ] **Step 2: Update dependencies — remove old packages, add absorbed deps**

Remove these two lines from `[project] dependencies`:
```
"zndraw-auth>=0.2.3",
"zndraw-joblib>=0.1.7",
```

Add these new dependencies (absorbed from auth and joblib):
```
"fastapi-users[sqlalchemy]>=14.0.0",
"taskiq>=0.12.1",
"taskiq-redis>=1.2.2",
```

Verify `aiosqlite` is covered: the existing `sqlalchemy[asyncio, aiosqlite]>=2.0.46` installs aiosqlite as an extra. Confirm with:
```bash
uv run python -c "import aiosqlite; print(aiosqlite.__version__)"
```

- [ ] **Step 3: Add `known-first-party` for isort**

Change:
```toml
known-first-party = ["zndraw"]
```
To:
```toml
known-first-party = ["zndraw", "zndraw_auth", "zndraw_joblib"]
```

- [ ] **Step 4: Add pytest configuration section**

Add before `[tool.codespell]`:
```toml
[tool.pytest.ini_options]
testpaths = ["tests/zndraw", "tests/zndraw_auth", "tests/zndraw_joblib"]
asyncio_mode = "auto"
asyncio_default_fixture_loop_scope = "function"
```

- [ ] **Step 5: Verify no editable-install references remain**

Check that `[tool.uv.sources]` in `pyproject.toml` has no relative paths to `../zndraw-auth` or `../zndraw-joblib`. Currently it is empty — confirm this is still the case.

- [ ] **Step 6: Run `uv sync` to resolve new dependency graph**

```bash
uv sync --all-extras --dev
```

Expected: resolves without conflicts.

- [ ] **Step 7: Verify build produces correct wheel**

```bash
uv build --wheel
uv run python -c "
import zipfile, glob
whl = glob.glob('dist/*.whl')[-1]
with zipfile.ZipFile(whl) as z:
    dirs = {n.split('/')[0] for n in z.namelist()}
    for pkg in ['zndraw', 'zndraw_auth', 'zndraw_joblib']:
        assert pkg in dirs, f'{pkg} missing from wheel'
    print('All three packages in wheel: OK')
"
```

Expected: `All three packages in wheel: OK`

- [ ] **Step 8: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "build: absorb auth+joblib deps, configure wheel for 3 packages"
```

---

### Task 7: Lint compliance

**Files:**
- Modify: potentially any file in `src/zndraw_auth/`, `src/zndraw_joblib/`, `tests/zndraw_auth/`, `tests/zndraw_joblib/`

- [ ] **Step 1: Run pre-commit on all files**

```bash
uvx prek --all-files
```

This applies ruff format + ruff check (with fixes) + import sorting across the entire repo, including the newly copied files which may not comply with zndraw's stricter ruff config.

- [ ] **Step 2: Review and commit any changes**

```bash
git diff --stat
git add src/zndraw_auth/ src/zndraw_joblib/ tests/zndraw_auth/ tests/zndraw_joblib/ tests/zndraw/
git commit -m "style: lint compliance for absorbed auth+joblib packages"
```

If no changes, skip this commit.

---

### Task 8: Update CI workflow

**Files:**
- Modify: `.github/workflows/test.yaml`

- [ ] **Step 1: Replace single pytest step with sequential per-suite runs**

Change the Pytest step from:
```yaml
      - name: Pytest
        run: |
          uv run python --version
          uv run pytest --cov --junitxml=junit.xml -o junit_family=legacy
```

To:
```yaml
      - name: Pytest (zndraw-auth)
        run: uv run pytest tests/zndraw_auth/ -v
      - name: Pytest (zndraw-joblib)
        run: uv run pytest tests/zndraw_joblib/ -v
      - name: Pytest (zndraw)
        run: uv run pytest tests/zndraw/ -v
      - name: Pytest (all — combined)
        run: uv run pytest --cov --junitxml=junit.xml -o junit_family=legacy
```

- [ ] **Step 2: Commit**

```bash
git add .github/workflows/test.yaml
git commit -m "ci: run auth/joblib/zndraw test suites sequentially"
```

---

### Task 9: Verify all tests pass

**Files:** None (verification only)

- [ ] **Step 1: Run auth tests in isolation**

```bash
uv run pytest tests/zndraw_auth/ -v
```

Expected: all auth tests pass.

- [ ] **Step 2: Run joblib tests in isolation**

```bash
uv run pytest tests/zndraw_joblib/ -v
```

Expected: all joblib tests pass.

- [ ] **Step 3: Run zndraw tests in isolation**

```bash
uv run pytest tests/zndraw/ -v
```

Expected: all zndraw tests pass.

- [ ] **Step 4: Run all tests combined**

```bash
uv run pytest tests/ -v
```

Expected: all tests pass, no fixture collisions.

- [ ] **Step 5: Verify pyright type checking**

```bash
uv run pyright src/zndraw_auth/ src/zndraw_joblib/
```

Expected: no new errors (existing baseline errors from redis stubs etc. are acceptable).
