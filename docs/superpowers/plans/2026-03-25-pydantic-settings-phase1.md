# Pydantic-Settings Phase 1: Server-Side Unification

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Change `Settings` env prefix from `ZNDRAW_` to `ZNDRAW_SERVER_`, add `pyproject.toml` support to all three settings classes, and adopt the thin-CLI pattern (Typer defaults to `None`, no `envvar=`, no `os.environ` writes).

**Architecture:** Each settings class (`Settings`, `AuthSettings`, `JobLibSettings`) independently adds `PyprojectTomlConfigSettingsSource` via `settings_customise_sources()`. The CLI uses the thin-Typer pattern where options default to `None` and non-None values are passed as `Settings(**overrides)`. Priority: init args > env vars > pyproject.toml > field defaults.

**Tech Stack:** pydantic-settings v2 (`PyprojectTomlConfigSettingsSource`, `settings_customise_sources`), Typer

**Spec:** `docs/superpowers/specs/2026-03-25-pydantic-settings-unification-design.md`

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `src/zndraw/config.py` | Modify | Change `env_prefix` to `ZNDRAW_SERVER_`, drop `env_nested_delimiter` (no longer needed — no nested fields), rename `server_url` → `internal_url`, add `settings_customise_sources()` |
| `src/zndraw/broker.py` | Modify | Update `settings.server_url` → `settings.internal_url`, update docstring/error msg |
| `src/zndraw/cli.py` | Modify | Remove `envvar=`, remove `os.environ` writes, thin-CLI pattern. Change `host` default from `"127.0.0.1"` to `None` (**behavioral change**: Settings default `"0.0.0.0"` takes over) |
| `zndraw-auth: settings.py` | Modify | Add `settings_customise_sources()` |
| `tests/conftest.py` | Modify | Rename `ZNDRAW_` server env vars to `ZNDRAW_SERVER_*`. Fix orphaned `ZNDRAW_PRESENCE_TTL` docstring. |
| `tests/test_config.py` | Modify | Rename env vars, add pyproject.toml tests |
| `tests/test_lifespan.py` | Modify | Rename env vars |
| `tests/test_socketio_scaling.py` | Modify | Rename env vars |
| `tests/test_cli.py` | Modify | Rewrite/remove tests that assert `os.environ` writes (old pattern) |
| `tests/test_template_room_isolation.py` | Modify | Rename env vars |
| `tests/worker/test_resilience.py` | Modify | Rename env vars |
| `docker/standalone/docker-compose.yaml` | Modify | Rename env vars |
| `docker/production/docker-compose.yaml` | Modify | Rename env vars |
| `docker/templates/.env` | Modify | Rename env vars |
| `docker/standalone/README.md` | Modify | Rename env vars in reference tables |
| `docker/production/README.md` | Modify | Rename env vars in reference tables |
| `README.md` | Modify | Update examples, add env var / pyproject.toml reference |

**Note on `server_url` → `internal_url` rename:** Only rename `settings.server_url` attribute accesses (pydantic field). Do NOT rename local variables or function parameters named `server_url` in `cli.py` or other files — those are unrelated.

**Note on `host` default behavioral change:** The CLI currently defaults `host` to `"127.0.0.1"` (localhost-only). After thin-CLI adoption, Typer defaults to `None`, and `Settings.host` default `"0.0.0.0"` (all interfaces) takes over. This is intentional — the Settings class is the single source of truth for defaults.

**Note on `env_nested_delimiter`:** The current `Settings` has `env_nested_delimiter="__"`. This is deliberately removed — `Settings` has no nested fields, and removing it avoids accidental collision with future `ClientSettings` which also uses `ZNDRAW_` prefix (Phase 2).

---

### Task 1: Update `Settings` class — prefix + pyproject.toml + field rename

**Files:**
- Modify: `src/zndraw/config.py`

- [ ] **Step 1: Read the current file**

Read `src/zndraw/config.py` to confirm exact current state.

- [ ] **Step 2: Change `env_prefix`, drop `env_nested_delimiter`, add `settings_customise_sources`, rename field**

Replace the full `Settings` class. Key changes:
- `env_prefix="ZNDRAW_"` → `env_prefix="ZNDRAW_SERVER_"`
- Remove `env_nested_delimiter="__"` (no nested fields)
- Add `pyproject_toml_table_header=("tool", "zndraw", "server")`
- Add `settings_customise_sources()` classmethod
- Rename `server_url` → `internal_url`
- Update docstring: `ZNDRAW_` → `ZNDRAW_SERVER_`

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/config.py
git commit -m "feat: change Settings env_prefix to ZNDRAW_SERVER_, add pyproject.toml source

Breaking change: all server env vars now use ZNDRAW_SERVER_ prefix.
Drop env_nested_delimiter (no nested fields).
Rename server_url → internal_url to avoid ZNDRAW_SERVER_SERVER_URL stutter.
Add PyprojectTomlConfigSettingsSource for [tool.zndraw.server] support."
```

---

### Task 2: Update `broker.py` — field rename + docstring

**Files:**
- Modify: `src/zndraw/broker.py`

- [ ] **Step 1: Read the current file**

Read `src/zndraw/broker.py` to find all `server_url` and `ZNDRAW_REDIS_URL` references.

- [ ] **Step 2: Update references**

Changes needed:
- Docstring: `ZNDRAW_REDIS_URL` → `ZNDRAW_SERVER_REDIS_URL`
- Error message: `ZNDRAW_REDIS_URL` → `ZNDRAW_SERVER_REDIS_URL`
- Field access: `settings.server_url` → `settings.internal_url`

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/broker.py
git commit -m "refactor: update broker.py for ZNDRAW_SERVER_ prefix and internal_url rename"
```

---

### Task 3: Rename `settings.server_url` attribute accesses across codebase

**Files:**
- Modify: any file referencing `settings.server_url` or `Settings.server_url`

- [ ] **Step 1: Search for `settings.server_url` attribute accesses only**

```bash
uv run ruff check . 2>&1 | head -20  # quick check
```

Then search specifically for the Settings field usage:

```bash
rg "settings\.server_url|\.server_url" src/zndraw/ --type py
```

**Important:** Only rename `settings.server_url` (attribute access on Settings instance). Do NOT rename local variables or function parameters named `server_url` in `cli.py` or elsewhere — those are unrelated to the Settings field.

- [ ] **Step 2: Update each `settings.server_url` → `settings.internal_url`**

- [ ] **Step 3: Run pyright on changed files**

```bash
uv run pyright src/zndraw/config.py src/zndraw/broker.py
```

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/broker.py src/zndraw/config.py
# add any other changed files by name
git commit -m "refactor: rename settings.server_url → settings.internal_url across codebase"
```

---

### Task 4: Rename env vars in test infrastructure

**Files:**
- Modify: `tests/conftest.py`

- [ ] **Step 1: Read conftest.py**

Read `tests/conftest.py` and identify every `ZNDRAW_` env var that maps to a `Settings` field (not `ZNDRAW_AUTH_*` — those stay).

- [ ] **Step 2: Rename env vars**

Apply these renames throughout `tests/conftest.py`:
- `ZNDRAW_DATABASE_URL` → `ZNDRAW_SERVER_DATABASE_URL`
- `ZNDRAW_REDIS_URL` → `ZNDRAW_SERVER_REDIS_URL`
- `ZNDRAW_HOST` → `ZNDRAW_SERVER_HOST`
- `ZNDRAW_PORT` → `ZNDRAW_SERVER_PORT`
- `ZNDRAW_STORAGE` → `ZNDRAW_SERVER_STORAGE` (if present)

Also fix the orphaned docstring at line ~328: `ZNDRAW_PRESENCE_TTL` is not a Settings field. `presence_ttl` does not exist on `Settings`. Replace with a valid field in the docstring example, e.g., `ZNDRAW_SERVER_EDIT_LOCK_TTL`.

Leave `ZNDRAW_AUTH_*` env vars unchanged.

- [ ] **Step 3: Commit**

```bash
git add tests/conftest.py
git commit -m "test: rename ZNDRAW_ → ZNDRAW_SERVER_ env vars in conftest.py"
```

---

### Task 5: Rename env vars in remaining test files

**Files:**
- Modify: `tests/test_config.py`, `tests/test_lifespan.py`, `tests/test_socketio_scaling.py`, `tests/test_template_room_isolation.py`, `tests/worker/test_resilience.py`

- [ ] **Step 1: Read each file and apply renames**

For each file, rename `ZNDRAW_` server env vars to `ZNDRAW_SERVER_*`. The mapping is the same as Task 4. Also rename `ZNDRAW_SERVER_URL` → `ZNDRAW_SERVER_INTERNAL_URL` where it refers to the Settings field.

**Important:** `test_config.py` tests specifically test env var reading. Update both the env var names AND the assertion comments/docstrings. The existing tests use class-based organization (`class TestStorageUri:` etc.) — the project CLAUDE.md prefers bare functions but don't refactor existing test structure, just update the env var names.

- [ ] **Step 2: Run the affected tests**

```bash
uv run pytest tests/test_config.py tests/test_lifespan.py tests/test_socketio_scaling.py -x -v
```

Expected: All pass with new env var names.

- [ ] **Step 3: Commit**

```bash
git add tests/test_config.py tests/test_lifespan.py tests/test_socketio_scaling.py tests/test_template_room_isolation.py tests/worker/test_resilience.py
git commit -m "test: rename ZNDRAW_ → ZNDRAW_SERVER_ env vars in test files"
```

---

### Task 6: Run full test suite to verify prefix change

- [ ] **Step 1: Run all tests**

```bash
uv run pytest tests/ -x --timeout=900
```

Expected: All tests pass. If any fail, the failure will be due to a missed env var rename — fix and re-run.

- [ ] **Step 2: Run pyright**

```bash
uv run pyright .
```

- [ ] **Step 3: Run ruff**

```bash
uv run ruff check .
uv run ruff format --check .
```

---

### Task 7: Update `cli.py` — thin-CLI pattern

This task makes significant behavioral changes. Tests that assert `os.environ` writes will break and need rewriting.

**Files:**
- Modify: `src/zndraw/cli.py`
- Modify: `tests/test_cli.py`

- [ ] **Step 1: Read cli.py and test_cli.py**

Read `src/zndraw/cli.py` focusing on:
- Lines ~466-472: `--port` (with `envvar="ZNDRAW_PORT"`) and `--host` (with `envvar="ZNDRAW_HOST"`, default `"127.0.0.1"`)
- Lines ~340-341: `os.environ["ZNDRAW_HOST"]` and `os.environ["ZNDRAW_PORT"]` writes
- Line ~72: `os.environ["ZNDRAW_DATABASE_URL"]` write
- Line ~55: help text referencing `ZNDRAW_DATABASE_URL`

Read `tests/test_cli.py` focusing on tests that assert env var writes:
- `test_cli_sets_host_and_port_env_vars` (line ~245)
- `test_cli_writes_default_port_to_env` (line ~270)
- `test_cli_reads_host_from_env` (line ~290)

- [ ] **Step 2: Remove `envvar=` from Typer options, change defaults to `None`**

```python
# Before:
port: Annotated[int | None, typer.Option("--port", ..., envvar="ZNDRAW_PORT")] = None,
host: Annotated[str, typer.Option(help="...", envvar="ZNDRAW_HOST")] = "127.0.0.1",

# After:
port: Annotated[int | None, typer.Option("--port",
    help="Server port [env: ZNDRAW_SERVER_PORT].")] = None,
host: Annotated[str | None, typer.Option(
    help="Server hostname or IP address [env: ZNDRAW_SERVER_HOST].")] = None,
```

**Behavioral change:** `host` default changes from `"127.0.0.1"` to `None`. When `None`, `Settings.host` default `"0.0.0.0"` takes over. This is intentional — Settings is the single source of truth.

- [ ] **Step 3: Remove `os.environ` writes in `resolve_server()`**

Remove `os.environ["ZNDRAW_HOST"] = host` and `os.environ["ZNDRAW_PORT"] = str(port)`. Instead, build `Settings(**overrides)` and pass it through to the lifespan (via `app.state.settings`).

- [ ] **Step 4: Update `db_init()` similarly**

Remove `os.environ["ZNDRAW_DATABASE_URL"] = database_url`. Instead:

```python
overrides = {}
if database_url:
    overrides["database_url"] = database_url
settings = Settings(**overrides)
```

- [ ] **Step 5: Update help text referencing old env var names**

```python
# Before:
help="Database URL (overrides ZNDRAW_DATABASE_URL)"
# After:
help="Database URL [env: ZNDRAW_SERVER_DATABASE_URL]."
```

- [ ] **Step 6: Rewrite `test_cli.py` tests that assert `os.environ` writes**

The following tests assert the old pattern (CLI writes env vars) and must be rewritten or removed:
- `test_cli_sets_host_and_port_env_vars` — rewrite to test that `Settings(**overrides)` receives the CLI values
- `test_cli_writes_default_port_to_env` — rewrite or remove (default now comes from Settings, not env)
- `test_cli_reads_host_from_env` — rewrite to test `Settings` reads `ZNDRAW_SERVER_HOST` from env

- [ ] **Step 7: Run CLI-related tests**

```bash
uv run pytest tests/test_cli.py -x -v
```

- [ ] **Step 8: Commit**

```bash
git add src/zndraw/cli.py tests/test_cli.py
git commit -m "refactor: adopt thin-CLI pattern — remove envvar= and os.environ writes

Typer options now default to None with no envvar= parameter.
Pydantic-settings is the single source of truth for defaults and env vars.
Behavioral change: host default is now 0.0.0.0 (from Settings) instead of 127.0.0.1."
```

---

### Task 8: Update Docker files

**Files:**
- Modify: `docker/standalone/docker-compose.yaml`
- Modify: `docker/production/docker-compose.yaml`
- Modify: `docker/templates/.env`
- Modify: `docker/standalone/README.md`
- Modify: `docker/production/README.md`

- [ ] **Step 1: Read all five files**

- [ ] **Step 2: Rename `ZNDRAW_` server env vars to `ZNDRAW_SERVER_*`**

Apply renames across all five files:
- `ZNDRAW_DATABASE_URL` → `ZNDRAW_SERVER_DATABASE_URL`
- `ZNDRAW_STORAGE` → `ZNDRAW_SERVER_STORAGE`
- `ZNDRAW_INIT_DB_ON_STARTUP` → `ZNDRAW_SERVER_INIT_DB_ON_STARTUP`
- `ZNDRAW_SERVER_URL` → `ZNDRAW_SERVER_INTERNAL_URL` (note: the old env var `ZNDRAW_SERVER_URL` mapped to field `server_url` with prefix `ZNDRAW_`. After prefix change to `ZNDRAW_SERVER_` and field rename to `internal_url`, the new env var is `ZNDRAW_SERVER_INTERNAL_URL`)
- `ZNDRAW_REDIS_URL` → `ZNDRAW_SERVER_REDIS_URL`
- `ZNDRAW_WORKER_ENABLED` → `ZNDRAW_SERVER_WORKER_ENABLED`
- `ZNDRAW_GUEST_PASSWORD` → `ZNDRAW_SERVER_GUEST_PASSWORD`
- `ZNDRAW_WORKER_PASSWORD` → `ZNDRAW_SERVER_WORKER_PASSWORD`

Leave unchanged: `ZNDRAW_AUTH_*`, `ZNDRAW_URL`, `ZNDRAW_USER`, `ZNDRAW_PASSWORD` (client/auth env vars).

- [ ] **Step 3: Commit**

```bash
git add docker/standalone/docker-compose.yaml docker/production/docker-compose.yaml docker/templates/.env docker/standalone/README.md docker/production/README.md
git commit -m "chore: rename ZNDRAW_ → ZNDRAW_SERVER_ in Docker files"
```

---

### Task 9: Add pyproject.toml source tests

**Files:**
- Modify: `tests/test_config.py`

- [ ] **Step 1: Write tests for pyproject.toml loading**

Add bare test functions (following CLAUDE.md convention — existing class-based tests in this file are legacy):

```python
def test_settings_from_pyproject_toml(tmp_path, monkeypatch):
    """Settings should load from [tool.zndraw.server] in pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
host = "192.168.1.1"
storage = "/data/frames.lmdb"
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)

    settings = Settings()
    assert settings.port == 9999
    assert settings.host == "192.168.1.1"
    assert settings.storage == "/data/frames.lmdb"


def test_env_overrides_pyproject_toml(tmp_path, monkeypatch):
    """Env vars should take priority over pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "7777")

    settings = Settings()
    assert settings.port == 7777


def test_init_overrides_env_and_pyproject(tmp_path, monkeypatch):
    """Init args should take priority over env and pyproject.toml."""
    toml_content = """\
[tool.zndraw.server]
port = 9999
"""
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text(toml_content)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("ZNDRAW_SERVER_PORT", "7777")

    settings = Settings(port=5555)
    assert settings.port == 5555


def test_missing_pyproject_toml_is_silent(tmp_path, monkeypatch):
    """Settings should work fine without a pyproject.toml."""
    monkeypatch.chdir(tmp_path)
    settings = Settings()
    assert settings.port == 8000  # default
```

- [ ] **Step 2: Run the new tests**

```bash
uv run pytest tests/test_config.py -x -v
```

Expected: All pass.

- [ ] **Step 3: Commit**

```bash
git add tests/test_config.py
git commit -m "test: add pyproject.toml source tests for Settings"
```

---

### Task 10: Add pyproject.toml source to `AuthSettings` (zndraw-auth repo)

**Files:**
- Modify: `/Users/fzills/tools/zndraw-auth/src/zndraw_auth/settings.py`

- [ ] **Step 1: Read current AuthSettings**

- [ ] **Step 2: Add `pyproject_toml_table_header` and `settings_customise_sources`**

Add to `model_config`:
- `pyproject_toml_table_header=("tool", "zndraw", "auth")`

Add classmethod `settings_customise_sources()` returning `(init_settings, env_settings, dotenv_settings, PyprojectTomlConfigSettingsSource(settings_cls))`. Keep `dotenv_settings` since AuthSettings already supports `.env` files.

- [ ] **Step 3: Add test in zndraw-auth repo**

- [ ] **Step 4: Commit in zndraw-auth repo**

```bash
cd /Users/fzills/tools/zndraw-auth
git add src/zndraw_auth/settings.py
git commit -m "feat: add pyproject.toml source for [tool.zndraw.auth]"
```

---

### Task 11: Add pyproject.toml source to `JobLibSettings` (zndraw-joblib repo)

**Files:**
- Modify: `zndraw-joblib` package's `settings.py` (in its own repo — not in zndraw-fastapi)

- [ ] **Step 1: Locate and read JobLibSettings source**

Check if zndraw-joblib is available as a local editable install or find its repo. Read the `settings.py` file.

- [ ] **Step 2: Add `pyproject_toml_table_header` and `settings_customise_sources`**

Same pattern as Task 10 but with:
- `pyproject_toml_table_header=("tool", "zndraw", "joblib")`
- Source chain: `(init_settings, env_settings, PyprojectTomlConfigSettingsSource(settings_cls))` (no dotenv)

- [ ] **Step 3: Add test and commit in zndraw-joblib repo**

---

### Task 12: Update README

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Read current README**

- [ ] **Step 2: Update Python API examples**

The current README shows `vis = ZnDraw(url="http://localhost:1234", room="my-room")`. Add a zero-config example and mention env var / pyproject.toml support. Note: `ZnDraw()` zero-arg auto-discovery is Phase 2 — for now, show the env var / pyproject.toml support.

- [ ] **Step 3: Update Self-Hosting section with new env var names**

If the Self-Hosting section references any `ZNDRAW_` server env vars, update them to `ZNDRAW_SERVER_*`.

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: update README for pydantic-settings unification (Phase 1)"
```

---

### Task 13: Final verification

- [ ] **Step 1: Run full test suite**

```bash
uv run pytest tests/ --timeout=900
```

- [ ] **Step 2: Run ruff format + check**

```bash
uv run ruff format .
uv run ruff check --select I --fix .
```

- [ ] **Step 3: Run pyright**

```bash
uv run pyright .
```

- [ ] **Step 4: Verify Docker configs parse**

```bash
cd docker/standalone && docker compose config > /dev/null
cd docker/production && docker compose config > /dev/null
```

- [ ] **Step 5: Commit any formatting fixes**

```bash
git add src/ tests/
git commit -m "style: format after Phase 1 settings unification"
```
