# CLI Auth Consistency Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Unify CLI authentication by making `--user`/`--password` explicit Typer flags and extracting `resolve_token()` into a shared utility.

**Architecture:** New `src/zndraw/auth_utils.py` module owns token resolution. CLI (`connection.py`) and Python client (`client/core.py`) both import it. Every CLI command gets `--user`/`--password` flags alongside existing `--token`.

**Tech Stack:** Python, Typer, httpx, Pydantic (`SecretStr`), pytest

**Spec:** `docs/superpowers/specs/2026-03-19-cli-auth-consistency-design.md`

---

## File Map

| File | Role | Action |
|---|---|---|
| `src/zndraw/auth_utils.py` | Shared token resolution (single source of truth) | **Create** |
| `tests/test_resolve_token.py` | Tests for shared `resolve_token()` | **Rewrite** |
| `src/zndraw/cli_agent/connection.py` | CLI connection helpers, option types | **Modify** — remove old `resolve_token()`, add `UserOpt`/`PasswordOpt`, update `get_connection()`/`get_zndraw()` |
| `src/zndraw/cli_agent/auth.py` | Auth subcommands | **Modify** — add `--user`/`--password` to `auth status` |
| `src/zndraw/cli_agent/rooms.py` | Room commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/frames.py` | Frame commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/sessions.py` | Session commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/step.py` | Step commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/selection.py` | Selection commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/selection_groups.py` | Selection group commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/bookmarks.py` | Bookmark commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/figures.py` | Figure commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/screenshots.py` | Screenshot commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/extensions.py` | Extension commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/chat.py` | Chat commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/geometries.py` | Geometry commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/gif.py` | GIF capture command | **Modify** — add auth params |
| `src/zndraw/cli_agent/presets.py` | Preset commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/jobs.py` | Job commands | **Modify** — add auth params |
| `src/zndraw/cli_agent/mount.py` | Mount command | **Modify** — add auth params |
| `src/zndraw/cli_agent/admin.py` | Admin commands | **Modify** — add auth params |
| `src/zndraw/client/core.py` | ZnDraw client class | **Modify** — use shared `resolve_token()` |
| `docker/standalone/docker-compose.yaml` | Docker standalone config | **Modify** — rename env var |
| `docker/production/docker-compose.yaml` | Docker production config | **Modify** — rename env var |
| `skills/zndraw/SKILL.md` | Skill documentation | **Modify** — update env var references |

---

## Task 1: Create shared `resolve_token()` with tests (TDD)

**Files:**
- Create: `src/zndraw/auth_utils.py`
- Rewrite: `tests/test_resolve_token.py`

### Step 1.1: Write validation tests

- [ ] **Write failing tests for input validation**

```python
# tests/test_resolve_token.py
"""Tests for shared resolve_token utility."""

from __future__ import annotations

from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import httpx
import pytest
from pydantic import SecretStr

from zndraw.server_manager import TokenEntry, TokenStore


@pytest.fixture
def token_store(tmp_path):
    return TokenStore(directory=tmp_path)


@pytest.fixture
def stored_entry():
    return TokenEntry(
        access_token="stored.jwt.token",
        email="user@example.com",
        stored_at=datetime(2026, 3, 1, tzinfo=UTC),
    )


@pytest.fixture
def mock_httpx_client():
    """Yield a mock for httpx.Client as context manager."""
    mock_client = MagicMock()
    mock_client.__enter__ = MagicMock(return_value=mock_client)
    mock_client.__exit__ = MagicMock(return_value=False)
    with patch("zndraw.auth_utils.httpx.Client", return_value=mock_client):
        yield mock_client


# --- Validation tests ---


def test_token_and_user_raises():
    """Cannot combine --token with --user/--password."""
    from zndraw.auth_utils import resolve_token

    with pytest.raises(ValueError, match="Cannot combine"):
        resolve_token("http://localhost:8000", token="t", user="u", password="p")


def test_token_and_password_raises():
    """Cannot combine --token with --password alone."""
    from zndraw.auth_utils import resolve_token

    with pytest.raises(ValueError, match="Cannot combine"):
        resolve_token("http://localhost:8000", token="t", password="p")


def test_user_without_password_raises():
    """--user without --password should fail."""
    from zndraw.auth_utils import resolve_token

    with pytest.raises(ValueError, match="Missing --password"):
        resolve_token("http://localhost:8000", user="u")


def test_password_without_user_raises():
    """--password without --user should fail."""
    from zndraw.auth_utils import resolve_token

    with pytest.raises(ValueError, match="Missing --user"):
        resolve_token("http://localhost:8000", password="p")


def test_secretstr_password_accepted(mock_httpx_client):
    """SecretStr password should be unwrapped automatically."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.token"}
    mock_httpx_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.get_token_store"):
        result = resolve_token(
            "http://localhost:8000",
            user="admin@example.com",
            password=SecretStr("secret"),
        )

    assert result == "login.token"
    mock_httpx_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "admin@example.com", "password": "secret"},
    )
```

- [ ] **Run tests to verify they fail**

Run: `uv run pytest tests/test_resolve_token.py -v -x`
Expected: FAIL — `zndraw.auth_utils` does not exist yet

### Step 1.2: Implement `resolve_token()` — validation and explicit credentials

- [ ] **Create `src/zndraw/auth_utils.py` with validation + tier 1**

```python
# src/zndraw/auth_utils.py
"""Shared token resolution for CLI and Python client.

Single source of truth for the authentication fallback chain:
1. Explicit credentials (--token OR --user/--password)
2. Stored token from ~/.zndraw/tokens.json
3. Guest session fallback
"""

from __future__ import annotations

import os
import sys
import warnings

import httpx
from pydantic import SecretStr

from zndraw.server_manager import TokenStore


def get_token_store() -> TokenStore:
    """Return the default TokenStore (testable seam)."""
    return TokenStore()


def resolve_token(
    base_url: str,
    token: str | None = None,
    user: str | None = None,
    password: SecretStr | str | None = None,
) -> str:
    """Resolve an auth token from explicit credentials, stored token, or guest.

    Resolution order (explicit always wins over implicit):

    1. Explicit: ``token`` returned as-is, OR ``user``+``password``
       login via ``POST /v1/auth/jwt/login`` (token NOT stored).
    2. Stored token from ``~/.zndraw/tokens.json`` (validated, deleted on 401).
    3. Guest fallback via ``POST /v1/auth/guest``.

    Parameters
    ----------
    base_url
        Server URL for network calls.
    token
        Explicit JWT token (from ``--token`` / ``ZNDRAW_TOKEN``).
    user
        User email (from ``--user`` / ``ZNDRAW_USER``).
    password
        Password (from ``--password`` / ``ZNDRAW_PASSWORD``).
        Accepts ``SecretStr`` or plain ``str``.

    Raises
    ------
    ValueError
        On invalid credential combinations.
    httpx.HTTPStatusError
        On failed login or guest session creation.
    """
    # --- Migration warning ---
    if os.environ.get("ZNDRAW_EMAIL") and not os.environ.get("ZNDRAW_USER"):
        warnings.warn(
            "ZNDRAW_EMAIL is deprecated, use ZNDRAW_USER instead.",
            DeprecationWarning,
            stacklevel=2,
        )

    # --- Validation (fail fast) ---
    if token is not None and (user is not None or password is not None):
        msg = "Cannot combine --token with --user/--password"
        raise ValueError(msg)
    if user is not None and password is None:
        msg = "Missing --password (required when --user is provided)"
        raise ValueError(msg)
    if password is not None and user is None:
        msg = "Missing --user (required when --password is provided)"
        raise ValueError(msg)

    # Unwrap SecretStr
    raw_password: str | None = None
    if isinstance(password, SecretStr):
        raw_password = password.get_secret_value()
    elif isinstance(password, str):
        raw_password = password

    # --- Tier 1: Explicit credentials ---
    if token is not None:
        return token

    if user is not None and raw_password is not None:
        with httpx.Client(base_url=base_url, timeout=10.0) as client:
            resp = client.post(
                "/v1/auth/jwt/login",
                data={"username": user, "password": raw_password},
            )
            resp.raise_for_status()
            return resp.json()["access_token"]

    # --- Tier 2: Stored token ---
    store = get_token_store()
    entry = store.get(base_url)
    if entry is not None:
        with httpx.Client(
            base_url=base_url,
            headers={"Authorization": f"Bearer {entry.access_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get("/v1/auth/users/me")
            if resp.status_code == 200:
                return entry.access_token
            store.delete(base_url)

    # --- Tier 3: Guest fallback ---
    with httpx.Client(base_url=base_url, timeout=10.0) as client:
        resp = client.post("/v1/auth/guest")
        resp.raise_for_status()
        return resp.json()["access_token"]
```

- [ ] **Run validation tests to verify they pass**

Run: `uv run pytest tests/test_resolve_token.py::test_token_and_user_raises tests/test_resolve_token.py::test_token_and_password_raises tests/test_resolve_token.py::test_user_without_password_raises tests/test_resolve_token.py::test_password_without_user_raises tests/test_resolve_token.py::test_secretstr_password_accepted -v`
Expected: all 5 PASS

### Step 1.3: Write and verify resolution chain tests

- [ ] **Add resolution chain tests to `tests/test_resolve_token.py`**

```python
# Append to tests/test_resolve_token.py

# --- Resolution chain tests ---


def test_explicit_token_takes_priority(token_store, stored_entry):
    """--token flag should always win over stored token."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000", token="explicit.flag.token")

    assert result == "explicit.flag.token"


def test_user_password_login_success(token_store, mock_httpx_client):
    """--user/--password should POST to /v1/auth/jwt/login."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status = MagicMock()
    mock_resp.json.return_value = {"access_token": "login.jwt.token"}
    mock_httpx_client.post.return_value = mock_resp

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token(
            "http://localhost:8000", user="admin@example.com", password="secret"
        )

    assert result == "login.jwt.token"
    mock_httpx_client.post.assert_called_once_with(
        "/v1/auth/jwt/login",
        data={"username": "admin@example.com", "password": "secret"},
    )


def test_stored_token_used_when_valid(token_store, stored_entry, mock_httpx_client):
    """Stored token should be used when GET /v1/auth/users/me returns 200."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_httpx_client.get.return_value = mock_response

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "stored.jwt.token"


def test_stored_token_deleted_on_401(token_store, stored_entry, mock_httpx_client):
    """Stored token should be removed on 401 and fall through to guest."""
    from zndraw.auth_utils import resolve_token

    token_store.set("http://localhost:8000", stored_entry)

    mock_401 = MagicMock()
    mock_401.status_code = 401

    mock_guest = MagicMock()
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}

    mock_httpx_client.get.return_value = mock_401
    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "guest.token"
    assert token_store.get("http://localhost:8000") is None


def test_no_stored_token_falls_through_to_guest(token_store, mock_httpx_client):
    """When no stored token exists, should create guest session."""
    from zndraw.auth_utils import resolve_token

    mock_guest = MagicMock()
    mock_guest.raise_for_status = MagicMock()
    mock_guest.json.return_value = {"access_token": "guest.token"}
    mock_httpx_client.post.return_value = mock_guest

    with patch("zndraw.auth_utils.get_token_store", return_value=token_store):
        result = resolve_token("http://localhost:8000")

    assert result == "guest.token"


def test_user_password_error_does_not_leak_password(mock_httpx_client):
    """Failed login error should not contain the password."""
    from zndraw.auth_utils import resolve_token

    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = httpx.HTTPStatusError(
        "401", request=MagicMock(), response=MagicMock()
    )
    mock_httpx_client.post.return_value = mock_resp

    with (
        patch("zndraw.auth_utils.get_token_store"),
        pytest.raises(httpx.HTTPStatusError),
    ):
        resolve_token(
            "http://localhost:8000",
            user="admin@example.com",
            password="super-secret-pass",
        )


def test_migration_warning_when_zndraw_email_set(monkeypatch):
    """Emits DeprecationWarning if ZNDRAW_EMAIL set but ZNDRAW_USER is not."""
    from zndraw.auth_utils import resolve_token

    monkeypatch.setenv("ZNDRAW_EMAIL", "old@example.com")
    monkeypatch.delenv("ZNDRAW_USER", raising=False)

    with pytest.warns(DeprecationWarning, match="ZNDRAW_EMAIL is deprecated"):
        resolve_token("http://localhost:8000", token="t")
```

- [ ] **Run all resolve_token tests**

Run: `uv run pytest tests/test_resolve_token.py -v`
Expected: all PASS

- [ ] **Commit**

```bash
git add src/zndraw/auth_utils.py tests/test_resolve_token.py
git commit -m "feat: add shared resolve_token() in auth_utils with validation and tests"
```

---

## Task 2: Update `connection.py` — remove old `resolve_token()`, add option types

**Files:**
- Modify: `src/zndraw/cli_agent/connection.py`

### Step 2.1: Add `UserOpt` and `PasswordOpt` type aliases

- [ ] **Add new option types after existing `RoomOpt` (line 39)**

Add at `src/zndraw/cli_agent/connection.py:39` after `RoomOpt`:

```python
UserOpt = Annotated[
    str | None,
    typer.Option(
        "--user", envvar="ZNDRAW_USER", help="User email for authentication"
    ),
]
PasswordOpt = Annotated[
    str | None,
    typer.Option(
        "--password", envvar="ZNDRAW_PASSWORD", help="Password for authentication"
    ),
]
```

- [ ] **Run format check**

Run: `uv run ruff format src/zndraw/cli_agent/connection.py`

### Step 2.2: Replace old `resolve_token()` and update `get_connection()`/`get_zndraw()`

- [ ] **Delete the old `resolve_token()` function** (lines 199–274)

Replace with an import-and-wrap pattern that converts `ValueError` from `auth_utils.resolve_token` into CLI-appropriate `die()` calls:

```python
def resolve_token(
    base_url: str,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> str:
    """Resolve auth token — CLI wrapper around shared auth_utils."""
    from zndraw.auth_utils import resolve_token as _resolve_token

    try:
        return _resolve_token(base_url, token=token, user=user, password=password)
    except ValueError as exc:
        die(str(exc), str(exc), 400, EXIT_CLIENT_ERROR)
    except httpx.HTTPStatusError as exc:
        detail = f"Authentication failed: {exc}"
        die("Authentication Failed", detail, 401, EXIT_CLIENT_ERROR)
    except (httpx.RequestError, KeyError) as exc:
        die(
            "Authentication Failed",
            f"Failed to authenticate: {exc}",
            401,
            EXIT_CONNECTION_ERROR,
        )
```

- [ ] **Remove the `import os` at the top of `connection.py`** (line 12) — no longer needed since env var handling moved to `auth_utils.py`. Keep `get_token_store()` in `connection.py` — it is imported by `auth.py` and `admin.py` as a simple factory and is unrelated to token resolution logic.

- [ ] **Update `get_connection()` signature** (line 310):

```python
def get_connection(
    url: str | None,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> Connection:
    """Create a Connection from resolved URL and token."""
    base_url = resolve_url(url)
    resolved_token = resolve_token(base_url, token, user=user, password=password)
    return Connection(base_url=base_url, token=resolved_token)
```

- [ ] **Update `get_zndraw()` signature** (line 325):

```python
def get_zndraw(
    url: str | None,
    token: str | None,
    room: str,
    user: str | None = None,
    password: str | None = None,
) -> ZnDraw:
    """Create a ZnDraw instance from CLI context."""
    from zndraw import ZnDraw

    base_url = resolve_url(url)
    resolved_token = resolve_token(base_url, token, user=user, password=password)
    return ZnDraw(
        url=base_url, room=room, token=resolved_token, create_if_missing=False
    )
```

- [ ] **Run existing tests to verify no breakage**

Run: `uv run pytest tests/test_resolve_token.py tests/test_cli_auth.py -v`
Expected: all PASS

- [ ] **Commit**

```bash
git add src/zndraw/cli_agent/connection.py
git commit -m "refactor: replace connection.py resolve_token with shared auth_utils wrapper"
```

---

## Task 3: Update `auth status` to accept `--user`/`--password`

**Files:**
- Modify: `src/zndraw/cli_agent/auth.py:103-166`

### Step 3.1: Update auth status signature and logic

- [ ] **Add `UserOpt`, `PasswordOpt` imports and update `status()` command**

Update imports at line 14:
```python
from .connection import (
    PasswordOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    cli_error_handler,
    get_token_store,
    resolve_url,
)
```

Replace the `status` function (lines 103–166) with:

```python
@auth_app.command("status")
def status(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
) -> None:
    """Show current authentication identity."""
    with cli_error_handler():
        from zndraw.auth_utils import resolve_token as _resolve_token

        resolved_url = resolve_url(url)
        store = get_token_store()

        # Validate mutual exclusion (raises ValueError -> caught by cli_error_handler)
        if token is not None and (user is not None or password is not None):
            msg = "Cannot combine --token with --user/--password"
            raise ValueError(msg)
        if user is not None and password is None:
            msg = "Missing --password (required when --user is provided)"
            raise ValueError(msg)
        if password is not None and user is None:
            msg = "Missing --user (required when --password is provided)"
            raise ValueError(msg)

        # Determine token and source — no guest fallback
        if token is not None:
            active_token = token
            token_source = "flag"
        elif user is not None and password is not None:
            with httpx.Client(base_url=resolved_url, timeout=10.0) as client:
                resp = client.post(
                    "/v1/auth/jwt/login",
                    data={"username": user, "password": password},
                )
                resp.raise_for_status()
                active_token = resp.json()["access_token"]
            token_source = "login"
        else:
            entry = store.get(resolved_url)
            if entry is not None:
                active_token = entry.access_token
                token_source = "stored"
            else:
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "none",
                    }
                )
                return

        with httpx.Client(
            base_url=resolved_url,
            headers={"Authorization": f"Bearer {active_token}"},
            timeout=10.0,
        ) as client:
            resp = client.get("/v1/auth/users/me")
            if resp.status_code != 200:
                if resp.status_code in (401, 403) and token_source == "stored":
                    store.delete(resolved_url)
                json_print(
                    {
                        "server": resolved_url,
                        "user_id": None,
                        "email": None,
                        "is_superuser": False,
                        "token_source": "expired"
                        if resp.status_code in (401, 403)
                        else "error",
                    }
                )
                return
            user_data = resp.json()

        json_print(
            {
                "server": resolved_url,
                "user_id": user_data.get("id"),
                "email": user_data.get("email"),
                "is_superuser": user_data.get("is_superuser", False),
                "token_source": token_source,
            }
        )
```

Note: `auth status` does NOT use `resolve_token()` — it has its own logic because it must NOT fall through to guest (it reports `"none"` instead). The mutual exclusion validation is duplicated here intentionally.

- [ ] **Run auth tests**

Run: `uv run pytest tests/test_cli_auth.py -v`
Expected: PASS

- [ ] **Commit**

```bash
git add src/zndraw/cli_agent/auth.py
git commit -m "feat: add --user/--password flags to auth status command"
```

---

## Task 4: Update all CLI subcommand signatures (mechanical)

**Files:** All files in `src/zndraw/cli_agent/` listed below.

This is a mechanical change: for every command that takes `token: TokenOpt`, add `user: UserOpt = None, password: PasswordOpt = None` after it, and thread them through to `get_connection()`/`get_zndraw()`/`resolve_token()`.

### Step 4.1: Update subcommands that call `get_connection(url, token)`

- [ ] **Update `rooms.py`** — functions `create_room` (line 58), `set_default_room` (line 137):
  - Add `UserOpt`, `PasswordOpt` to imports
  - Add `user: UserOpt = None, password: PasswordOpt = None` params
  - Change `get_connection(url, token)` → `get_connection(url, token, user, password)`

- [ ] **Update `jobs.py`** — function `status` (line 23):
  - Add `UserOpt`, `PasswordOpt` to imports
  - Add `user: UserOpt = None, password: PasswordOpt = None` param
  - Change `get_connection(url, token)` → `get_connection(url, token, user, password)`

### Step 4.2: Update subcommands that call `get_zndraw(url, token, room)`

For each file below, add `UserOpt`, `PasswordOpt` to imports and update every command function:
- Add `user: UserOpt = None, password: PasswordOpt = None` after `token: TokenOpt`
- Change `get_zndraw(url, token, room)` → `get_zndraw(url, token, room, user, password)`

- [ ] **Update `rooms.py`** — functions `room_info` (78), `lock_room` (92), `unlock_room` (107), `open_room` (122, uses `_token` — add `_user: UserOpt = None, _password: PasswordOpt = None` for `--help` consistency), `list_rooms` (34, calls `resolve_token` directly)
- [ ] **Update `frames.py`** — functions `count` (23), `get` (38), `list_frames` (90), `extend` (139), `export` (180), `delete` (245)
- [ ] **Update `sessions.py`** — functions `list_sessions` (23), `get_camera` (38), `set_camera` (60)
- [ ] **Update `step.py`** — functions `get_step` (23), `set_step` (38)
- [ ] **Update `selection.py`** — functions `get_selection` (26), `set_selection` (42), `clear_selection` (63)
- [ ] **Update `selection_groups.py`** — functions `list_selection_groups` (28), `get_selection_group` (46), `set_selection_group` (68), `delete_selection_group` (96)
- [ ] **Update `bookmarks.py`** — functions `list_bookmarks` (28), `set_bookmark` (46), `delete_bookmark` (72)
- [ ] **Update `figures.py`** — functions `list_figures` (30), `get` (48), `set_figure` (69), `delete` (95)
- [ ] **Update `screenshots.py`** — functions `list_screenshots` (29), `request_screenshot` (48), `get` (83)
- [ ] **Update `extensions.py`** — functions `list_extensions` (27), `describe` (42), `run_extension` (69)
- [ ] **Update `chat.py`** — functions `list_messages` (23), `send_message` (44)
- [ ] **Update `geometries.py`** — functions `list_geometries` (29), `get` (87), `set_geometry` (108), `toggle_geometry` (153), `set_prop` (182), `delete` (230)
- [ ] **Update `gif.py`** — function `capture` (146)
- [ ] **Update `presets.py`** — functions `list_presets` (30), `get_preset` (49), `load_preset` (66), `apply_preset` (86), `save_preset` (104), `reset_preset` (159), `export_preset` (174), `delete_preset` (197)

### Step 4.3: Update subcommands that call `resolve_token()` directly

- [ ] **Update `mount.py`** — function `mount_cmd` (line 18):
  - Add `UserOpt`, `PasswordOpt` to imports
  - Add params, thread through to `resolve_token(resolved_url, token, user, password)`

- [ ] **Update `admin.py`** — functions `list_users` (line 29), `login_as_user` (line 53):
  - Add `UserOpt`, `PasswordOpt` to imports
  - Add params, thread through to `resolve_token(resolved_url, token, user, password)`

- [ ] **Update `jobs.py`** — function `list_jobs` (line 36) also calls `get_zndraw`:
  - Thread `user`, `password` through

### Step 4.4: Verify and commit

- [ ] **Run format and import sort**

Run: `uv run ruff format src/zndraw/cli_agent/ && uv run ruff check --select I --fix src/zndraw/cli_agent/`

- [ ] **Run type check on CLI module**

Run: `uv run pyright src/zndraw/cli_agent/`
Expected: no new errors

- [ ] **Run full test suite**

Run: `uv run pytest tests/ -x -q`
Expected: all PASS

- [ ] **Commit**

```bash
git add src/zndraw/cli_agent/
git commit -m "feat: add --user/--password flags to all CLI subcommands"
```

---

## Task 5: Update `ZnDraw` client to use shared `resolve_token()`

**Files:**
- Modify: `src/zndraw/client/core.py:130-212` (`__post_init__`, `_try_stored_token`)
- Modify: `src/zndraw/client/core.py:541-598` (`list_rooms`, `login`)

### Step 5.1: Replace `__post_init__` auth chain

- [ ] **Replace lines 130–212** of `src/zndraw/client/core.py`

Replace the `__post_init__` method body (keep everything up to `# Authenticate` the same, then replace the auth block):

```python
    def __post_init__(self) -> None:
        """Initialize the client (REST-only, socket connects lazily)."""
        import atexit

        from zndraw.auth_utils import resolve_token

        # Normalize password to SecretStr
        if isinstance(self.password, str):
            self.password = SecretStr(self.password)

        # Generate room ID if not provided
        if self.room is None:
            self.room = str(uuid.uuid4())

        # Resolve URL (auto-discover if None)
        self.url = self._resolve_url(self.url)

        # Create API manager (token may be None initially)
        self.api = APIManager(url=self.url, room_id=self.room, token=self.token)

        # Authenticate via shared resolve_token
        self.api.token = resolve_token(
            self.url,
            token=self.token,
            user=self.user,
            password=self.password,
        )

        # Populate self.user for guest/stored-token sessions
        if self.user is None:
            resp = self.api.http.get(
                "/v1/auth/users/me",
                headers={"Authorization": f"Bearer {self.api.token}"},
            )
            if resp.status_code == 200:
                self.user = resp.json().get("email")

        # Create socket manager (no connection yet -- connects lazily)
        self.socket = SocketManager(zndraw=self)

        # Create job manager (zero-cost until first register())
        self._jobs = JobManager(
            api=self.api,
            tsio=self.socket.tsio,
            execute=self._execute_task if self.auto_pickup else None,
            heartbeat_interval=self.heartbeat_interval,
            polling_interval=self.polling_interval,
        )

        # Verify/create room via REST and seed frame count cache
        try:
            info = self.api.get_room_info()
            self.cached_length = info.get("frame_count", 0)
        except KeyError:
            if not self.create_if_missing:
                raise
            self.api.create_room(copy_from=self.copy_from)
            self.cached_length = 0

        # Ensure cleanup on interpreter exit
        atexit.register(self.disconnect)
```

- [ ] **Delete `_try_stored_token()` method** (lines 189–212)

**Breaking change:** The old code fell back to `Settings().guest_password` when `user` was provided without `password`. The new code raises `ValueError` instead. Callers must now pass `password` explicitly. This is intentional — clean implementation, no backwards compat shims.

- [ ] **Update `list_rooms()` classmethod** (lines 541–567):

```python
    @classmethod
    def list_rooms(
        cls,
        url: str | None = None,
        *,
        token: str | None = None,
        search: str | None = None,
    ) -> list[dict[str, Any]]:
        """List all rooms on the server."""
        from zndraw.auth_utils import resolve_token

        resolved = cls._resolve_url(url)
        resolved_token = resolve_token(resolved, token=token)
        api = APIManager(url=resolved, room_id="", token=resolved_token)
        try:
            return api.list_rooms(search=search)
        finally:
            api.close()
```

- [ ] **Update `login()` classmethod** (lines 569–598):

```python
    @classmethod
    def login(
        cls,
        url: str | None = None,
        username: str = "",
        password: str = "",
    ) -> str:
        """Authenticate and return a JWT token."""
        from zndraw.auth_utils import resolve_token

        resolved = cls._resolve_url(url)
        return resolve_token(resolved, user=username, password=password)
```

### Step 5.2: Verify and commit

- [ ] **Run format**

Run: `uv run ruff format src/zndraw/client/core.py && uv run ruff check --select I --fix src/zndraw/client/core.py`

- [ ] **Run tests**

Run: `uv run pytest tests/ -x -q`
Expected: all PASS

- [ ] **Commit**

```bash
git add src/zndraw/client/core.py
git commit -m "refactor: use shared resolve_token in ZnDraw client, remove _try_stored_token"
```

---

## Task 6: Rename environment variables in Docker and docs

**Files:**
- Modify: `docker/standalone/docker-compose.yaml:54`
- Modify: `docker/production/docker-compose.yaml:106`
- Modify: `skills/zndraw/SKILL.md:76`

### Step 6.1: Update Docker Compose files

- [ ] **`docker/standalone/docker-compose.yaml` line 54:**
  Change `ZNDRAW_EMAIL` → `ZNDRAW_USER`

- [ ] **`docker/production/docker-compose.yaml` line 106:**
  Change `ZNDRAW_EMAIL` → `ZNDRAW_USER`

### Step 6.2: Update skill docs

- [ ] **`skills/zndraw/SKILL.md` line 76:**
  Change `ZNDRAW_EMAIL` → `ZNDRAW_USER` and update the description to match the new resolution chain (3 tiers, not 4)

### Step 6.3: Verify no remaining references

- [ ] **Search for any remaining `ZNDRAW_EMAIL` references**

Run: `uv run ruff format . && grep -r "ZNDRAW_EMAIL" src/ tests/ docker/ skills/ --include="*.py" --include="*.yaml" --include="*.yml" --include="*.md" | grep -v auth_utils.py`
Expected: no results (only `auth_utils.py` migration warning should reference it)

- [ ] **Commit**

```bash
git add docker/ skills/
git commit -m "chore: rename ZNDRAW_EMAIL to ZNDRAW_USER in Docker and docs"
```

---

## Task 7: Final verification

### Step 7.1: Full test suite

- [ ] **Run complete test suite**

Run: `uv run pytest tests/ -q`
Expected: all tests PASS

### Step 7.2: Type check

- [ ] **Run pyright**

Run: `uv run pyright .`
Expected: no new errors

### Step 7.3: Verify CLI help shows new flags

- [ ] **Spot-check a command's help output**

Run: `uv run zndraw-cli rooms create --help`
Expected: should show `--user` and `--password` options alongside `--token`
