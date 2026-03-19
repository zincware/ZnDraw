# CLI Auth Consistency

## Problem

`--token` is a visible Typer flag on every CLI command, but `ZNDRAW_EMAIL`/`ZNDRAW_PASSWORD` are consumed silently inside `resolve_token()` from environment variables — invisible in `--help` and undiscoverable without reading source code. Additionally, token resolution logic is duplicated between `cli_agent/connection.py` and `ZnDraw.__post_init__`.

## Solution

1. Make all auth methods explicit CLI flags with envvar fallbacks.
2. Extract token resolution into a single shared utility.

## Shared Utility: `src/zndraw/auth_utils.py`

New module with a single public function:

```python
def resolve_token(
    base_url: str,
    token: str | None = None,
    user: str | None = None,
    password: SecretStr | str | None = None,
) -> str:
```

Takes `base_url` and creates short-lived `httpx.Client` instances internally for network calls. Accepts `SecretStr | str | None` for `password` so the `ZnDraw` client can pass its field directly without unwrapping.

### Validation (fail fast, before any network call)

- `token` + (`user` or `password`) -> error: "Cannot combine --token with --user/--password"
- `user` without `password` -> error: "Missing --password (required when --user is provided)"
- `password` without `user` -> error: "Missing --user (required when --password is provided)"

### Resolution Chain (3 tiers)

Explicit credentials always take priority over implicit state (stored tokens). This is intentional: when a user passes `--user`/`--password` or `--token`, that should win regardless of what is stored locally.

1. **Explicit credentials** — `token` provided: return it. OR `user`+`password`: login via `POST /v1/auth/jwt/login`, return access token. Token is NOT stored.
2. **Stored token** — look up `~/.zndraw/tokens.json` via `TokenStore`, validate against `GET /v1/auth/users/me`, return on 200 or delete on 401 and fall through.
3. **Guest fallback** — `POST /v1/auth/guest`, return auto-generated guest token.

### Error Handling

Raises exceptions (e.g. `ValueError`, `httpx.HTTPStatusError`). Each consumer wraps errors in its own way:
- CLI uses `die()` for RFC 9457 problem JSON to stderr.
- `ZnDraw` client raises its own exceptions.

## CLI Option Types

New shared Typer annotated types in `connection.py` alongside existing `UrlOpt`/`TokenOpt`/`RoomOpt`:

```python
UserOpt = Annotated[
    str | None,
    typer.Option("--user", envvar="ZNDRAW_USER", help="User email for authentication"),
]
PasswordOpt = Annotated[
    str | None,
    typer.Option("--password", envvar="ZNDRAW_PASSWORD", help="Password for authentication"),
]
```

## Command Signatures

Every command that currently takes `url`/`token` also gets `user`/`password`. This adds flags to ~60 command functions across 18 files — a mechanical change, but intentional for discoverability: every command shows `--user`/`--password` in `--help`.

```python
def some_command(
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    ...
):
    conn = get_connection(url, token, user, password)
```

### Auth Subcommands

- `auth login` — `--url` only. Browser device-code flow. Stores token to `~/.zndraw/tokens.json`.
- `auth status` — `--url`, `--token`, `--user`/`--password`. Validates credential against server. No guest fallback. Same mutual exclusion rules apply. When no credentials are provided and no stored token exists, reports `{"token_source": "none"}` (current behavior preserved).
- `auth logout` — `--url` only. Clears stored token.

### `get_connection()` / `get_zndraw()` Updated Signatures

```python
def get_connection(
    url: str | None,
    token: str | None,
    user: str | None = None,
    password: str | None = None,
) -> Connection:
```

```python
def get_zndraw(
    url: str | None,
    token: str | None,
    room: str,
    user: str | None = None,
    password: str | None = None,
) -> ZnDraw:
```

## ZnDraw Client

`ZnDraw.__post_init__` replaces its bespoke auth chain and `_try_stored_token()` with a call to shared `resolve_token()`. The `_try_stored_token()` method is removed entirely.

The `ZnDraw.list_rooms()` and `ZnDraw.login()` class methods also use the shared `resolve_token()` instead of their ad-hoc auth patterns.

## Environment Variables

Hard rename, no backwards compatibility:

- `ZNDRAW_EMAIL` -> `ZNDRAW_USER`
- `ZNDRAW_PASSWORD` unchanged

As a migration aid, `resolve_token()` emits a warning to stderr if `ZNDRAW_EMAIL` is set but `ZNDRAW_USER` is not, telling users to rename the variable.

Docker Compose files updated accordingly.

## Files Touched

| File | Change |
|---|---|
| `src/zndraw/auth_utils.py` | **New** — shared `resolve_token()` |
| `src/zndraw/cli_agent/connection.py` | Remove old `resolve_token()`, add `UserOpt`/`PasswordOpt`, update `get_connection()`/`get_zndraw()` |
| `src/zndraw/cli_agent/auth.py` | Add `--user`/`--password` to `auth status` |
| `src/zndraw/cli_agent/*.py` (all subcommands) | Add `user: UserOpt`, `password: PasswordOpt` to signatures |
| `src/zndraw/client/core.py` | Replace bespoke auth chain + `_try_stored_token()` with shared `resolve_token()`, including `list_rooms()` and `login()` |
| `docker/*/docker-compose.yaml` | `ZNDRAW_EMAIL` -> `ZNDRAW_USER` |
| `skills/zndraw/SKILL.md` | Update `ZNDRAW_EMAIL` references |
| `tests/test_resolve_token.py` | Rewrite for new `resolve_token(base_url, token, user, password)` signature and validation tests |
| `tests/test_cli_auth.py` | Add `--user`/`--password` coverage for `auth status` |
