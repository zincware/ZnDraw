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
    client: httpx.Client,
    token: str | None = None,
    user: str | None = None,
    password: str | None = None,
) -> str:
```

Takes an `httpx.Client` (callers already have one) — no unions, each consumer passes its own client.

### Validation (fail fast, before any network call)

- `token` + (`user` or `password`) -> error: "Cannot combine --token with --user/--password"
- `user` without `password` -> error: "Missing --password (required when --user is provided)"
- `password` without `user` -> error: "Missing --user (required when --password is provided)"

### Resolution Chain (3 tiers)

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

Every command that currently takes `url`/`token` also gets `user`/`password`:

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
- `auth status` — `--url`, `--token`, `--user`/`--password`. Validates credential against server. No guest fallback. Same mutual exclusion rules apply.
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

`ZnDraw.__post_init__` replaces its bespoke auth chain and `_try_stored_token()` with a call to shared `resolve_token()` via `self.api.http`. The `_try_stored_token()` method is removed entirely.

## Environment Variables

Hard rename, no backwards compatibility:

- `ZNDRAW_EMAIL` -> `ZNDRAW_USER`
- `ZNDRAW_PASSWORD` unchanged

Docker Compose files updated accordingly.

## Files Touched

| File | Change |
|---|---|
| `src/zndraw/auth_utils.py` | **New** — shared `resolve_token()` |
| `src/zndraw/cli_agent/connection.py` | Remove old `resolve_token()`, add `UserOpt`/`PasswordOpt`, update `get_connection()`/`get_zndraw()` |
| `src/zndraw/cli_agent/auth.py` | Add `--user`/`--password` to `auth status` |
| `src/zndraw/cli_agent/*.py` (all subcommands) | Add `user: UserOpt`, `password: PasswordOpt` to signatures |
| `src/zndraw/client/core.py` | Replace bespoke auth chain + `_try_stored_token()` with shared `resolve_token()` |
| `docker-compose*.yml` | `ZNDRAW_EMAIL` -> `ZNDRAW_USER` |
| `tests/test_resolve_token.py` | Update for new API and validation tests |
| `tests/test_cli_auth.py` | Add `--user`/`--password` coverage for `auth status` |
