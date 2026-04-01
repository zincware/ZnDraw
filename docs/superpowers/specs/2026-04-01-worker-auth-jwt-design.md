# Worker Auth JWT Redesign

**Date:** 2026-04-01
**Branch:** `fix/worker-auth-jwt`
**Problem:** The internal worker user (`worker@internal.user`) is created with a
well-known default password (`zndraw-worker`). Anyone who reads the source code
can log in as a superuser on any deployment that hasn't changed the default.

## Decision

Replace the shared-password model with per-task JWT tokens minted by the server
at dispatch time. The worker password becomes an auto-generated UUID that never
needs to be configured or shared. External TaskIQ workers receive a short-lived
JWT in each job payload instead of holding a static credential.

## Design

### 1. Worker user creation

On startup (or via `db-init` in production), `ensure_internal_worker()` creates
the `worker@internal.user` superuser with a **random UUID password**. The
password is write-only — it satisfies the DB schema but is never read back or
used for authentication.

```python
# database.py — ensure_internal_worker
password = SecretStr(str(uuid.uuid4()))
```

Multi-replica safety: only the process with `init_db_on_startup=true` (the
`db-init` service in Docker, or the single server in standalone) calls
`ensure_internal_worker()`. Replicas never touch the worker user row. JWT
minting at dispatch time only needs the user's DB record + the shared JWT
signing secret — no password involved.

### 2. Block worker login and registration

Add explicit guards that reject the internal worker email on public auth
endpoints:

- **Login** (`POST /v1/auth/jwt/login`): reject with 403 if username matches
  `settings.internal_worker_email`.
- **Register** (`POST /v1/auth/register`): reject with 403 if email matches.
  Registration would already fail with "user exists", but an explicit check
  avoids leaking that the email is taken and makes the intent clear.

The internal worker email moves from a module-level constant in `database.py` to
a config field:

```python
# config.py — Settings
internal_worker_email: str = "worker@internal.user"
```

### 3. JWT minting at dispatch time (with DI)

A new FastAPI dependency mints a fresh JWT for the worker user on each task
submission:

```python
async def get_worker_token(
    session: SessionDep,
    strategy: JWTStrategyDep,
    settings: SettingsDep,
) -> str:
    user = await lookup_worker_user(session, settings.internal_worker_email)
    return await strategy.write_token(user)

WorkerTokenDep = Annotated[str, Depends(get_worker_token)]
```

The `submit_task` route in `router.py` declares `worker_token: WorkerTokenDep`
and passes it through to `kiq()`:

```python
await internal_registry.tasks[job.full_name].kiq(
    task_id=str(task.id),
    room_id=room_id,
    payload=request.payload,
    token=worker_token,
)
```

Token properties:
- **TTL:** matches user token lifetime (default 1 hour via
  `ZNDRAW_AUTH_TOKEN_LIFETIME_SECONDS`)
- **Scope:** standard fastapi-users JWT with `sub` = worker user ID
- **Signing:** HS256 with `ZNDRAW_AUTH_SECRET_KEY` (same as all user tokens)

### 4. Executor changes

`InternalExtensionExecutor` is simplified:

**Removed fields:** `worker_email`, `worker_password`
**Kept fields:** `base_url`
**New `__call__` parameter:** `token: str`

```python
@dataclass
class InternalExtensionExecutor:
    base_url: str

    async def __call__(
        self,
        extension_cls: type,
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
        token: str,
    ) -> None:
        def _run() -> None:
            vis = ZnDraw(url=base_url, room=room_id, token=token)
            # ... run extension, report status ...
```

The task function in `registry.py` adds `token: str` to its signature and passes
it through to the executor.

`broker.py` simplifies — the executor only needs `base_url` at creation time.
No auth settings are read at module level.

### 5. Config cleanup — full sweep

**Remove from code:**
- `Settings.worker_password` field in `config.py`
- `WORKER_EMAIL` constant in `database.py` (replaced by
  `Settings.internal_worker_email`)
- `worker_password` parameter on `ensure_internal_worker()`
- `worker_email` / `worker_password` fields on `InternalExtensionExecutor`
- Any references to `worker_password` in `broker.py`, `database.py`, tests

**Remove from Docker:**
- `ZNDRAW_SERVER_WORKER_PASSWORD=zndraw-worker` from `docker/templates/.env`
- Any references in `docker/standalone/README.md`,
  `docker/production/README.md`, `docker/standalone/docker-compose.yaml`,
  `docker/production/docker-compose.yaml`

**Remove from docs:**
- Any references in `docs/superpowers/specs/` or `docs/superpowers/plans/`

**Grep verification:** after all changes, `grep -r worker_password` and
`grep -r WORKER_PASSWORD` across the entire repo must return zero hits (excluding
git history and this spec).

**Add:**
- `Settings.internal_worker_email: str = "worker@internal.user"` in `config.py`

### 6. Test updates

- Existing tests that set `ZNDRAW_SERVER_WORKER_PASSWORD` env vars: remove those
  env vars. The worker user is created with a random password automatically.
- Add a test that `POST /v1/auth/jwt/login` with the worker email returns 403.
- Add a test that `POST /v1/auth/register` with the worker email returns 403.
- Update `InternalExtensionExecutor` tests to pass `token` instead of
  `worker_email` / `worker_password`.

## Security properties

| Property | Before | After |
|----------|--------|-------|
| Worker password | Static default `"zndraw-worker"`, configurable | Random UUID per startup, never exposed |
| Public login | Worker user loginable via `/auth/jwt/login` | Blocked with 403 |
| Registration | Worker email registrable (fails with "exists") | Blocked with 403 (no email leak) |
| Redis exposure | Password never in Redis | JWT in Redis per task (1-hour TTL) |
| TaskIQ worker config | Needs `ZNDRAW_SERVER_WORKER_PASSWORD` env var | Needs no auth config |
| Multi-replica | All replicas need same password | Stateless — any replica mints JWTs |

## Files changed

| File | Change |
|------|--------|
| `src/zndraw/config.py` | Remove `worker_password`, add `internal_worker_email` |
| `src/zndraw/database.py` | Auto-gen UUID password, use `settings.internal_worker_email`, export lookup helper |
| `src/zndraw/executor.py` | Drop email/password fields, accept `token` in `__call__` |
| `src/zndraw/broker.py` | Simplify — executor only needs `base_url` |
| `src/zndraw_joblib/registry.py` | Add `token: str` to task function signature, pass to executor |
| `src/zndraw_joblib/router.py` | Add `WorkerTokenDep`, pass token to `kiq()` |
| `src/zndraw/routes/auth.py` | Add login/register guards for internal worker email |
| `docker/templates/.env` | Remove `ZNDRAW_SERVER_WORKER_PASSWORD` |
| `docker/*/README.md` | Remove worker password from config tables |
| `tests/` | Update env vars, add login/register block tests |
