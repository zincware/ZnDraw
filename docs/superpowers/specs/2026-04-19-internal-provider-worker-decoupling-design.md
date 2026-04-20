# `@internal` Provider / Worker Decoupling — Design Spec

**Date:** 2026-04-19
**Branch:** `spec/dockview-ui-redesign`
**PR:** [zincware/ZnDraw#920](https://github.com/zincware/ZnDraw/pull/920)
**Supersedes (corrects):** the sweeper workaround in commit `f34bad0c` (Phase 5). That commit patched a symptom; this spec fixes the root cause.

---

## 1. Problem

`@internal` providers are seeded at server startup and dispatched by the in-process taskiq worker. Taskiq has no heartbeat — brokers are just queues that in-process consumers poll. `@internal` is "the server itself," not a remote participant.

`ProviderRecord.worker_id: Mapped[UUID]` is a **non-nullable FK** to `worker.id`. When Phase 2 seeded the default `@internal:filesystem:FilesystemRead` row, `ensure_internal_worker_row` manufactured a fake `Worker` row just to satisfy the FK. That ghost Worker never heartbeats. After `worker_timeout_seconds=60`, `cleanup_stale_workers` finds it, `cleanup_worker` cascade-deletes every provider owned by it — including the `@internal` row the server itself registered. End-to-end symptom: the filesystem browser works for the first minute of a server session and then 404s forever.

Commit `f34bad0c` patched this by excluding workers that own `@internal` providers from the stale-worker scan. It works, but the root cause is that `@internal` providers shouldn't have a `worker_id` at all — mirror how `@internal` jobs have no `WorkerJobLink`.

### Architectural asymmetry (current state, post Phase 5)

| | Remote | `@internal` |
|---|---|---|
| Job ownership | `WorkerJobLink` rows | None — dispatched from in-memory taskiq registry |
| Provider ownership | `ProviderRecord.worker_id → worker.id` | Ghost `Worker` row whose only purpose is to satisfy the FK |

Jobs are coherent. Providers are not. This spec removes the asymmetry.

---

## 2. Design

### 2.1 Schema

`ProviderRecord.worker_id` becomes nullable:

```python
# src/zndraw_joblib/models.py
class ProviderRecord(Base):
    ...
    worker_id: Mapped[UUID | None] = mapped_column(
        ForeignKey("worker.id", ondelete="CASCADE"), index=True, nullable=True
    )
```

The `ondelete="CASCADE"` semantics are preserved for remote providers: when a real Worker row is deleted, its providers are cascade-deleted at the DB level. For `@internal` rows, `worker_id IS NULL` never cascades.

No Alembic migration file. Tables are created via `SQLModel.metadata.create_all` at startup — recreating the DB on each test run / fresh deploy. Acceptable because `@internal` providers are a brand-new concept on this branch and have never shipped in a stable release (per the "no compat weight on unshipped" rule in the project's feedback memory).

### 2.2 Bootstrap (`src/zndraw/database.py`)

**Remove** `ensure_internal_worker_row` — the function and all call sites. The internal-worker User row is still created by `ensure_internal_worker` (it's needed for `get_worker_token` JWT minting via the cached `app.state.internal_worker_user`).

`init_database` simplifies to:

```python
if settings.filebrowser_path is not None:
    from zndraw_joblib.registry import ensure_internal_providers

    async with session_maker() as session:
        result = await session.exec(
            select(User).where(User.email == settings.internal_worker_email)
        )
        internal_user = result.one()
        internal_user_id = internal_user.id

    await ensure_internal_providers(
        _collect_providers(),
        session_maker,
        user_id=internal_user_id,
    )
else:
    # stale-row cleanup — unchanged
    ...
```

### 2.3 Registry (`src/zndraw_joblib/registry.py`)

`ensure_internal_providers` drops the `worker_id` kwarg:

```python
async def ensure_internal_providers(
    providers: list[type[Provider]],
    session_factory: ...,
    *,
    user_id: UUID,
) -> None:
```

Seeded rows set `worker_id=None`:

```python
ProviderRecord(
    room_id="@internal",
    category=category,
    name=name,
    schema_=schema,
    content_type=content_type,
    user_id=user_id,
    worker_id=None,
)
```

The IntegrityError requery-update path (from Phase 2 Task 2.4) is preserved verbatim — just also clears `worker_id` on the existing-row branch (in case a stale non-null value is present).

### 2.4 Sweeper (`src/zndraw_joblib/sweeper.py`)

**Revert commit `f34bad0c`'s edit to `cleanup_stale_workers`** — restore the simple query:

```python
result = await session.exec(select(Worker).where(Worker.last_heartbeat < cutoff))
```

`cleanup_worker`'s provider-delete query is already correctly scoped:

```python
select(ProviderRecord).where(ProviderRecord.worker_id == worker.id)
```

`NULL == <uuid>` is false in SQL — `@internal` rows with `worker_id IS NULL` are naturally excluded. No special case needed.

**Invariant preserved**: remote providers (those with a real `worker_id`) continue to be cleaned up when their owning Worker goes stale. Heartbeat-based sweep remains correct for the remote case — the only path that changes semantics is `@internal`, which now has no Worker to sweep in the first place.

### 2.5 Callsite audit

Anywhere that reads `provider.worker_id` must handle `None`:

- `src/zndraw_joblib/router.py::upload_provider_result` — authorisation check uses `provider.user_id`, not `worker_id`. Nothing to change.
- `src/zndraw_joblib/router.py::delete_provider` — already refuses `@internal`. Remote providers still have valid `worker_id`. Nothing to change.
- **`src/zndraw_joblib/schemas.py::ProviderResponse.worker_id: UUID` — load-bearing change.** Must become `worker_id: UUID | None = None`. Without this, `ProviderResponse.from_record(record)` raises `ValidationError` when passed an `@internal` row with `worker_id=None`. Every GET / PUT that returns a `ProviderResponse` (register, list, read) hits this.
- `cleanup_worker`: the delete is filtered by concrete `worker_id == worker.id`; NULL rows unaffected.

An end-of-task grep (`grep -rn "provider.worker_id\|\.worker_id" src/zndraw_joblib src/zndraw`) runs during implementation to surface anything missed.

### 2.6 Tests

- **Rewrite** `tests/zndraw/test_internal_worker_sweeper.py`. New shape: create an `@internal` provider with `worker_id=None` AND an unrelated stale `Worker` row. Assert the sweeper cleans up the stale Worker normally but the `@internal` provider survives. This catches both directions of the invariant (remote sweeping still works, `@internal` is not touched).
- **Update** all three `ensure_internal_providers(..., worker_id=...)` callsites in `tests/zndraw_joblib/test_registry.py` (~ lines 267–282, 304–310, 358–365) — drop the `worker_id=` kwarg.
- **Add** `tests/zndraw_joblib/test_registry.py::test_ensure_internal_providers_stores_null_worker_id` — assert seeded row has `worker_id is None`.
- **Add** `tests/zndraw_joblib/test_providers.py` (or adjacent) coverage: serializing an `@internal` ProviderRecord via `ProviderResponse.from_record` succeeds with `worker_id=None` in the payload. Catches the Pydantic schema regression directly.
- **Hygiene (optional)**: existing `@internal` seeds in `tests/zndraw_joblib/test_providers.py` (~ lines 171, 218, 452, 672) construct `ProviderRecord(..., worker_id=<real worker>)`. These continue to pass but drift from how the bootstrap now seeds (`worker_id=None`). Flip them to `None` for consistency if the diff stays small; otherwise leave a follow-up ticket.

---

## 3. What stays as-is

- `ensure_internal_worker` (the User row). Still required for JWT minting.
- Cached `app.state.internal_worker_user`. Still needed; decoupling doesn't touch auth.
- Phase 1–4 work. Untouched.
- Phase 5 LoadFile registry addition. Unrelated to this fix; keeps value.
- Remote provider flow. Unchanged — `worker_id` stays set, heartbeat still enforces ownership, sweeper still cleans up stale remote workers and their providers via the existing query.

---

## 4. Commit plan

Single commit on `spec/dockview-ui-redesign`:

```text
fix(providers): decouple @internal providers from the Worker table

- ProviderRecord.worker_id: nullable; @internal rows store NULL.
- ensure_internal_worker_row deleted; @internal providers have no Worker
  row — they are server-owned, mirroring how @internal jobs have no
  WorkerJobLink.
- Revert the sweeper special-case from f34bad0c. cleanup_worker's
  WHERE worker_id == :id delete naturally excludes NULL rows.
- Update tests: rewrite test_internal_worker_sweeper.py to cover both
  invariants (remote still swept; @internal untouched).
```

---

## 5. Acceptance criteria

1. `uv run pytest tests/zndraw tests/zndraw_joblib -q` — all green (same count as Phase 5, +/- the new/rewritten tests).
2. `grep -rn "ensure_internal_worker_row\|worker_id=internal_worker_id" src/` — zero matches.
3. The rewritten `test_internal_worker_sweeper.py` unit test passes: it seeds a stale remote Worker AND an `@internal` provider with `worker_id=None`, calls `cleanup_stale_workers(session, timedelta(seconds=0))`, and asserts the stale remote worker IS swept while the `@internal` provider IS preserved. No 90-second end-to-end wait required.
4. Remote provider path unchanged: a test client registering a regular filesystem provider still has it cleaned up when its Worker goes stale.
5. `ProviderResponse.from_record` round-trips an `@internal` row (`worker_id=None`) without raising `ValidationError`.

---

## 6. Out of scope

- Making `ProviderRecord.user_id` nullable. `@internal` providers still have a user (the internal worker user, for JWT POST-back auth). Remote providers still have a user (the owning human).
- Broader refactor of Worker heartbeat semantics. The remote ownership model is sound; only `@internal`'s misuse of it was wrong.
- Frontend "FS goes offline but sidebar still visible" follow-up. Made moot by this fix — once `@internal` providers no longer disappear on their own, the stale icon scenario can't reproduce from the sweeper path. Any other scenario is a separate ticket.

---

## References

- PR: https://github.com/zincware/ZnDraw/pull/920
- Phase 5 workaround commit: `f34bad0c fix(providers): survive sweeper; register LoadFile as @internal modifier`
- Phase 2 spec (original @internal design): `docs/superpowers/specs/2026-04-17-internal-providers-design.md`
- Sweeper source: `src/zndraw_joblib/sweeper.py:140-152` (provider delete query), `:190` (stale worker query)
