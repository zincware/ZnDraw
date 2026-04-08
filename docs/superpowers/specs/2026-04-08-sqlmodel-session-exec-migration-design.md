# Migrate `src/` to SQLModel `session.exec()`

**Date:** 2026-04-08
**Related issue:** [zincware/ZnDraw#890](https://github.com/zincware/ZnDraw/issues/890)
**Status:** Approved, ready for implementation

## Problem

SQLModel's `AsyncSession.execute()` is decorated with `@typing_extensions.deprecated(...)` and emits a `DeprecationWarning` on every call. Across `src/`, there are **77 call sites in 17 files**, producing 15,000+ warnings per test run and drowning out real warnings in CI.

The GitHub issue names `zndraw_joblib/registry.py:108` as the exemplar, but patching only that line leaves the vast majority of warnings untouched. The joblib, zndraw, and zndraw_auth packages are now all part of `src/` in the same repo; any migration must cover all three or a CI warning filter would fail on unmigrated files.

Beyond noise reduction, the SQLModel upstream recommendation (verified via context7 against sqlmodel.tiangolo.com) is unambiguous: **"It is recommended to always use `session.exec()`."** It auto-handles `.scalars()` for SELECT statements, provides proper editor autocompletion, typed results, and inline error checking. This is both a deprecation fix *and* an SOTA alignment.

## Scope

Full inventory of `session.execute()` call sites under `src/` (as of 2026-04-08):

| File | Call sites |
|---|---|
| `src/zndraw_joblib/router.py` | 40 |
| `src/zndraw_joblib/sweeper.py` | 9 |
| `src/zndraw_joblib/registry.py` | 1 |
| `src/zndraw_joblib/dependencies.py` | 1 |
| `src/zndraw/routes/rooms.py` | 5 |
| `src/zndraw/routes/presets.py` | 4 |
| `src/zndraw/routes/chat.py` | 3 |
| `src/zndraw/routes/admin.py` | 2 |
| `src/zndraw/routes/screenshots.py` | 2 |
| `src/zndraw/routes/geometries.py` | 1 |
| `src/zndraw/routes/figures.py` | 1 |
| `src/zndraw/routes/selection_groups.py` | 1 |
| `src/zndraw/routes/bookmarks.py` | 1 |
| `src/zndraw/routes/frames.py` | 1 |
| `src/zndraw/database.py` | 1 |
| `src/zndraw_auth/db.py` | 1 |
| `src/zndraw_auth/cli_login.py` | 3 |
| **Total** | **77 across 17 files** |

**Out of scope:** `tests/**`, everything outside `src/`. `src/zndraw/socketio.py` already uses `session.exec()` and serves as a reference pattern.

## Source-verified foundations

Reading `sqlmodel/ext/asyncio/session.py` directly, `AsyncSession.exec()` has three overloads:

```python
@overload
async def exec(self, statement: Select[_T], ...) -> TupleResult[_T]: ...
@overload
async def exec(self, statement: SelectOfScalar[_T], ...) -> ScalarResult[_T]: ...
@overload
async def exec(self, statement: UpdateBase, ...) -> CursorResult[Any]: ...
```

Consequences:

1. SELECT-of-model returns a `ScalarResult` with `.one()`, `.one_or_none()`, `.first()`, `.all()`, and direct iteration — no manual `.scalars()`.
2. **Bulk `update()` / `delete()` are first-class** via the `UpdateBase` overload and return `CursorResult[Any]` with `.rowcount`. No type ignore, no mixing with raw `session.execute()`, no fallback needed.
3. The deprecation warning is raised by `@typing_extensions.deprecated` on the `execute()` method (line 109 of `sqlmodel/ext/asyncio/session.py`). It does **not** fire on `exec()`.
4. `typing_extensions.deprecated` uses `stacklevel=2`, which makes the warning's `module` attribute resolve to the **caller's** module (e.g. `zndraw_joblib.registry`), not to `sqlmodel.*`. This has implications for the CI warning filter — see below.

## Mechanical rewrite rules

All patterns below are derived from actual call sites in the codebase.

### Rule 1 — Imports

Swap `select` to the SQLModel re-export. Keep everything else from sqlalchemy — SQLModel does not re-export `func`, `and_`, `or_`, `update`, `delete`, etc., and they interop transparently.

```python
# Before
from sqlalchemy import func, select, update

# After
from sqlalchemy import func, update
from sqlmodel import select
```

### Rule 2 — SELECT returning model instances

| Before | After |
|---|---|
| `result = await session.execute(select(Model).where(...))`<br>`obj = result.scalar_one_or_none()` | `result = await session.exec(select(Model).where(...))`<br>`obj = result.one_or_none()` |
| `result = await session.execute(select(Model).where(...))`<br>`obj = result.scalar_one()` | `result = await session.exec(select(Model).where(...))`<br>`obj = result.one()` |
| `result = await session.execute(select(Model).where(...))`<br>`rows = result.scalars().all()` | `result = await session.exec(select(Model).where(...))`<br>`rows = result.all()` |
| `rows = (await session.execute(stmt)).scalars().all()` | `rows = (await session.exec(stmt)).all()` |
| `total = (await session.execute(count_stmt)).scalar_one()` | `total = (await session.exec(count_stmt)).one()` |

### Rule 3 — SELECT aggregates (`func.count()`)

```python
# Before
total_result = await session.execute(select(func.count()).select_from(Worker))
total = total_result.scalar()

# After
total_result = await session.exec(select(func.count()).select_from(Worker))
total = total_result.one()
```

`.one()` on a count query returns the integer directly.

### Rule 4 — Bulk `update()` / `delete()`

There is exactly one site in the codebase (`src/zndraw_joblib/router.py:733`) that uses a bulk `update(Task).where(...).values(...)` statement for optimistic task-claiming. SQLModel's `exec()` accepts `UpdateBase` statements natively via the third overload; the return type is `CursorResult[Any]` with `.rowcount` intact.

```python
# Before
stmt = update(Task).where(Task.id == task_id, Task.status == TaskStatus.PENDING).values(...)
cursor_result = await session.execute(stmt)
if cursor_result.rowcount == 1:
    ...

# After
stmt = update(Task).where(Task.id == task_id, Task.status == TaskStatus.PENDING).values(...)
cursor_result = await session.exec(stmt)
if cursor_result.rowcount == 1:
    ...
```

No type ignore, no `# noqa`, no SQLAlchemy/SQLModel mixing. This preserves the atomic optimistic-locking semantics the task-claim path depends on.

### Rule 5 — Iteration patterns

Where code iterates over a result:

```python
# Before
for row in (await session.execute(select(Model))).scalars():
    ...

# After
for row in await session.exec(select(Model)):
    ...
```

## Regression guard

Add to `[tool.pytest.ini_options]` in `pyproject.toml`:

```toml
filterwarnings = [
    'error:[\s\S]*You probably want to use .session\.exec:DeprecationWarning',
]
```

Notes:

- **Why message-based, not module-based.** `typing_extensions.@deprecated` uses `stacklevel=2`, so the warning's `module` attribute points to the calling code (e.g. `zndraw_joblib.registry`), not to `sqlmodel.*`. A `:sqlmodel.*:` filter would silently match nothing. Message-regex is the only reliable scoping.
- **Why `[\s\S]*` at the start.** The SQLModel deprecation message is a triple-quoted string that begins with `\n    🚨 You probably want to use...`. Python's `warnings.filterwarnings` compiles the regex without DOTALL, so `.` doesn't cross newlines. `[\s\S]*` matches the leading whitespace/emoji without needing to literal-match them.
- **Why TOML single-quoted.** Avoids double-backslashing the regex.
- **Effect.** Any reintroduction of `session.execute()` anywhere in `src/` (or in tests) becomes a hard CI failure with a clear message pointing at the culprit.

## Verification

- Run the full test suite via `uv run pytest tests/`.
- The existing suite already exercises every migrated code path — the 15k warning count is proof of coverage.
- Expected outcome: warning count goes from 15,000+ to 0.
- The `filterwarnings = error` gate converts any missed site into a hard test failure, pointing at the exact file and line.
- No new unit tests. The migration is mechanical and adding per-file tests would duplicate existing coverage.

## Risks & mitigations

| Risk | Mitigation |
|---|---|
| Missed call site silently slips through | `filterwarnings = error` converts it to a hard test failure |
| Aggregate `.one()` misbehaves on a count query | Admin/router count endpoints are under test; any divergence surfaces as a failing test |
| Complex selects (joins, subqueries, `select_from`, window functions) interact badly with sqlmodel's `select` | `sqlmodel.select` is a thin wrapper over sqlalchemy's — joins/subqueries/aggregates pass through identically per source inspection |
| Missing `and_`/`or_`/`update` imports after the `select` swap | Grep each file post-edit; the existing imports of these helpers stay put on the sqlalchemy import line |
| Test-order dependency caused by warning escalation | Warning is emitted synchronously on each call; no cross-test state |

## Execution shape

- Single branch: `fix/sqlmodel-session-exec-migration` in a dedicated worktree at `.claude/worktrees/fix+sqlmodel-session-exec-migration`.
- **Commit shape:** per-file commits (17), each `refactor(<pkg>): migrate <file> to session.exec()`, plus one final `chore(ci): error on sqlmodel session.execute deprecation warning`. Keeps bisection simple if something regresses.
- Final step: `uv run pytest tests/` must pass cleanly, then open a PR against `main`.

## Non-goals

- Touching test code (tests may still use `session.execute()` — the `filterwarnings` filter catches that too, so any test regressions show up naturally).
- Converting any `session.execute()` that lives outside `src/`.
- Migrating ORM-style `session.add()` / `session.delete()` / `session.get()` calls (these are already the SOTA pattern).
- Adding new tests — the migration is mechanical and existing coverage is adequate.
