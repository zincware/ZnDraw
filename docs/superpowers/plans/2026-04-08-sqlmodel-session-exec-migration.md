# SQLModel session.exec() Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate 15,000+ SQLModel `DeprecationWarning` instances from CI by migrating all 77 `session.execute()` call sites in `src/` to `session.exec()`, and add a pytest warning filter that prevents regressions.

**Architecture:** Mechanical, deterministic rewrite across 17 files. Each file becomes an independent atomic commit so bisection is trivial if a regression slips through. A single final commit adds a `pyproject.toml` `filterwarnings` entry that escalates any future reintroduction of the deprecated API to a hard test failure. No behavioral changes, no new code, no new tests.

**Tech Stack:** Python 3.11, SQLModel 0.0.24+, SQLAlchemy 2.x async, pytest, uv.

**Related spec:** `docs/superpowers/specs/2026-04-08-sqlmodel-session-exec-migration-design.md`
**Related issue:** [zincware/ZnDraw#890](https://github.com/zincware/ZnDraw/issues/890)
**Worktree:** `/Users/fzills/tools/zndraw-fastapi/.claude/worktrees/fix+sqlmodel-session-exec-migration`

---

## Rewrite rules (apply to every task)

These rules are derived from SQLModel's `AsyncSession.exec()` overloads in `.venv/lib/python3.11/site-packages/sqlmodel/ext/asyncio/session.py` (verified directly). They are the SOTA recommendation per sqlmodel.tiangolo.com.

### Rule 1 — Imports

Swap `select` to the SQLModel re-export. Keep `func`, `and_`, `or_`, `update`, `delete`, `text`, `selectinload`, etc. from sqlalchemy — SQLModel doesn't re-export all of them, and they interop transparently.

```python
# Before
from sqlalchemy import func, select, update

# After
from sqlalchemy import func, update
from sqlmodel import select
```

If a file already imports `select` from `sqlmodel`, do nothing to imports.

### Rule 2 — SELECT-one returning a model instance

| Before | After |
|---|---|
| `result = await session.execute(select(Model).where(...))`<br>`obj = result.scalar_one_or_none()` | `result = await session.exec(select(Model).where(...))`<br>`obj = result.one_or_none()` |
| `result = await session.execute(select(Model).where(...))`<br>`obj = result.scalar_one()` | `result = await session.exec(select(Model).where(...))`<br>`obj = result.one()` |

One-liner variants:

| Before | After |
|---|---|
| `obj = (await session.execute(stmt)).scalar_one_or_none()` | `obj = (await session.exec(stmt)).one_or_none()` |
| `obj = (await session.execute(stmt)).scalar_one()` | `obj = (await session.exec(stmt)).one()` |

### Rule 3 — SELECT-many returning model instances

| Before | After |
|---|---|
| `result = await session.execute(select(Model).where(...))`<br>`rows = result.scalars().all()` | `result = await session.exec(select(Model).where(...))`<br>`rows = result.all()` |
| `rows = (await session.execute(stmt)).scalars().all()` | `rows = (await session.exec(stmt)).all()` |
| `rows = list(result.scalars().all())` | `rows = list(result.all())` |

### Rule 4 — Iteration over scalars

| Before | After |
|---|---|
| `for row in result.scalars().all(): ...` | `for row in result.all(): ...` |
| `for row in result.scalars(): ...` | `for row in result: ...` |

### Rule 5 — SELECT aggregates (`func.count()`)

```python
# Before
total_result = await session.execute(select(func.count()).select_from(X))
total = total_result.scalar()

# After
total_result = await session.exec(select(func.count()).select_from(X))
total = total_result.one()
```

`.one()` on an aggregate result returns the scalar int directly.

### Rule 6 — Bulk `update()` / `delete()`

Per SQLModel source, `AsyncSession.exec()` has an `UpdateBase` overload returning `CursorResult[Any]` — bulk updates/deletes are fully typed. No type-ignore needed.

```python
# Before
stmt = update(Task).where(...).values(...)
cursor_result = await session.execute(stmt)
if cursor_result.rowcount == 1:
    ...

# After
stmt = update(Task).where(...).values(...)
cursor_result = await session.exec(stmt)
if cursor_result.rowcount == 1:
    ...
```

### Rule 7 — Multi-column SELECT (tuple result)

`select(Model.col1, Model.col2)` returns a `TupleResult` from `session.exec()`. Rows are tuple-like with attribute access by column key — the consuming code rarely needs to change.

```python
# Before
result = await session.execute(select(Job.id, Job.room_id).where(...))
job_rooms = {row.id: row.room_id for row in result.all()}

# After
result = await session.exec(select(Job.id, Job.room_id).where(...))
job_rooms = {row.id: row.room_id for row in result.all()}
```

Note: do **not** call `.scalars()` on a multi-column select — it would drop all columns but the first.

---

## Special notes

- **sqlmodel re-exports `func`** (`from sqlmodel import func, select` works), but most files in this codebase import `func` from sqlalchemy. Keep the existing style — don't churn imports beyond the `select` swap.
- **`src/zndraw/routes/frames.py`** uses `from sqlalchemy import select as sa_select` and references `sa_select(...)` at line 88. There is only one call site in the file. Replace the alias with a normal `from sqlmodel import select` and drop the `sa_` prefix.
- **`src/zndraw_joblib/registry.py`** has a **local** `from sqlalchemy import select` inside a function (line 98), not a module-level import. Change the local import to `from sqlmodel import select`.
- **`src/zndraw/routes/selection_groups.py`** similarly has a local `from sqlmodel import select` inside a function (line 48). It already imports from sqlmodel — only the `session.execute` → `session.exec` + `.scalars()` swap is needed.
- **`src/zndraw/routes/rooms.py`** uses a pattern where `for row in result.scalars().all(): ...` appears four times (lines 205, 223, 236, 250). Apply Rule 4.
- **`src/zndraw_joblib/router.py:730`** is the one bulk `update()` call site in the entire `src/` tree. Rule 6 applies.
- **`src/zndraw_joblib/sweeper.py:127-130`** is the one multi-column select. Rule 7 applies — the consuming dict comprehension stays identical.

---

## File structure

| File | Sites | Rules invoked | Notes |
|---|---|---|---|
| `src/zndraw_joblib/registry.py` | 1 | R1 (local), R2 | Local import inside function |
| `src/zndraw_joblib/dependencies.py` | 1 | R2 | Already imports sqlmodel select |
| `src/zndraw_joblib/sweeper.py` | 9 | R1, R2, R3, R4, R7 | Multi-column select at L127-130 |
| `src/zndraw_joblib/router.py` | 40 | R1, R2, R3, R5, R6 | The big file, incl. bulk update |
| `src/zndraw/database.py` | 1 | R2 | Already imports sqlmodel select |
| `src/zndraw/routes/admin.py` | 2 | R2, R3, R5 | Already imports sqlmodel select+func |
| `src/zndraw/routes/rooms.py` | 5 | R3, R4 | 4× iteration pattern |
| `src/zndraw/routes/presets.py` | 4 | R3 | One-liner patterns |
| `src/zndraw/routes/chat.py` | 3 | R3, R5 | Mixed select + count |
| `src/zndraw/routes/screenshots.py` | 2 | R3, R5 | One-liner count + select |
| `src/zndraw/routes/geometries.py` | 1 | R3 | |
| `src/zndraw/routes/figures.py` | 1 | R3 | |
| `src/zndraw/routes/selection_groups.py` | 1 | R4 | Local sqlmodel select already |
| `src/zndraw/routes/bookmarks.py` | 1 | R3 | |
| `src/zndraw/routes/frames.py` | 1 | R1, R2 | Drop `sa_select` alias |
| `src/zndraw_auth/db.py` | 1 | R1, R2 | |
| `src/zndraw_auth/cli_login.py` | 3 | R1, R2 | All R2 pattern |

**Total: 77 sites across 17 files.**

Final change: `pyproject.toml` — add `filterwarnings` entry under `[tool.pytest.ini_options]`.

---

## Task list

### Task 1: Baseline verification

**Purpose:** Confirm the working tree matches the spec's assumptions before migrating anything. Establishes green baseline for bisection.

**Files:**
- Read: `pyproject.toml` (verify no existing `filterwarnings`)
- Read: `src/zndraw/socketio.py:200` (existing `session.exec()` reference — sanity check)

- [ ] **Step 1: Verify worktree and branch**

Run:
```bash
pwd
git status --short
git branch --show-current
```

Expected: pwd ends in `.claude/worktrees/fix+sqlmodel-session-exec-migration`, branch is `worktree-fix+sqlmodel-session-exec-migration`, working tree clean.

- [ ] **Step 2: Confirm call site count matches the spec**

Run:
```bash
grep -rc "session\.execute(" src/ | grep -v ":0$" | sort
```

Expected output (17 files, summing to 77):
```
src/zndraw/database.py:1
src/zndraw/routes/admin.py:2
src/zndraw/routes/bookmarks.py:1
src/zndraw/routes/chat.py:3
src/zndraw/routes/figures.py:1
src/zndraw/routes/frames.py:1
src/zndraw/routes/geometries.py:1
src/zndraw/routes/presets.py:4
src/zndraw/routes/rooms.py:5
src/zndraw/routes/screenshots.py:2
src/zndraw/routes/selection_groups.py:1
src/zndraw_auth/cli_login.py:3
src/zndraw_auth/db.py:1
src/zndraw_joblib/dependencies.py:1
src/zndraw_joblib/registry.py:1
src/zndraw_joblib/router.py:40
src/zndraw_joblib/sweeper.py:9
```

If any count differs, stop and reconcile against the spec before proceeding.

- [ ] **Step 3: Run the full existing test suite as a green baseline**

Run:
```bash
uv sync
uv run --active pytest tests/ -q 2>&1 | tail -20
```

Expected: all tests pass. Note the "X passed" and "Y warnings" counts for later comparison.

If tests are failing on `main`, stop — the migration's verification step won't distinguish new failures from pre-existing ones.

---

### Task 2: Migrate `src/zndraw_joblib/registry.py` (1 site — worked example)

**Files:**
- Modify: `src/zndraw_joblib/registry.py:98,108,115`

This is the canonical single-site migration. Use it as a template for understanding the other tasks.

- [ ] **Step 1: Read the current state of the function**

Read lines 95-125 of `src/zndraw_joblib/registry.py`. The import is local inside `register_internal_jobs()` at line 98. The `session.execute()` call is at line 108 with `scalar_one_or_none()` at line 115.

- [ ] **Step 2: Swap the local import**

Change line 98:

```python
# Before
    from sqlalchemy import select

# After
    from sqlmodel import select
```

- [ ] **Step 3: Swap `session.execute` → `session.exec` and `scalar_one_or_none` → `one_or_none`**

Change line 108 and line 115:

```python
# Before
            result = await session.execute(
                select(Job).where(
                    Job.room_id == "@internal",
                    Job.category == category,
                    Job.name == name,
                )
            )
            existing = result.scalar_one_or_none()

# After
            result = await session.exec(
                select(Job).where(
                    Job.room_id == "@internal",
                    Job.category == category,
                    Job.name == name,
                )
            )
            existing = result.one_or_none()
```

- [ ] **Step 4: Verify no `session.execute` or `scalar_one_or_none` remains in the file**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_joblib/registry.py
```

Expected: no output.

- [ ] **Step 5: Run the test that covers registry.py**

Run:
```bash
uv run --active pytest tests/zndraw_joblib/test_registry.py -q
```

Expected: `8 passed`.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_joblib/registry.py
git commit -m "$(cat <<'EOF'
refactor(joblib): migrate registry.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Migrate `src/zndraw_joblib/dependencies.py` (1 site)

**Files:**
- Modify: `src/zndraw_joblib/dependencies.py:153,156`

Module already imports `from sqlmodel import select` (line 8) — no import change.

- [ ] **Step 1: Read lines 150-165 of `src/zndraw_joblib/dependencies.py`**

Confirm: line 153 calls `session.execute(...)`, line 156 calls `result.scalar_one_or_none()`.

- [ ] **Step 2: Rewrite the call**

```python
# Before (lines 153-156)
    result = await session.execute(
        select(User).where(User.email == settings.internal_worker_email)  # type: ignore[arg-type]
    )
    user = result.scalar_one_or_none()

# After
    result = await session.exec(
        select(User).where(User.email == settings.internal_worker_email)  # type: ignore[arg-type]
    )
    user = result.one_or_none()
```

- [ ] **Step 3: Verify the file is clean**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_joblib/dependencies.py
```

Expected: no output.

- [ ] **Step 4: Run relevant tests**

Run:
```bash
uv run --active pytest tests/zndraw_joblib/ -q -x
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_joblib/dependencies.py
git commit -m "$(cat <<'EOF'
refactor(joblib): migrate dependencies.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Migrate `src/zndraw_joblib/sweeper.py` (9 sites)

**Files:**
- Modify: `src/zndraw_joblib/sweeper.py` — line 11 (import), and all 9 call sites.

Call sites map (from `grep -n "session\.execute"`):

| Line | Pattern | Rule |
|---|---|---|
| 49 | `session.execute(select(Job).where(Job.id == job_id))` → `.scalar_one_or_none()` | R2 |
| 55 | `session.execute(select(WorkerJobLink)...limit(1))` → `if result.scalar_one_or_none()` | R2 |
| 62 | `session.execute(select(Task)...limit(1))` → `if result.scalar_one_or_none()` | R2 |
| 100 | `session.execute(select(Task).options(selectinload)...)` → `.scalars().all()` | R3 |
| 119 | `session.execute(select(WorkerJobLink).where(...))` → `.scalars().all()` | R3 |
| 127 | `session.execute(select(Job.id, Job.room_id).where(...))` → `result.all()` for tuple rows | R7 |
| 139 | `session.execute(select(ProviderRecord).where(...))` → `.scalars().all()` | R3 |
| 191 | `session.execute(select(Worker).where(...))` → `.scalars().all()` | R3 |
| 229 | `session.execute(select(Task).join(Job)...)` → `.scalars().all()` | R3 |

- [ ] **Step 1: Swap the import**

Change line 11:

```python
# Before
from sqlalchemy import func as sa_func, select

# After
from sqlalchemy import func as sa_func
from sqlmodel import select
```

Keep `sa_func` as-is — it's used at line 236 and renaming would churn an unrelated line.

- [ ] **Step 2: Rewrite line 49**

```python
# Before
    result = await session.execute(select(Job).where(Job.id == job_id))
    job = result.scalar_one_or_none()

# After
    result = await session.exec(select(Job).where(Job.id == job_id))
    job = result.one_or_none()
```

- [ ] **Step 3: Rewrite lines 55-58**

```python
# Before
    result = await session.execute(
        select(WorkerJobLink).where(WorkerJobLink.job_id == job_id).limit(1)
    )
    if result.scalar_one_or_none():

# After
    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.job_id == job_id).limit(1)
    )
    if result.one_or_none():
```

- [ ] **Step 4: Rewrite lines 62-70**

```python
# Before
    result = await session.execute(
        select(Task)
        .where(
            Task.job_id == job_id,
            Task.status == TaskStatus.PENDING,
        )
        .limit(1)
    )
    if result.scalar_one_or_none():

# After
    result = await session.exec(
        select(Task)
        .where(
            Task.job_id == job_id,
            Task.status == TaskStatus.PENDING,
        )
        .limit(1)
    )
    if result.one_or_none():
```

- [ ] **Step 5: Rewrite lines 100-108**

```python
# Before
    result = await session.execute(
        select(Task)
        .options(selectinload(Task.job))
        .where(
            Task.worker_id == worker.id,
            Task.status.in_({TaskStatus.CLAIMED, TaskStatus.RUNNING}),
        )
    )
    worker_tasks = result.scalars().all()

# After
    result = await session.exec(
        select(Task)
        .options(selectinload(Task.job))
        .where(
            Task.worker_id == worker.id,
            Task.status.in_({TaskStatus.CLAIMED, TaskStatus.RUNNING}),
        )
    )
    worker_tasks = result.all()
```

- [ ] **Step 6: Rewrite lines 119-122**

```python
# Before
    result = await session.execute(
        select(WorkerJobLink).where(WorkerJobLink.worker_id == worker.id)
    )
    links = result.scalars().all()

# After
    result = await session.exec(
        select(WorkerJobLink).where(WorkerJobLink.worker_id == worker.id)
    )
    links = result.all()
```

- [ ] **Step 7: Rewrite lines 127-130 (multi-column select — Rule 7)**

Only `execute` → `exec` changes. The dict comprehension stays identical — tuple rows keep attribute access.

```python
# Before
        result = await session.execute(
            select(Job.id, Job.room_id).where(Job.id.in_(job_ids))
        )
        job_rooms = {row.id: row.room_id for row in result.all()}

# After
        result = await session.exec(
            select(Job.id, Job.room_id).where(Job.id.in_(job_ids))
        )
        job_rooms = {row.id: row.room_id for row in result.all()}
```

- [ ] **Step 8: Rewrite lines 139-142**

```python
# Before
    result = await session.execute(
        select(ProviderRecord).where(ProviderRecord.worker_id == worker.id)
    )
    providers = result.scalars().all()

# After
    result = await session.exec(
        select(ProviderRecord).where(ProviderRecord.worker_id == worker.id)
    )
    providers = result.all()
```

- [ ] **Step 9: Rewrite lines 191-192**

```python
# Before
    result = await session.execute(select(Worker).where(Worker.last_heartbeat < cutoff))
    stale_workers = result.scalars().all()

# After
    result = await session.exec(select(Worker).where(Worker.last_heartbeat < cutoff))
    stale_workers = result.all()
```

- [ ] **Step 10: Rewrite lines 229-239**

```python
# Before
    result = await session.execute(
        select(Task)
        .join(Job)
        .options(selectinload(Task.job))
        .where(
            Job.room_id == "@internal",
            Task.status.in_({TaskStatus.RUNNING, TaskStatus.CLAIMED}),
            sa_func.coalesce(Task.started_at, Task.created_at) < cutoff,
        )
    )
    stuck_tasks = result.scalars().all()

# After
    result = await session.exec(
        select(Task)
        .join(Job)
        .options(selectinload(Task.job))
        .where(
            Job.room_id == "@internal",
            Task.status.in_({TaskStatus.RUNNING, TaskStatus.CLAIMED}),
            sa_func.coalesce(Task.started_at, Task.created_at) < cutoff,
        )
    )
    stuck_tasks = result.all()
```

- [ ] **Step 11: Verify the file is clean**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_joblib/sweeper.py
```

Expected: no output.

- [ ] **Step 12: Run sweeper tests**

Run:
```bash
uv run --active pytest tests/zndraw_joblib/test_sweeper.py -q
```

Expected: `16 passed`.

- [ ] **Step 13: Commit**

```bash
git add src/zndraw_joblib/sweeper.py
git commit -m "$(cat <<'EOF'
refactor(joblib): migrate sweeper.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Migrate `src/zndraw_joblib/router.py` (40 sites)

**Files:**
- Modify: `src/zndraw_joblib/router.py` — line 11 (import), plus 40 call sites across ~25 functions.

This is the largest file. Call sites are grouped by function. Apply the rules rule-by-rule; most sites are Rule 2 or Rule 3 mechanical rewrites. **One bulk update site at line 733 uses Rule 6. Four count aggregates use Rule 5.**

Full call site map (from `grep -n`):

| Line | Function | Rule | Pattern |
|---|---|---|---|
| 137 | `_resolve_job` | R2 | `.scalar_one_or_none()` |
| 151 | `_task_response` | R2 | `.scalar_one_or_none()` |
| 184 | `_bulk_queue_positions` | R3 | `.scalars().all()` (iterated) |
| 223 | `_task_status_emission` | R2 | `.scalar_one_or_none()` |
| 294 | `register_job` | R2 | `.scalar_one_or_none()` |
| 331 | `register_job` | R2 | `.scalar_one_or_none()` |
| 349 | `register_job` | R2 | `.scalar_one_or_none()` |
| 365 | `register_job` | R3 | `.scalars().all()` |
| 396 | `list_jobs` | R5 | `select(...).scalar()` count |
| 402 | `list_jobs` | R3 | `.scalars().all()` |
| 432 | `list_workers_for_room` | R3 | `.scalars().all()` |
| 448 | `list_workers_for_room` | R5 | count |
| 454 | `list_workers_for_room` | R3 | `.scalars().all()` |
| 494 | `list_tasks_for_room` | R5 | count |
| 500 | `list_tasks_for_room` | R3 | `.scalars().all()` |
| 534 | `list_tasks_for_job` | R5 | count |
| 540 | `list_tasks_for_job` | R3 | `.scalars().all()` |
| 563 | `get_job` | R3 | `.scalars().all()` |
| 613 | `submit_task` | R5 | one-liner `.scalar_one()` count |
| 687 | `claim_task` | R2 | `.scalar_one_or_none()` |
| 697 | `claim_task` | R3 | `.scalars().all()` |
| 714 | `claim_task` | R2 | `.scalar_one_or_none()` (inside retry loop) |
| **733** | `claim_task` | **R6** | **bulk `update(Task)` — keep `.rowcount`** |
| 765 | `claim_task` | R2 | `.scalar_one()` |
| 784 | `get_task_status` | R2 | `.scalar_one_or_none()` |
| 805 | `get_task_status` | R2 | `.scalar_one_or_none()` (inside poll loop) |
| 814 | `get_task_status` | R2 | `.scalar_one()` |
| 828 | `update_task_status` | R2 | `.scalar_one_or_none()` |
| 837 | `update_task_status` | R2 | `.scalar_one_or_none()` |
| 891 | `list_workers` | R5 | one-liner-ish count |
| 895 | `list_workers` | R3 | `.scalars().all()` |
| 922 | `worker_heartbeat` | R2 | `.scalar_one_or_none()` |
| 950 | `delete_worker` | R2 | `.scalar_one_or_none()` |
| 997 | `_resolve_provider` | R2 | `.scalar_one_or_none()` |
| 1043 | `register_provider` | R2 | `.scalar_one_or_none()` |
| 1060 | `register_provider` | R2 | `.scalar_one_or_none()` |
| 1114 | `list_providers` | R5 | count |
| 1119 | `list_providers` | R3 | `.scalars().all()` |
| 1226 | `delete_provider` | R2 | `.scalar_one_or_none()` |
| 1258 | `upload_provider_result` | R2 | `.scalar_one_or_none()` |

- [ ] **Step 1: Swap the import**

Change line 11:

```python
# Before
from sqlalchemy import func, select, update

# After
from sqlalchemy import func, update
from sqlmodel import select
```

- [ ] **Step 2: Apply Rule 2 rewrites (SELECT-one)**

For each of these sites, change `session.execute` → `session.exec` and the result accessor as shown. Use `Edit` with enough context to make each change unique.

Lines 137 + 144, 151 + 152, 223 + 224, 294 + 301, 331 + 337, 349 + 355, 687 + 688, 714 + 722, 765 + 766, 784 + 785, 805 + 806, 814 + 815, 828 + 829, 837 + 840, 922 + 923, 950 + 951, 997 + 1004, 1043 + 1049, 1060 + 1067, 1226 + 1229, 1258 + 1261:

Rewrite each:
- `await session.execute(` → `await session.exec(`
- `result.scalar_one_or_none()` → `result.one_or_none()`
- `result.scalar_one()` → `result.one()` (at lines 766 and 815)

- [ ] **Step 3: Apply Rule 3 rewrites (SELECT-many)**

Lines 184, 365, 402, 432, 454, 500, 540, 563, 697, 895, 1119:

Rewrite each:
- `await session.execute(` → `await session.exec(`
- `result.scalars().all()` → `result.all()`

- [ ] **Step 4: Apply Rule 5 rewrites (aggregates/count)**

Lines 396 + 399, 448 + 451, 494 + 497, 534 + 537, 891 + 892, 1114 + 1117:

Rewrite each:
- `total_result = await session.execute(` → `total_result = await session.exec(`
- `total = total_result.scalar()` → `total = total_result.one()`

Line 613-616 is a one-liner count in `submit_task`:

```python
# Before
        worker_count = (
            await session.execute(
                select(func.count()).where(WorkerJobLink.job_id == job.id)
            )
        ).scalar_one()

# After
        worker_count = (
            await session.exec(
                select(func.count()).where(WorkerJobLink.job_id == job.id)
            )
        ).one()
```

Line 891 is also one-liner style:

```python
# Before
    total_result = await session.execute(select(func.count()).select_from(Worker))
    total = total_result.scalar()

# After
    total_result = await session.exec(select(func.count()).select_from(Worker))
    total = total_result.one()
```

- [ ] **Step 5: Apply Rule 6 rewrite (bulk UPDATE at line 733)**

This is the optimistic-locking task-claim path. The `CursorResult.rowcount` semantics must be preserved.

```python
# Before (lines 727-737)
            # Atomically update only if still PENDING (optimistic locking)
            stmt = (
                update(Task)
                .where(Task.id == task_id, Task.status == TaskStatus.PENDING)
                .values(status=TaskStatus.CLAIMED, worker_id=request.worker_id)
            )
            cursor_result = await session.execute(stmt)
            await session.commit()

            # Check if we actually claimed it (rowcount == 1 means success)
            if cursor_result.rowcount == 1:

# After
            # Atomically update only if still PENDING (optimistic locking)
            stmt = (
                update(Task)
                .where(Task.id == task_id, Task.status == TaskStatus.PENDING)
                .values(status=TaskStatus.CLAIMED, worker_id=request.worker_id)
            )
            cursor_result = await session.exec(stmt)
            await session.commit()

            # Check if we actually claimed it (rowcount == 1 means success)
            if cursor_result.rowcount == 1:
```

SQLModel's `exec()` has an `UpdateBase` overload returning `CursorResult[Any]`, so `.rowcount` is fully typed. No ignore comment needed.

- [ ] **Step 6: Verify the file is clean**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_joblib/router.py
```

Expected: no output.

- [ ] **Step 7: Run the full joblib test suite**

Run:
```bash
uv run --active pytest tests/zndraw_joblib/ -q
```

Expected: all tests pass. Pay special attention to `test_resilience.py` which exercises the claim-task / optimistic-locking path.

- [ ] **Step 8: Commit**

```bash
git add src/zndraw_joblib/router.py
git commit -m "$(cat <<'EOF'
refactor(joblib): migrate router.py to session.exec() (#890)

40 sites across ~25 endpoints, including the bulk update() in the
optimistic-locking claim_task path (uses the UpdateBase overload on
AsyncSession.exec for fully typed CursorResult.rowcount).

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Migrate `src/zndraw/database.py` (1 site)

**Files:**
- Modify: `src/zndraw/database.py:113,116`

Module already imports `select` from sqlmodel (line 24) — no import change.

- [ ] **Step 1: Read lines 108-125 of `src/zndraw/database.py`**

- [ ] **Step 2: Rewrite the call**

```python
# Before (lines 113-116)
    result = await session.execute(
        select(User).where(User.email == internal_worker_email)
    )
    existing = result.scalar_one_or_none()

# After
    result = await session.exec(
        select(User).where(User.email == internal_worker_email)
    )
    existing = result.one_or_none()
```

Note: the exact `.where(...)` expression may differ — read the file and preserve it verbatim. Only change `execute` → `exec` and `scalar_one_or_none` → `one_or_none`.

- [ ] **Step 3: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/database.py
```

Expected: no output.

- [ ] **Step 4: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/ -q -x -k "database or room"
```

Expected: tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/database.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate database.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Migrate `src/zndraw/routes/admin.py` (2 sites)

**Files:**
- Modify: `src/zndraw/routes/admin.py:75,76,79,80`

Module already imports `from sqlmodel import func, select` (line 14) — no import change.

- [ ] **Step 1: Rewrite lines 75-80**

```python
# Before
    count_result = await session.execute(select(func.count()).select_from(User))
    total = count_result.scalar_one()

    # Get paginated users
    result = await session.execute(select(User).offset(offset).limit(limit))
    users = list(result.scalars().all())

# After
    count_result = await session.exec(select(func.count()).select_from(User))
    total = count_result.one()

    # Get paginated users
    result = await session.exec(select(User).offset(offset).limit(limit))
    users = list(result.all())
```

- [ ] **Step 2: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/admin.py
```

Expected: no output.

- [ ] **Step 3: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_admin.py -q
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/admin.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/admin.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Migrate `src/zndraw/routes/rooms.py` (5 sites)

**Files:**
- Modify: `src/zndraw/routes/rooms.py:202-205, 220-223, 233-236, 247-250, 420-421`

Module already imports `from sqlmodel import select` (line 14) — no import change.

All four sites at 202-250 use the same iterate-over-scalars pattern (Rule 3 + Rule 4 combined). Site at 420 uses a list-wrap pattern.

- [ ] **Step 1: Rewrite lines 202-205**

```python
# Before
    result = await session.execute(
        select(RoomGeometry).where(RoomGeometry.room_id == source_room_id)
    )
    for row in result.scalars().all():

# After
    result = await session.exec(
        select(RoomGeometry).where(RoomGeometry.room_id == source_room_id)
    )
    for row in result.all():
```

- [ ] **Step 2: Rewrite lines 220-223**

```python
# Before
    result = await session.execute(
        select(RoomBookmark).where(RoomBookmark.room_id == source_room_id)
    )
    for row in result.scalars().all():

# After
    result = await session.exec(
        select(RoomBookmark).where(RoomBookmark.room_id == source_room_id)
    )
    for row in result.all():
```

- [ ] **Step 3: Rewrite lines 233-236**

```python
# Before
    result = await session.execute(
        select(RoomFigure).where(RoomFigure.room_id == source_room_id)
    )
    for row in result.scalars().all():

# After
    result = await session.exec(
        select(RoomFigure).where(RoomFigure.room_id == source_room_id)
    )
    for row in result.all():
```

- [ ] **Step 4: Rewrite lines 247-250**

Read the 4 lines before the call to capture the exact model reference — it's one of `SelectionGroup`, `RoomFigure`, `RoomBookmark`, `RoomGeometry`. Apply the same pattern:

```python
# Before
    result = await session.execute(
        select(<Model>).where(<Model>.room_id == source_room_id)
    )
    for row in result.scalars().all():

# After
    result = await session.exec(
        select(<Model>).where(<Model>.room_id == source_room_id)
    )
    for row in result.all():
```

- [ ] **Step 5: Rewrite lines 420-421**

```python
# Before
    result = await session.execute(statement)
    rooms = list(result.scalars().all())

# After
    result = await session.exec(statement)
    rooms = list(result.all())
```

- [ ] **Step 6: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/rooms.py
```

Expected: no output.

- [ ] **Step 7: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_admin.py tests/zndraw/test_rooms.py -q 2>&1 | tail -10
```

Expected: all tests pass. (test_admin covers room-copy via `test_create_room_uses_default_when_no_copy_from`.)

- [ ] **Step 8: Commit**

```bash
git add src/zndraw/routes/rooms.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/rooms.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 9: Migrate `src/zndraw/routes/presets.py` (4 sites)

**Files:**
- Modify: `src/zndraw/routes/presets.py:100-103, 260, 269, 288`

Module already imports `from sqlmodel import select` (line 9) — no import change.

- [ ] **Step 1: Rewrite lines 100-103**

```python
# Before
    result = await session.execute(
        <statement>
    )
    rows = result.scalars().all()

# After
    result = await session.exec(
        <statement>
    )
    rows = result.all()
```

Read the file to capture the exact statement — only change `execute` → `exec` and `result.scalars().all()` → `result.all()`.

- [ ] **Step 2: Rewrite line 260 (one-liner)**

```python
# Before
        existing = (await session.execute(stmt)).scalars().all()

# After
        existing = (await session.exec(stmt)).all()
```

- [ ] **Step 3: Rewrite line 269 (one-liner inside comprehension)**

```python
# Before
        new_keys = [g.key for g in (await session.execute(new_stmt)).scalars().all()]

# After
        new_keys = [g.key for g in (await session.exec(new_stmt)).all()]
```

- [ ] **Step 4: Rewrite line 288 (one-liner)**

```python
# Before
    geometries = (await session.execute(stmt)).scalars().all()

# After
    geometries = (await session.exec(stmt)).all()
```

- [ ] **Step 5: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/presets.py
```

Expected: no output.

- [ ] **Step 6: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_presets.py -q
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/zndraw/routes/presets.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/presets.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Migrate `src/zndraw/routes/chat.py` (3 sites)

**Files:**
- Modify: `src/zndraw/routes/chat.py:78-79, 88, 94-97`

Module already imports `from sqlmodel import col, select` (line 8). Keep `col` unchanged.

- [ ] **Step 1: Rewrite lines 78-79**

```python
# Before
    result = await session.execute(stmt)
    rows = list(result.scalars().all())

# After
    result = await session.exec(stmt)
    rows = list(result.all())
```

- [ ] **Step 2: Rewrite line 88 (count one-liner)**

```python
# Before
    total_count = (await session.execute(count_stmt)).scalar_one()

# After
    total_count = (await session.exec(count_stmt)).one()
```

- [ ] **Step 3: Rewrite lines 94-97**

```python
# Before
        users_result = await session.execute(
            <statement>
        )
        for user in users_result.scalars().all():

# After
        users_result = await session.exec(
            <statement>
        )
        for user in users_result.all():
```

Read the file to preserve the exact `<statement>` text.

- [ ] **Step 4: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/chat.py
```

Expected: no output.

- [ ] **Step 5: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_chat.py -q 2>&1 | tail -10
```

Expected: all tests pass (or "no tests ran" if no chat tests exist — still acceptable).

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/routes/chat.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/chat.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 11: Migrate `src/zndraw/routes/screenshots.py` (2 sites)

**Files:**
- Modify: `src/zndraw/routes/screenshots.py:240-241, 248`

Module already imports `from sqlmodel import col, select` (line 10). No import change.

- [ ] **Step 1: Rewrite lines 240-241**

```python
# Before
    result = await session.execute(stmt)
    rows = list(result.scalars().all())

# After
    result = await session.exec(stmt)
    rows = list(result.all())
```

- [ ] **Step 2: Rewrite line 248**

```python
# Before
    total = (await session.execute(count_stmt)).scalar_one()

# After
    total = (await session.exec(count_stmt)).one()
```

- [ ] **Step 3: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/screenshots.py
```

Expected: no output.

- [ ] **Step 4: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_screenshots.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/screenshots.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/screenshots.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 12: Migrate `src/zndraw/routes/geometries.py` (1 site)

**Files:**
- Modify: `src/zndraw/routes/geometries.py:107-110`

Module already imports `from sqlmodel import select` (line 9). No import change.

- [ ] **Step 1: Rewrite lines 107-110**

```python
# Before
    result = await session.execute(
        <statement>
    )
    rows = result.scalars().all()

# After
    result = await session.exec(
        <statement>
    )
    rows = result.all()
```

Read the file to preserve the exact `<statement>`.

- [ ] **Step 2: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/geometries.py
```

Expected: no output.

- [ ] **Step 3: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_geometries.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/geometries.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/geometries.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 13: Migrate `src/zndraw/routes/figures.py` (1 site)

**Files:**
- Modify: `src/zndraw/routes/figures.py:46-49`

Module already imports `from sqlmodel import select` (line 4). No import change.

- [ ] **Step 1: Rewrite lines 46-49**

```python
# Before
    result = await session.execute(
        <statement>
    )
    keys = list(result.scalars().all())

# After
    result = await session.exec(
        <statement>
    )
    keys = list(result.all())
```

- [ ] **Step 2: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/figures.py
```

Expected: no output.

- [ ] **Step 3: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_figures.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/figures.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/figures.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 14: Migrate `src/zndraw/routes/selection_groups.py` (1 site)

**Files:**
- Modify: `src/zndraw/routes/selection_groups.py:50-54`

This file has a **local** `from sqlmodel import select` inside the function (line 48). No module-level import change.

- [ ] **Step 1: Rewrite lines 50-54**

```python
# Before
    result = await session.execute(
        select(SelectionGroup).where(SelectionGroup.room_id == room_id)
    )
    groups: dict[str, dict[str, list[int]]] = {}
    for row in result.scalars().all():

# After
    result = await session.exec(
        select(SelectionGroup).where(SelectionGroup.room_id == room_id)
    )
    groups: dict[str, dict[str, list[int]]] = {}
    for row in result.all():
```

- [ ] **Step 2: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/selection_groups.py
```

Expected: no output.

- [ ] **Step 3: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_selection_groups.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/selection_groups.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/selection_groups.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 15: Migrate `src/zndraw/routes/bookmarks.py` (1 site)

**Files:**
- Modify: `src/zndraw/routes/bookmarks.py:44-47`

Module already imports `from sqlmodel import select` (line 4). No import change.

- [ ] **Step 1: Rewrite lines 44-47**

```python
# Before
    result = await session.execute(
        <statement>
    )
    rows = result.scalars().all()

# After
    result = await session.exec(
        <statement>
    )
    rows = result.all()
```

Read the file to preserve the exact `<statement>`.

- [ ] **Step 2: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/bookmarks.py
```

Expected: no output.

- [ ] **Step 3: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_bookmarks.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/bookmarks.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/bookmarks.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 16: Migrate `src/zndraw/routes/frames.py` (1 site — drop `sa_select` alias)

**Files:**
- Modify: `src/zndraw/routes/frames.py:17, 87-93, 88`

This file has a `sa_select` alias. There is only one call site (line 88), so the alias can be removed entirely.

- [ ] **Step 1: Swap the import at line 17**

```python
# Before
from sqlalchemy import select as sa_select

# After
from sqlmodel import select
```

- [ ] **Step 2: Rewrite line 88**

```python
# Before
        sa_select(ProviderRecord).where(

# After
        select(ProviderRecord).where(
```

- [ ] **Step 3: Rewrite lines 87 + 93**

```python
# Before (lines 87-93)
    result = await session.execute(
        select(ProviderRecord).where(
            <filters>
        )
    )
    return result.scalar_one_or_none()

# After
    result = await session.exec(
        select(ProviderRecord).where(
            <filters>
        )
    )
    return result.one_or_none()
```

- [ ] **Step 4: Verify no `sa_select` or old patterns remain**

Run:
```bash
grep -n "sa_select\|session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw/routes/frames.py
```

Expected: no output.

- [ ] **Step 5: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw/test_frames.py -q 2>&1 | tail -10
```

Expected: tests pass (or "no tests ran").

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/routes/frames.py
git commit -m "$(cat <<'EOF'
refactor(zndraw): migrate routes/frames.py to session.exec() (#890)

Drops the sa_select alias — the only call site migrates to sqlmodel's
select via session.exec().

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 17: Migrate `src/zndraw_auth/db.py` (1 site)

**Files:**
- Modify: `src/zndraw_auth/db.py:12, 161-164`

Module imports `from sqlalchemy import select` (line 12). Swap to sqlmodel.

- [ ] **Step 1: Swap the import at line 12**

```python
# Before
from sqlalchemy import select

# After
from sqlmodel import select
```

- [ ] **Step 2: Rewrite lines 161-164**

```python
# Before
    result = await session.execute(
        <statement>
    )
    existing = result.scalar_one_or_none()

# After
    result = await session.exec(
        <statement>
    )
    existing = result.one_or_none()
```

Read the file to preserve the exact `<statement>`.

- [ ] **Step 3: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_auth/db.py
```

Expected: no output.

- [ ] **Step 4: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw_auth/ -q
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_auth/db.py
git commit -m "$(cat <<'EOF'
refactor(auth): migrate db.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 18: Migrate `src/zndraw_auth/cli_login.py` (3 sites)

**Files:**
- Modify: `src/zndraw_auth/cli_login.py:12, 68, 71, 121, 124, 155, 158`

Module imports `from sqlalchemy import select` (line 12). Swap to sqlmodel.

- [ ] **Step 1: Swap the import at line 12**

```python
# Before
from sqlalchemy import select

# After
from sqlmodel import select
```

- [ ] **Step 2: Rewrite lines 68-71**

```python
# Before
    result = await session.execute(
        <statement>
    )
    challenge = result.scalar_one_or_none()

# After
    result = await session.exec(
        <statement>
    )
    challenge = result.one_or_none()
```

- [ ] **Step 3: Rewrite lines 121-124**

Same pattern as Step 2 — different `<statement>`. Read the file to capture it verbatim.

```python
# Before
    result = await session.execute(
        <statement>
    )
    challenge = result.scalar_one_or_none()

# After
    result = await session.exec(
        <statement>
    )
    challenge = result.one_or_none()
```

- [ ] **Step 4: Rewrite lines 155-158**

Same pattern as Steps 2-3 — different `<statement>`.

```python
# Before
    result = await session.execute(
        <statement>
    )
    challenge = result.scalar_one_or_none()

# After
    result = await session.exec(
        <statement>
    )
    challenge = result.one_or_none()
```

- [ ] **Step 5: Verify**

Run:
```bash
grep -n "session\.execute\|scalar_one_or_none\|scalar_one\b\|\.scalars(" src/zndraw_auth/cli_login.py
```

Expected: no output.

- [ ] **Step 6: Run tests**

Run:
```bash
uv run --active pytest tests/zndraw_auth/ -q
```

Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/zndraw_auth/cli_login.py
git commit -m "$(cat <<'EOF'
refactor(auth): migrate cli_login.py to session.exec() (#890)

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 19: Add the `filterwarnings` gate

**Files:**
- Modify: `pyproject.toml:222-226` (`[tool.pytest.ini_options]` block)

Adds the regression guard. This commit lands **last** so every earlier commit bisects against a green baseline.

- [ ] **Step 1: Read the current `[tool.pytest.ini_options]` block**

Run:
```bash
grep -nA 10 "\[tool.pytest.ini_options\]" pyproject.toml
```

Expected: the current block has `testpaths`, `asyncio_mode`, `asyncio_default_fixture_loop_scope`, `markers`. No existing `filterwarnings`.

- [ ] **Step 2: Add the `filterwarnings` entry**

Append to the block (after the `markers` line):

```toml
[tool.pytest.ini_options]
testpaths = ["tests/zndraw", "tests/zndraw_auth", "tests/zndraw_joblib"]
asyncio_mode = "auto"
asyncio_default_fixture_loop_scope = "function"
markers = ["integration: marks tests as integration tests (require running server)"]
filterwarnings = [
    'error:[\s\S]*You probably want to use .session\.exec:DeprecationWarning',
]
```

**Why single quotes:** avoids double-backslashing the regex in TOML.
**Why `[\s\S]*`:** the SQLModel deprecation message begins with `\n    🚨 You probably...`. Python's `warnings.filterwarnings` compiles the regex without DOTALL, so `.` doesn't cross newlines. `[\s\S]*` matches across newlines without needing literal emoji or whitespace in the pattern.
**Why message-regex (not module):** `typing_extensions.@deprecated` uses `stacklevel=2`, so the warning's `module` attribute points to the calling code (e.g. `zndraw_joblib.registry`), not to `sqlmodel.*`. Filtering by module would silently match nothing.

- [ ] **Step 3: Verify the filter works on a single small test file**

Run:
```bash
uv run --active pytest tests/zndraw_joblib/test_registry.py -q
```

Expected: `8 passed`, no warnings (confirming no leftover src/ call sites are exercised by this file).

- [ ] **Step 4: Run the full test suite to catch any missed call site**

Run:
```bash
uv run --active pytest tests/ 2>&1 | tail -30
```

Expected: all tests pass, **0 warnings**. If a test fails with `DeprecationWarning` escalated to error, the failure traceback will point at the file and line of the missed `session.execute()`. Fix it in the already-committed file-level commit (amend via `git commit --fixup=<sha>` + interactive rebase, OR add a follow-up `fix:` commit).

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml
git commit -m "$(cat <<'EOF'
chore(ci): error on sqlmodel session.execute deprecation warning (#890)

Escalate the SQLModel session.execute() DeprecationWarning to an error
in pytest so any regression into the deprecated API fails CI loudly.

Scoped by message regex (not module) because typing_extensions.@deprecated
uses stacklevel=2, which attributes the warning to the calling code.
The [\\s\\S]* prefix handles the leading newline/emoji of the SQLModel
message since Python's warnings filter compiles regex without DOTALL.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 20: Final verification, push, and PR

**Files:**
- None modified.

- [ ] **Step 1: Full test suite pass with zero warnings**

Run:
```bash
uv run --active pytest tests/ 2>&1 | tail -5
```

Expected: all tests pass, warning count shows `0 warnings` or no warning summary at all.

- [ ] **Step 2: Confirm no leftover `session.execute()` in `src/`**

Run:
```bash
grep -rn "session\.execute(" src/ || echo "CLEAN"
```

Expected: `CLEAN` (no matches).

- [ ] **Step 3: Confirm no leftover `.scalars()` cruft in `src/`**

Run:
```bash
grep -rn "\.scalars(" src/ || echo "CLEAN"
```

Expected: `CLEAN`.

- [ ] **Step 4: Review the commit log**

Run:
```bash
git log --oneline main..HEAD
```

Expected: 19 commits (1 docs spec, 17 per-file migrations, 1 filterwarnings gate).

- [ ] **Step 5: Push the branch**

```bash
git push -u origin worktree-fix+sqlmodel-session-exec-migration
```

- [ ] **Step 6: Open the PR**

```bash
gh pr create --title "fix: migrate to SQLModel session.exec() and gate regressions (#890)" --body "$(cat <<'EOF'
## Summary

- Migrate all 77 `session.execute()` call sites across 17 files in `src/` to SQLModel's `session.exec()` API — the SOTA recommendation per sqlmodel.tiangolo.com.
- Eliminates 15,000+ `DeprecationWarning` instances per test run that were drowning out real warnings in CI.
- Adds a `pyproject.toml` `filterwarnings` gate that escalates any future `session.execute()` reintroduction (on a SQLModel session) to a hard CI failure.

Closes #890.

## Approach

- Per-file atomic commits so bisection is trivial if a regression slips through.
- No behavioral changes. `session.exec()` returns the same data via `.one()` / `.one_or_none()` / `.all()` instead of `.scalar_one()` / `.scalar_one_or_none()` / `.scalars().all()`.
- The one bulk `update(Task).where(...).values(...)` call in `router.py` (optimistic-locking task-claim path) uses SQLModel's `UpdateBase` overload on `exec()`, which returns a fully-typed `CursorResult` with `.rowcount` — no type-ignore, no mixing with raw SQLAlchemy execute.
- Test files are untouched. Empirical verification showed test fixtures create pure `sqlalchemy.ext.asyncio.AsyncSession` instances, which never trigger SQLModel's deprecated wrapper. The 15k warnings all originate from `src/` code exercised via integration tests that go through the production `SQLModelAsyncSession` factory.

Spec: `docs/superpowers/specs/2026-04-08-sqlmodel-session-exec-migration-design.md`
Plan: `docs/superpowers/plans/2026-04-08-sqlmodel-session-exec-migration.md`

## Test plan

- [x] `uv run pytest tests/` passes with 0 warnings.
- [x] `grep -rn "session\\.execute(" src/` returns empty.
- [x] `grep -rn "\\.scalars(" src/` returns empty.
- [x] The `filterwarnings` gate is in effect (a deliberate `session.execute()` in any `src/` file now fails the suite with a clear traceback).

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 7: Return PR URL**

Capture the URL printed by `gh pr create` and report it to the user.

---

## Self-review

- **Spec coverage:** Every file in the spec's 17-file inventory has a dedicated task (Tasks 2-18). The `pyproject.toml` filter is Task 19. Final verification + PR is Task 20. The baseline verification (Task 1) covers the spec's "Verification" section (run tests, establish green baseline). All rewrite rules from the spec (Rules 1-7) are inlined in the "Rewrite rules" reference at the top of the plan and invoked by name in each task's "Rules invoked" column.
- **Placeholder scan:** The `<statement>` placeholders in Tasks 6, 9, 10, 12, 13, 15, 17, 18 are intentional and marked with an instruction to read the file and preserve the exact text. This is not a "TODO" — the mechanical transformation (`execute` → `exec`, `scalar_one_or_none` → `one_or_none`, etc.) is fully specified; only the enclosed statement (which is a 2-5 line SELECT expression) is to be copied verbatim from the existing code. An engineer can perform this without judgment calls.
- **Type consistency:** `session.exec()`, `one_or_none()`, `one()`, `all()` are spelled identically in every task. `.rowcount` is only used in Task 5 Step 5 (Rule 6, bulk update) — consistent with SQLModel's `UpdateBase` overload signature from the spec.
- **Commit atomicity:** Every task produces exactly one commit except Task 1 (no commit — read-only baseline) and Task 20 (no commit — push + PR only). Total commits: 19 (1 spec docs already committed + 17 file migrations + 1 filterwarnings gate).

## Execution handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-08-sqlmodel-session-exec-migration.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session using executing-plans, batch execution with checkpoints.

**Which approach?**
