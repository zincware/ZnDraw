# Global Extension Registration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Enable global (`@global`) and per-room extension registration through `ZnDraw.register_job()`, with deprecation shims for v0.6.0's `register_extension()`.

**Architecture:** Client-side only. Add `room="@global"` support and a deprecated `public` param to `register_job()`, add a deprecated `register_extension()` shim, export `GLOBAL_ROOM` constant. Server enforcement already exists.

**Tech Stack:** Python, typing_extensions, warnings, pytest

**Spec:** `docs/superpowers/specs/2026-03-20-global-extension-registration-design.md`

---

### Task 1: Export `GLOBAL_ROOM` constant

**Files:**
- Modify: `src/zndraw/__init__.py`

- [ ] **Step 1: Add constant and export**

In `src/zndraw/__init__.py`, add `GLOBAL_ROOM = "@global"` and include it in `__all__`:

```python
GLOBAL_ROOM = "@global"
```

Add `"GLOBAL_ROOM"` to the `__all__` list.

- [ ] **Step 2: Verify import works**

Run: `uv run python -c "from zndraw import GLOBAL_ROOM; print(GLOBAL_ROOM)"`
Expected: `@global`

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/__init__.py
git commit -m "feat: export GLOBAL_ROOM constant"
```

---

### Task 2: Update `register_job` with `room="@global"` and deprecated `public` param

**Files:**
- Modify: `src/zndraw/client/core.py:592-603`

The current `register_job` always goes through `_resolve_room()`, which defaults to `self.room`. We need to:
1. Accept `room="@global"` and pass it through directly (skip `_resolve_room`)
2. Add deprecated `public` param that maps `True` → `room="@global"`
3. Emit runtime `DeprecationWarning` when `public=True` is used

- [ ] **Step 1: Write the failing test for `room="@global"`**

Create `tests/worker/test_register_job_api.py`:

```python
"""Tests for ZnDraw.register_job() high-level API.

Covers: room='@global', deprecated public param, deprecated register_extension,
and validation of conflicting arguments.
"""

import warnings

import pytest

from zndraw import GLOBAL_ROOM, ZnDraw


def test_register_job_global(server, Echo, get_job_list):
    """register_job(cls, room='@global') registers a global job."""
    worker = ZnDraw(url=server)
    try:
        worker.register_job(Echo, room=GLOBAL_ROOM)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_job_global -v`
Expected: FAIL — `register_job` passes `"@global"` to `_resolve_room` which treats it as a room name

- [ ] **Step 3: Implement the updated `register_job`**

In `src/zndraw/client/core.py`, replace the current `register_job` method (lines 592-603) with:

```python
    def register_job(
        self,
        cls: type,
        *,
        room: Literal["@global"] | str | None = None,
        public: Annotated[
            bool | None,
            typing_extensions.deprecated(
                "Use room='@global' instead of public=True"
            ),
        ] = None,
    ) -> None:
        """Register an extension as a job. Connects the socket if needed.

        Parameters
        ----------
        cls
            Extension subclass to register.
        room
            Room scope. Use ``"@global"`` for global registration (admin-only).
            Defaults to ``self.room``.
        public
            .. deprecated::
                Use ``room='@global'`` instead.
        """
        if public and room is not None:
            raise ValueError("Cannot specify both 'room' and 'public'")
        if public:
            warnings.warn(
                "public=True is deprecated, use room='@global' instead",
                DeprecationWarning,
                stacklevel=2,
            )
            room = GLOBAL_ROOM
        elif room != GLOBAL_ROOM:
            room = self._resolve_room(room)

        self._ensure_socket_connected()
        self.jobs.register(cls, room=room)
```

Add necessary imports at the top of the file. `warnings` is already imported. Add:

```python
from typing import Annotated, Literal
```

(`Literal` needs to be added; `Annotated` may need to be added — check existing imports.)

Also import `GLOBAL_ROOM`:

```python
from zndraw import GLOBAL_ROOM
```

**Note:** Since `core.py` is inside the `zndraw` package, use a relative import or import the constant value directly. The simplest approach: define `GLOBAL_ROOM = "@global"` at module level in `core.py` and also in `__init__.py`. Alternatively, use a late import or just inline `"@global"`. The cleanest option: since `__init__.py` imports from `client`, use the string literal `"@global"` in `core.py` and reference the `GLOBAL_ROOM` constant only in `__init__.py` and tests. The implementation code can use the string directly since it's a simple sentinel.

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_job_global -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/client/core.py tests/worker/test_register_job_api.py
git commit -m "feat: register_job supports room='@global' with deprecated public param"
```

---

### Task 3: Add deprecated `register_extension` method

**Files:**
- Modify: `src/zndraw/client/core.py` (add method after `register_job`)

- [ ] **Step 1: Write the failing test**

Append to `tests/worker/test_register_job_api.py`:

```python
def test_register_extension_public(server, Echo, get_job_list):
    """Deprecated register_extension(cls, public=True) registers globally."""
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_extension(Echo, public=True)
            assert any(issubclass(x.category, DeprecationWarning) for x in w)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_extension_public -v`
Expected: FAIL — `AttributeError: 'ZnDraw' object has no attribute 'register_extension'`

- [ ] **Step 3: Implement `register_extension`**

In `src/zndraw/client/core.py`, add after the `register_job` method:

```python
    @typing_extensions.deprecated(
        "Use register_job(cls, room='@global') for global, "
        "or register_job(cls) for room-scoped"
    )
    def register_extension(
        self, cls: type, *, public: bool = False, **kwargs: Any
    ) -> None:
        """Register an extension.

        .. deprecated::
            Use :meth:`register_job` instead.
        """
        room = "@global" if public else kwargs.get("room")
        self.register_job(cls, room=room)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_extension_public -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/client/core.py tests/worker/test_register_job_api.py
git commit -m "feat: add deprecated register_extension shim"
```

---

### Task 4: Test deprecated `public=True` param on `register_job`

**Files:**
- Modify: `tests/worker/test_register_job_api.py`

- [ ] **Step 1: Write the test**

Append to `tests/worker/test_register_job_api.py`:

```python
def test_register_job_public_deprecated(server, Echo, get_job_list):
    """register_job(cls, public=True) works but emits DeprecationWarning."""
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_job(Echo, public=True)
            dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(dep_warnings) >= 1
            assert "room='@global'" in str(dep_warnings[0].message)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
```

- [ ] **Step 2: Run test to verify it passes**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_job_public_deprecated -v`
Expected: PASS (implementation already done in Task 2)

- [ ] **Step 3: Commit**

```bash
git add tests/worker/test_register_job_api.py
git commit -m "test: register_job public=True emits deprecation warning"
```

---

### Task 5: Test `ValueError` when both `room` and `public` are set

**Files:**
- Modify: `tests/worker/test_register_job_api.py`

- [ ] **Step 1: Write the test**

Append to `tests/worker/test_register_job_api.py`:

```python
def test_register_job_room_and_public_raises(server, Echo):
    """Passing both room= and public=True raises ValueError."""
    vis = ZnDraw(url=server)
    try:
        with pytest.raises(ValueError, match="Cannot specify both"):
            vis.register_job(Echo, room=vis.room, public=True)
    finally:
        vis.disconnect()
```

- [ ] **Step 2: Run test to verify it passes**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_register_job_room_and_public_raises -v`
Expected: PASS (implementation already done in Task 2)

- [ ] **Step 3: Commit**

```bash
git add tests/worker/test_register_job_api.py
git commit -m "test: ValueError when room and public both specified"
```

---

### Task 6: Test global extension visible in all rooms

**Files:**
- Modify: `tests/worker/test_register_job_api.py`

This is the critical cross-room visibility test using the high-level `register_job` API.

- [ ] **Step 1: Write the test**

Append to `tests/worker/test_register_job_api.py`:

```python
def test_global_extension_visible_in_all_rooms(server, Echo, get_job_list):
    """A @global extension registered via register_job is visible from any room."""
    registrar = ZnDraw(url=server)
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        # Register globally from registrar's context
        registrar.register_job(Echo, room=GLOBAL_ROOM)

        # Verify visible from room_a
        jobs_a = get_job_list(room_a, room_id=room_a.room)
        global_a = [j for j in jobs_a if j.full_name.startswith("@global")]
        names_a = {j.name for j in global_a}
        assert "Echo" in names_a

        # Verify visible from room_b (different room)
        jobs_b = get_job_list(room_b, room_id=room_b.room)
        global_b = [j for j in jobs_b if j.full_name.startswith("@global")]
        names_b = {j.name for j in global_b}
        assert "Echo" in names_b
    finally:
        registrar.jobs.disconnect()
        registrar.disconnect()
        room_a.disconnect()
        room_b.disconnect()
```

- [ ] **Step 2: Run test to verify it passes**

Run: `uv run pytest tests/worker/test_register_job_api.py::test_global_extension_visible_in_all_rooms -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tests/worker/test_register_job_api.py
git commit -m "test: global extension visible in all rooms"
```

---

### Task 7: Test admin enforcement via `register_job`

**Files:**
- Modify: `tests/worker/test_register_job_api.py`

- [ ] **Step 1: Write the tests**

Append to `tests/worker/test_register_job_api.py`:

```python
def test_guest_cannot_register_global_via_register_job(server_auth, Echo):
    """Guest user gets 403 when using register_job(room='@global')."""
    guest = ZnDraw(url=server_auth)
    try:
        with pytest.raises(PermissionError):
            guest.register_job(Echo, room=GLOBAL_ROOM)
    finally:
        guest.disconnect()


def test_admin_can_register_global_via_register_job(server_auth, Echo, get_job_list):
    """Admin user can register global jobs via register_job."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    try:
        admin.register_job(Echo, room=GLOBAL_ROOM)
        jobs = get_job_list(admin, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        admin.jobs.disconnect()
        admin.disconnect()
```

- [ ] **Step 2: Run tests to verify they pass**

Run: `uv run pytest tests/worker/test_register_job_api.py -k "admin or guest" -v`
Expected: PASS

- [ ] **Step 3: Commit**

```bash
git add tests/worker/test_register_job_api.py
git commit -m "test: admin enforcement for register_job(room='@global')"
```

---

### Task 8: Run full test suite and fix any issues

- [ ] **Step 1: Run all worker tests**

Run: `uv run pytest tests/worker/ -v`
Expected: All pass (existing tests + new tests)

- [ ] **Step 2: Run type checker**

Run: `uv run pyright src/zndraw/client/core.py src/zndraw/__init__.py`
Expected: No new errors

- [ ] **Step 3: Run formatter and import sorter**

Run: `uv run ruff format . && uv run ruff check --select I --fix .`

- [ ] **Step 4: Final commit if any formatting changes**

```bash
git add -u
git commit -m "style: format and fix imports"
```

---

### Task 9: Run complete test suite

- [ ] **Step 1: Run all tests**

Run: `uv run pytest tests/ -v`
Expected: All 499+ tests pass. No regressions.

Note: Tests can run up to 15 minutes — be patient.
