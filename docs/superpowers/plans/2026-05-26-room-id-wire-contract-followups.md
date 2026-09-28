# Room-ID Wire Contract — Review Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the 11 findings from the code review against `4f21d6cc..43b143ad` so the room-ID wire convention refactor lands without contract gaps.

**Architecture:** Each finding becomes a small, test-first task. Critical findings (#1, #2) close contract violations the original refactor was meant to prevent. Important findings (#3–#7) extend coverage and tighten ergonomics. Minor findings (#8, #10) are polish. Frequent commits, one task per commit.

**Tech Stack:** Python 3.12 + FastAPI + Pydantic v2 + Socket.IO (`zndraw-socketio`); SQLModel/SQLAlchemy async; pytest + pytest-asyncio; React/TanStack Query + bun-built frontend.

---

## Background — what the reviewer flagged

| # | Severity | One-line |
|---|----------|---------|
| 1 | Critical | `FrameSelectionUpdate` is emitted to `room:{id}` but is not a `RoomScopedEvent`. |
| 2 | Critical | `JobsInvalidate`/`ProvidersInvalidate`/`ProviderResultReady` direct construction can silently emit `NIL_ROOM_UUID` for what should be a real room. |
| 3 | Important | `RoomRenamed` notification is skipped when the previous owner was a group — members never get the rename signal. |
| 4 | Important | Runtime contract sweep covers ~6 of ~18 event families. |
| 5 | Important | `_room_address_for` duplicated and divergent between `joblib/router.py` and `joblib/sweeper.py`. |
| 6 | Important | Frontend non-`FramesInvalidate` handlers don't compare `room_address`. |
| 7 | Important | Schema audit relies on `__subclasses__()`; only sees imported classes. |
| 8 | Minor | `MockSioServer._skip_sid` parameter doesn't match production helper. |
| 9 | Minor | Joblib `JobsInvalidate` sigil/real branching is repeated and could be a helper. (Folded into Task 3.) |
| 10 | Minor | Six working-tree files contain pure-format drift — should be committed cleanly. |
| 11 | Minor | `room_id: string` field in `FramesInvalidateEvent` is unused. *(Kept for documentation; no task.)* |

---

## File Structure

**Python — source**
- `src/zndraw/socket_events.py` — promote `FrameSelectionUpdate` to `RoomScopedEvent`; add a `model_validator` rejecting NIL `room_id` with a composed `room_address`.
- `src/zndraw/broadcast.py` — accept an iterable of user-ids for fan-out (replaces singular `also_notify_user`).
- `src/zndraw/routes/utility.py` — route the `frame-selection` PUT through `broadcast_to_room`.
- `src/zndraw/routes/rooms.py` — capture `previous_owner_group_id` and expand `RoomRenamed` fan-out to group members.
- `src/zndraw_joblib/room_lookup.py` *(new)* — single home for `fetch_room` / `room_address_for`.
- `src/zndraw_joblib/events.py` — new `build_room_scoped_emission(...)` helper that returns `(event, channel)` for any `RoomScopedEvent` subclass, encapsulating the sigil-vs-real branch.
- `src/zndraw_joblib/router.py` — call the new helper in `register_job`, `register_provider`, `delete_provider`, `upload_provider_result`; import shared `room_lookup`.
- `src/zndraw_joblib/sweeper.py` — import shared `room_lookup`; delete local `_room_address_for`.

**Python — tests**
- `tests/zndraw/test_frame_selection_address.py` *(new)* — red regression for finding #1.
- `tests/zndraw/test_room_scoped_validator.py` *(new)* — red regression for finding #2 (validator rejects NIL+composed).
- `tests/zndraw/test_room_renamed_group_fanout.py` *(new)* — red regression for finding #3.
- `tests/zndraw/test_broadcast_contract.py` — extend matrix (finding #4).
- `tests/zndraw/test_event_schema_audit.py` — add docstring explaining the `__subclasses__` limitation (finding #7).
- `tests/zndraw/helpers.py` — rename `_skip_sid` → `skip_sid` (finding #8).
- `tests/zndraw_joblib/test_room_lookup.py` *(new)* — covers the consolidated module (finding #5).

**Frontend**
- `frontend/src/hooks/socketHandlers/utils.ts` — gate `createInvalidateHandler` on `data.room_address` (finding #6).
- `frontend/src/hooks/socketHandlers/geometryHandlers.ts` — explicit `room_address` gate on `onGeometriesInvalidate` (finding #6).

---

## Task 1: Commit working-tree drift

**Files:**
- Modify: `frontend/src/hooks/socketHandlers/frameHandlers.ts` (already drift)
- Modify: `src/zndraw/routes/edit_lock.py` (already drift)
- Modify: `src/zndraw/routes/trajectory.py` (already drift)
- Modify: `src/zndraw/routes/utility.py` (already drift)
- Modify: `src/zndraw/socket_events.py` (already drift — `Room` annotation requote)
- Modify: `tests/zndraw_joblib/test_events.py` (already drift — wrapping)

The drift is pure formatting (import-order, line-wrap, `from __future__ import annotations`-safe quote removal) introduced by ruff/biome runs. Verified safe in the review.

- [ ] **Step 1: Confirm drift is whitespace-only**

```bash
git diff HEAD --stat
git diff HEAD | grep -E "^[+-]" | grep -v "^[+-]\{3\}" | grep -vE "^[+-]\s*$" | head -40
```

Expected: only import-order swaps, line-wraps, and the `"Room"` → `Room` requote on `socket_events.py:29`. No semantic identifier or operator changes.

- [ ] **Step 2: Stage the six files and commit**

```bash
git add frontend/src/hooks/socketHandlers/frameHandlers.ts \
        src/zndraw/routes/edit_lock.py \
        src/zndraw/routes/trajectory.py \
        src/zndraw/routes/utility.py \
        src/zndraw/socket_events.py \
        tests/zndraw_joblib/test_events.py
git commit -m "$(cat <<'EOF'
chore(format): absorb pending whitespace drift

Pure import-order, line-wrap, and runtime-annotation quote cleanup.
No semantic changes; surfaces in the working tree from the previous
ruff/biome runs.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: Confirm clean tree**

```bash
git status --short
```

Expected: empty.

---

## Task 2: FrameSelectionUpdate joins RoomScopedEvent (Critical #1)

**Files:**
- Create: `tests/zndraw/test_frame_selection_address.py`
- Modify: `src/zndraw/socket_events.py:147-150`
- Modify: `src/zndraw/routes/utility.py:41,136-139`
- Modify: `frontend/src/hooks/socketHandlers/frameHandlers.ts:19-21`

- [ ] **Step 1: Write the failing regression test**

Create `tests/zndraw/test_frame_selection_address.py`:

```python
"""Regression: FrameSelectionUpdate carries room_address and routes via broadcast_to_room (review #1)."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession


@pytest.mark.asyncio
async def test_frame_selection_update_carries_room_address(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.public_address}/frame-selection",
        json={"indices": [0, 1, 2]},
        headers=auth_header(token),
    )
    assert response.status_code == 200, response.text

    emits = [e for e in mock_sio.emitted if e["event"] == "frame_selection_update"]
    assert len(emits) == 1, f"expected 1 frame_selection_update, got {len(emits)}"
    captured = emits[0]
    data = captured["data"]
    assert "room_id" in data, (
        f"FrameSelectionUpdate missing room_id; got keys: {list(data)}"
    )
    assert "room_address" in data, (
        f"FrameSelectionUpdate missing room_address; got keys: {list(data)}"
    )
    assert str(data["room_id"]) == room.id
    assert data["room_address"] == room.public_address
    assert captured["room"] == f"room:{room.id}"
```

- [ ] **Step 2: Run the test to confirm RED**

```bash
uv run pytest tests/zndraw/test_frame_selection_address.py -v
```

Expected: FAIL (`room_id`/`room_address` missing from payload).

- [ ] **Step 3: Promote `FrameSelectionUpdate` to `RoomScopedEvent`**

Edit `src/zndraw/socket_events.py:147-150`:

```python
class FrameSelectionUpdate(RoomScopedEvent):
    """Broadcast when frame selection changes."""

    indices: list[int]
```

- [ ] **Step 4: Route the emit through `broadcast_to_room`**

Edit `src/zndraw/routes/utility.py`:

Replace line 12 import:
```python
from zndraw.broadcast import broadcast_to_room
```

Drop `room_channel` from the same import (it's no longer used in this file).

Replace lines 136-139:
```python
    await broadcast_to_room(
        sio,
        FrameSelectionUpdate.for_room(room, indices=body.indices),
        room,
    )
```

- [ ] **Step 5: Update frontend event interface**

Edit `frontend/src/hooks/socketHandlers/frameHandlers.ts:19-21`:

```typescript
export interface FrameSelectionUpdateEvent {
	room_id: string;
	room_address: string;
	indices: number[] | null;
}
```

- [ ] **Step 6: Run tests + schema audit + contract sweep**

```bash
uv run pytest tests/zndraw/test_frame_selection_address.py \
              tests/zndraw/test_event_schema_audit.py \
              tests/zndraw/test_broadcast_contract.py -v
```

Expected: all PASS. The static audit picks up `FrameSelectionUpdate` automatically via `__subclasses__()`.

- [ ] **Step 7: Frontend type check**

```bash
cd frontend && bun run tsc --noEmit
```

Expected: no errors.

- [ ] **Step 8: Commit**

```bash
git add tests/zndraw/test_frame_selection_address.py \
        src/zndraw/socket_events.py \
        src/zndraw/routes/utility.py \
        frontend/src/hooks/socketHandlers/frameHandlers.ts
git commit -m "$(cat <<'EOF'
fix(events): FrameSelectionUpdate joins RoomScopedEvent

Closes review finding #1: the event was emitted to `room:{id}` but
carried no `room_id`/`room_address`, breaking the contract the refactor
was meant to enforce. Route through `broadcast_to_room` and update the
frontend payload type.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: NIL_UUID consistency validator + joblib emit helper (Critical #2 + Minor #9)

**Files:**
- Create: `tests/zndraw/test_room_scoped_validator.py`
- Modify: `src/zndraw/socket_events.py:22-34` (add `model_validator`)
- Modify: `src/zndraw_joblib/events.py` (new `build_room_scoped_emission` helper)
- Modify: `src/zndraw_joblib/router.py:421-431,1205-1217,1432-1444,1505-1520` (use the helper)

- [ ] **Step 1: Write the failing validator test**

Create `tests/zndraw/test_room_scoped_validator.py`:

```python
"""Regression: RoomScopedEvent rejects NIL room_id paired with a composed address (review #2)."""

from uuid import UUID

import pytest
from pydantic import ValidationError

from zndraw.socket_events import FramesInvalidate


def test_nil_room_id_with_composed_address_rejected() -> None:
    """Direct construction with NIL room_id but a composed address is a footgun."""
    composed = "11111111-1111-1111-1111-111111111111/demo"
    with pytest.raises(ValidationError, match="room_id"):
        FramesInvalidate(
            room_id=UUID(int=0),
            room_address=composed,
            action="add",
        )


def test_nil_room_id_with_sigil_address_allowed() -> None:
    """Sigil addresses (``@global``/``@internal``) intentionally pair with NIL room_id."""
    event = FramesInvalidate(
        room_id=UUID(int=0),
        room_address="@global",
        action="clear",
    )
    assert event.room_id == UUID(int=0)
    assert event.room_address == "@global"


def test_real_room_id_with_composed_address_allowed() -> None:
    real = UUID("22222222-2222-2222-2222-222222222222")
    composed = f"{real}/demo"
    event = FramesInvalidate(
        room_id=real,
        room_address=composed,
        action="modify",
        indices=[0],
    )
    assert event.room_id == real
```

- [ ] **Step 2: Run the test to confirm RED**

```bash
uv run pytest tests/zndraw/test_room_scoped_validator.py -v
```

Expected: FAIL — first test does not raise because the validator is absent.

- [ ] **Step 3: Add the validator to `RoomScopedEvent`**

Edit `src/zndraw/socket_events.py:22-34`. Replace the existing class body with:

```python
class RoomScopedEvent(BaseModel):
    """Base class for room-scoped broadcast events."""

    room_id: UUID
    room_address: str

    @classmethod
    def for_room(cls, room: Room, /, **kwargs: Any) -> Self:
        return cls(
            room_id=UUID(room.id),
            room_address=room.public_address,
            **kwargs,
        )

    @model_validator(mode="after")
    def _validate_room_id_address_consistency(self) -> Self:
        if self.room_id == _NIL_ROOM_UUID and _is_composed_address(self.room_address):
            raise ValueError(
                f"room_id must not be NIL when room_address looks composed: "
                f"room_id={self.room_id}, room_address={self.room_address!r}"
            )
        return self
```

Add the helpers and import near the top of the file (after the existing imports, before `class RoomScopedEvent`):

```python
from pydantic import BaseModel, model_validator

_NIL_ROOM_UUID = UUID(int=0)


def _is_composed_address(address: str) -> bool:
    """Composed addresses look like ``<uuid>/<name>``; sigils start with ``@``."""
    return "/" in address and not address.startswith("@")
```

Note: the existing `from pydantic import BaseModel` line (currently line 14) must be replaced with the line above.

- [ ] **Step 4: Run the validator test — expect GREEN**

```bash
uv run pytest tests/zndraw/test_room_scoped_validator.py -v
```

Expected: all three pass.

- [ ] **Step 5: Add the joblib emit helper**

Edit `src/zndraw_joblib/events.py`. Append before the existing `def emit(...)`:

```python
async def build_room_scoped_emission(
    session: AsyncSession,
    event_cls: type[RoomScopedEvent],
    room_id: str,
    **fields: Any,
) -> Emission:
    """Construct an Emission for a joblib room-scoped event.

    Handles the sigil/real-room split centrally:
      * for real rooms, defers to ``event_cls.for_room(room, **fields)``;
      * for sigils or unknown ids, builds the event by hand with
        ``room_id=NIL_ROOM_UUID`` and ``room_address=room_id`` (passes the
        ``RoomScopedEvent`` validator because sigil addresses do not look
        composed).

    The channel is always ``f"room:{room_id}"``.
    """
    from zndraw_joblib.room_lookup import fetch_room  # local import — see Task 4

    room = await fetch_room(session, room_id)
    if room is not None:
        return Emission(event_cls.for_room(room, **fields), f"room:{room.id}")
    return Emission(
        event_cls(room_id=NIL_ROOM_UUID, room_address=room_id, **fields),
        f"room:{room_id}",
    )
```

Also add the `AsyncSession` import at the top:

```python
from sqlmodel.ext.asyncio.session import AsyncSession
```

(Task 4 creates `room_lookup`; Task 3 takes a local import to avoid a chicken-and-egg, then Task 4 makes it valid.)

- [ ] **Step 6: Replace direct constructions in `router.py`**

There are three callsites to convert. For each, replace the inline `Emission(JobsInvalidate(...), f"room:{room_id}")` (and friends) with the helper.

In `register_job` (`router.py:421-431`), replace the block:

```python
await session.commit()
emission = await build_room_scoped_emission(session, JobsInvalidate, room_id)
await emit(tsio, {emission})
```

In `register_provider` (`router.py:1205-1217`):

```python
    await session.commit()
    await session.refresh(provider)
    emission = await build_room_scoped_emission(
        session, ProvidersInvalidate, provider.room_id
    )
    await emit(tsio, {emission})
```

In `delete_provider` (`router.py:1432-1444`):

```python
room_id = provider.room_id
await session.delete(provider)
await session.commit()
emission = await build_room_scoped_emission(session, ProvidersInvalidate, room_id)
await emit(tsio, {emission})
```

In `upload_provider_result` (`router.py:1505-1520`):

```python
    async with session_maker() as session_addr:
        emission = await build_room_scoped_emission(
            session_addr,
            ProviderResultReady,
            provider.room_id,
            provider_name=provider.full_name,
            request_hash=x_request_hash,
        )
    await emit(tsio, {emission})
```

Add `build_room_scoped_emission` to the existing `from zndraw_joblib.events import (...)` block at the top of `router.py`. Drop now-unused `event_room_uuid` from this same import if no other site in router.py references it (verify with `grep -n event_room_uuid src/zndraw_joblib/router.py`; sweeper.py and the test layer still use it).

- [ ] **Step 7: Run joblib tests + contract sweep**

```bash
uv run pytest tests/zndraw_joblib tests/zndraw/test_broadcast_contract.py \
              tests/zndraw/test_room_scoped_validator.py -v
```

Expected: all PASS.

- [ ] **Step 8: Commit**

```bash
git add tests/zndraw/test_room_scoped_validator.py \
        src/zndraw/socket_events.py \
        src/zndraw_joblib/events.py \
        src/zndraw_joblib/router.py
git commit -m "$(cat <<'EOF'
fix(events): reject NIL room_id paired with composed room_address

Closes review findings #2 and #9. Adds a model_validator on
RoomScopedEvent that rejects ``room_id == NIL`` when ``room_address``
looks composed (``<uuid>/<name>``), and centralizes joblib
``JobsInvalidate``/``ProvidersInvalidate``/``ProviderResultReady``
construction through ``build_room_scoped_emission`` so the sigil/real
branch lives in one place.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Shared room_lookup module (Important #5)

**Files:**
- Create: `src/zndraw_joblib/room_lookup.py`
- Create: `tests/zndraw_joblib/test_room_lookup.py`
- Modify: `src/zndraw_joblib/router.py:107-157` (delete local `_fetch_room` + `_room_address_for`, import from new module)
- Modify: `src/zndraw_joblib/sweeper.py:39-51` (delete local `_room_address_for`, import from new module)

- [ ] **Step 1: Write the consolidated module**

Read the canonical `_fetch_room` in `src/zndraw_joblib/router.py:107-148` first to capture its exact composed-address fallback. Then create `src/zndraw_joblib/room_lookup.py`:

```python
"""Shared room-lookup helpers for joblib emit machinery.

Both ``router.py`` and ``sweeper.py`` need to map a joblib ``room_id``
(real UUID, composed ``<uuid>/<name>``, ``@global``, ``@internal``, or
unknown garbage from the wire) to either the persisted ``Room`` row or
the address that should be surfaced on the wire when no row exists.

Joblib-only test environments may not have the ``Room`` table — the
``OperationalError`` / ``ProgrammingError`` catch keeps these helpers
usable in those harnesses.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from uuid import UUID

from sqlalchemy.exc import OperationalError, ProgrammingError

if TYPE_CHECKING:
    from sqlmodel.ext.asyncio.session import AsyncSession

    from zndraw.models import Room


async def fetch_room(session: AsyncSession, room_id: str) -> Room | None:
    """Return the Room for ``room_id``, or None.

    Accepts a UUID string, a composed ``<uuid>/<name>`` address, or a
    sigil. Returns None for sigils, unknown ids, and missing rooms.
    """
    if room_id in ("@global", "@internal"):
        return None
    from zndraw.models import Room  # local import — joblib must work without zndraw

    try:
        try:
            UUID(room_id)
        except ValueError:
            # Composed address fallback: split into owner / name.
            if "/" not in room_id:
                return None
            owner_part, _, name_part = room_id.partition("/")
            try:
                owner_uuid = UUID(owner_part)
            except ValueError:
                return None
            from sqlmodel import col, or_, select as sql_select

            result = await session.exec(
                sql_select(Room).where(
                    or_(
                        Room.owner_user_id == owner_uuid,
                        Room.owner_group_id == owner_uuid,
                    ),
                    col(Room.room_name) == name_part,
                )
            )
            return result.one_or_none()
        return await session.get(Room, room_id)
    except (OperationalError, ProgrammingError):
        return None


async def room_address_for(session: AsyncSession, room_id: str) -> str:
    """Return ``room.public_address`` for a known room, else echo ``room_id``."""
    if room_id in ("@global", "@internal"):
        return room_id
    room = await fetch_room(session, room_id)
    if room is None:
        return room_id
    return room.public_address
```

- [ ] **Step 2: Write tests for the helper**

Create `tests/zndraw_joblib/test_room_lookup.py`:

```python
"""Smoke tests for the consolidated room_lookup helpers."""

import pytest
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_joblib.room_lookup import fetch_room, room_address_for


@pytest.mark.asyncio
async def test_sigils_return_self(session: AsyncSession) -> None:
    assert await fetch_room(session, "@global") is None
    assert await fetch_room(session, "@internal") is None
    assert await room_address_for(session, "@global") == "@global"
    assert await room_address_for(session, "@internal") == "@internal"


@pytest.mark.asyncio
async def test_unknown_uuid_returns_self(session: AsyncSession) -> None:
    assert await fetch_room(session, "99999999-9999-9999-9999-999999999999") is None
    assert (
        await room_address_for(session, "99999999-9999-9999-9999-999999999999")
        == "99999999-9999-9999-9999-999999999999"
    )


@pytest.mark.asyncio
async def test_garbage_string_returns_self(session: AsyncSession) -> None:
    assert await fetch_room(session, "not-a-uuid") is None
    assert await room_address_for(session, "not-a-uuid") == "not-a-uuid"
```

Note: the existing `tests/zndraw_joblib/conftest.py` provides `session`. If those tests run in a joblib-only DB (no Room table), the `OperationalError` catch keeps them green.

- [ ] **Step 3: Run new tests — expect GREEN**

```bash
uv run pytest tests/zndraw_joblib/test_room_lookup.py -v
```

- [ ] **Step 4: Delete the duplicates and rewire imports**

In `src/zndraw_joblib/router.py`:
- Delete the local `_fetch_room` definition (`router.py:107-148`) and the local `_room_address_for` (`router.py:150-157`).
- Add to the imports at the top:

```python
from zndraw_joblib.room_lookup import fetch_room, room_address_for
```

- Replace every call to `_fetch_room(...)` with `fetch_room(...)` (one site at `router.py:421` for `register_job`; verify with `grep -n _fetch_room src/zndraw_joblib/router.py`).
- Replace every call to `_room_address_for(...)` with `room_address_for(...)` (sites around `router.py:283`, `1205`, `1432`, `1506`).

In `src/zndraw_joblib/sweeper.py`:
- Delete the local `_room_address_for` (`sweeper.py:39-51`).
- Replace its import block:

```python
from zndraw_joblib.room_lookup import room_address_for
```

- Replace every call to `_room_address_for(...)` with `room_address_for(...)` (sites at `sweeper.py:96`, `145`, `167`, `191`, `304`).

Also update `src/zndraw_joblib/events.py` (added in Task 3 Step 5) to use the new public name:

```python
    from zndraw_joblib.room_lookup import fetch_room
```

- [ ] **Step 5: Run full joblib + zndraw tests**

```bash
uv run pytest tests/zndraw_joblib tests/zndraw -v
```

Expected: all PASS.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_joblib/room_lookup.py \
        tests/zndraw_joblib/test_room_lookup.py \
        src/zndraw_joblib/router.py \
        src/zndraw_joblib/sweeper.py \
        src/zndraw_joblib/events.py
git commit -m "$(cat <<'EOF'
refactor(joblib): consolidate room_lookup into shared module

Closes review finding #5. ``router.py`` and ``sweeper.py`` carried
divergent copies of ``_room_address_for``; sweeper's lacked the
composed-address fallback. Moves both helpers into
``zndraw_joblib.room_lookup`` with the canonical (router-side)
behavior and rewires every callsite.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: RoomRenamed group fanout (Important #3)

**Files:**
- Create: `tests/zndraw/test_room_renamed_group_fanout.py`
- Modify: `src/zndraw/broadcast.py:15-32` (broaden `also_notify_user` → `also_notify_user_ids`)
- Modify: `src/zndraw/routes/rooms.py:645,713-719` (capture previous group, fan out)

The current API has `also_notify_user: UUID | str | None`. Replace it with an iterable parameter so the same helper handles user-owned and group-owned previous-owners. Only one existing call passes the parameter (rooms.py:718), so the migration is contained.

- [ ] **Step 1: Write the failing fan-out test**

Create `tests/zndraw/test_room_renamed_group_fanout.py`:

```python
"""Regression: RoomRenamed fans out to every previous group member (review #3)."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.access import GroupRole, Visibility
from zndraw.models import Group, GroupMembership, Room


@pytest.mark.asyncio
async def test_room_renamed_fans_out_to_previous_group_members(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    alice, alice_token = await create_test_user_in_db(
        session, email="alice@local.test", is_superuser=True
    )
    bob, _ = await create_test_user_in_db(session, email="bob@local.test")
    carol, _ = await create_test_user_in_db(session, email="carol@local.test")

    # Group G owns room "foo" via alice (admin), bob, and carol.
    group = Group(name="G", created_by_id=alice.id)
    session.add(group)
    await session.flush()
    session.add(
        GroupMembership(user_id=alice.id, group_id=group.id, role=GroupRole.ADMIN)
    )
    session.add(
        GroupMembership(user_id=bob.id, group_id=group.id, role=GroupRole.MEMBER)
    )
    session.add(
        GroupMembership(user_id=carol.id, group_id=group.id, role=GroupRole.MEMBER)
    )
    await session.commit()

    room = Room(
        room_name="foo",
        owner_group_id=group.id,
        visibility=Visibility.GROUP,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    # Transfer from group G back to alice (a user).
    response = await client.patch(
        f"/v1/rooms/{room.public_address}",
        json={"new_owner_id": str(alice.id), "visibility": Visibility.PRIVATE.value},
        headers=auth_header(alice_token),
    )
    assert response.status_code == 200, response.text

    renamed = [e for e in mock_sio.emitted if e["event"] == "room_renamed"]
    rooms_seen = {e["room"] for e in renamed}
    for member in (alice, bob, carol):
        assert f"user:{member.id}" in rooms_seen, (
            f"RoomRenamed missing for previous group member {member.email}; "
            f"saw {rooms_seen}"
        )
```

- [ ] **Step 2: Run — expect RED**

```bash
uv run pytest tests/zndraw/test_room_renamed_group_fanout.py -v
```

Expected: FAIL — bob and carol miss the rename event.

- [ ] **Step 3: Broaden `broadcast_to_room` API**

Edit `src/zndraw/broadcast.py:15-32`:

```python
from collections.abc import Iterable
from uuid import UUID

from zndraw_socketio import AsyncServerWrapper

from zndraw.models import Room
from zndraw.socket_events import RoomScopedEvent


def room_channel(room_id: str | UUID) -> str:
    return f"room:{room_id}"


async def broadcast_to_room(
    sio: AsyncServerWrapper,
    event: RoomScopedEvent,
    room: Room,
    *,
    also_notify_user_ids: Iterable[UUID | str] | None = None,
    skip_sid: str | None = None,
) -> None:
    """Emit ``event`` on the room channel; optionally fan out to user channels."""
    assert UUID(room.id) == event.room_id, (
        f"event.room_id {event.room_id} does not match room.id {room.id}"
    )
    if skip_sid is not None:
        await sio.emit(event, room=room_channel(room.id), skip_sid=skip_sid)
    else:
        await sio.emit(event, room=room_channel(room.id))
    if also_notify_user_ids:
        for uid in also_notify_user_ids:
            await sio.emit(event, room=f"user:{uid}")
```

- [ ] **Step 4: Wire up `routes/rooms.py`**

Edit `src/zndraw/routes/rooms.py:645`. Add a line right after the existing `previous_owner_user_id` capture:

```python
    previous_owner_user_id: UUID | None = room.owner_user_id
    previous_owner_group_id: UUID | None = room.owner_group_id
```

Edit `src/zndraw/routes/rooms.py:713-719`. Replace the existing block with:

```python
    if updates.new_owner_id is not None:
        from zndraw.models import GroupMembership

        prev_user_ids: list[UUID] = []
        if previous_owner_group_id is not None:
            result = await session.exec(
                select(GroupMembership.user_id).where(
                    GroupMembership.group_id == previous_owner_group_id
                )
            )
            prev_user_ids = list(result.all())
        elif previous_owner_user_id is not None:
            prev_user_ids = [previous_owner_user_id]

        await broadcast_to_room(
            sio,
            RoomRenamed.for_room(room, old_address=old_address),
            room,
            also_notify_user_ids=prev_user_ids,
        )
```

- [ ] **Step 5: Run the new test + the existing user-channel regression**

```bash
uv run pytest tests/zndraw/test_room_renamed_group_fanout.py \
              tests/zndraw/test_room_renamed_previous_owner.py -v
```

Expected: both PASS. (The existing user-channel test still uses single-user previous-owner, which the new API handles via the `[previous_owner_user_id]` branch.)

- [ ] **Step 6: Run the broader suite**

```bash
uv run pytest tests/zndraw -v
```

Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add tests/zndraw/test_room_renamed_group_fanout.py \
        src/zndraw/broadcast.py \
        src/zndraw/routes/rooms.py
git commit -m "$(cat <<'EOF'
fix(rooms): RoomRenamed fans out to previous group members

Closes review finding #3. The previous owner could be a group, in
which case the singular ``also_notify_user`` argument silently dropped
the rename signal for every member. Broaden the helper to accept an
iterable of user ids and resolve group membership in ``update_room``
before fanning out.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Extend contract-sweep coverage (Important #4)

**Files:**
- Modify: `tests/zndraw/test_broadcast_contract.py:75-118` (add cases)

The existing sweep drives 6 endpoints. Add the missing families. The sweep's final loop already validates everything captured — adding endpoints is enough.

- [ ] **Step 1: Add more endpoint calls to the sweep**

In `tests/zndraw/test_broadcast_contract.py`, between the existing "joblib family" block and the `rs_event_names = ...` line, insert:

```python
# step family — FrameUpdate
r = await client.put(f"{base}/step", json={"step": 0}, headers=headers)
assert r.status_code == 200, r.text

# edit lock family — LockUpdate
r = await client.put(f"{base}/edit-lock", json={"action": "acquire"}, headers=headers)
assert r.status_code == 200, r.text

# geometry family — GeometryInvalidate
r = await client.put(
    f"{base}/geometries/g1",
    json={"geometry": {"type": "Sphere", "data": {}}},
    headers=headers,
)
assert r.status_code in (200, 201), r.text

# selection-groups family — SelectionGroupsInvalidate
r = await client.put(
    f"{base}/selection-groups/sg1",
    json={"selection": {}},
    headers=headers,
)
assert r.status_code in (200, 201), r.text

# frame-selection family — FrameSelectionUpdate (added in Task 2)
r = await client.put(
    f"{base}/frame-selection",
    json={"indices": [0]},
    headers=headers,
)
assert r.status_code == 200, r.text
```

Note: the exact endpoint shapes must match the routes. Verify by reading each route's request schema before running.

- [ ] **Step 2: Run the sweep**

```bash
uv run pytest tests/zndraw/test_broadcast_contract.py -v
```

Expected: PASS. Some routes may have different payload shapes — adjust the JSON until the route returns 2xx, but never weaken the final loop's assertions.

- [ ] **Step 3: If a route's request schema differs from the example, adjust**

For each 4xx the test returns, read the relevant route module's request schema (`src/zndraw/routes/<name>.py`) and update the JSON body. Do not weaken the post-loop assertions.

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw/test_broadcast_contract.py
git commit -m "$(cat <<'EOF'
test(broadcast): extend contract sweep across remaining families

Closes review finding #4. Drives PUT step, PUT edit-lock, PUT
geometry, PUT selection-group, and PUT frame-selection through the
existing contract sweep so every room-scoped event family is checked
end-to-end against ``room_id``, ``room_address``, and channel.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Frontend room_address gate (Important #6)

**Files:**
- Modify: `frontend/src/hooks/socketHandlers/utils.ts` (gate `createInvalidateHandler` on `data.room_address`)
- Modify: `frontend/src/hooks/socketHandlers/geometryHandlers.ts` (explicit gate on `onGeometriesInvalidate`)

- [ ] **Step 1: Update `createInvalidateHandler`**

Replace the body of `frontend/src/hooks/socketHandlers/utils.ts`:

```typescript
/**
 * Factory for creating consistent invalidate handlers.
 *
 * Fetches data from the server and updates the store.
 * Uses a `getRoomId` getter so the handler always reads the current roomId
 * at event-fire time (same behavior as the original closure).
 *
 * If the incoming event payload carries a `room_address` and it does not
 * match the current room, the handler is a no-op. This prevents stale
 * events from a previously-joined room mutating the current room's cache
 * during navigation races.
 */
export function createInvalidateHandler<T>(
	fetchFn: (roomId: string) => Promise<T>,
	updateStoreFn: (data: T) => void,
	eventName: string,
	getRoomId: () => string | undefined,
): (data: unknown) => Promise<void> {
	return async (data) => {
		const roomId = getRoomId();
		if (!roomId) return;
		if (data && typeof data === "object" && "room_address" in data) {
			const evtAddr = (data as { room_address?: unknown }).room_address;
			if (typeof evtAddr === "string" && evtAddr !== roomId) return;
		}
		try {
			const response = await fetchFn(roomId);
			updateStoreFn(response);
		} catch (error) {
			console.error(`Error fetching ${eventName}:`, error);
		}
	};
}
```

- [ ] **Step 2: Gate `onGeometriesInvalidate`**

Edit `frontend/src/hooks/socketHandlers/geometryHandlers.ts`. First widen the event interface near line 14:

```typescript
export interface GeometryInvalidateEvent {
	room_id?: string;
	room_address?: string;
	operation?: "set" | "delete";
	key?: string;
}
```

Then add the gate at the top of `onGeometriesInvalidate` (right after `if (!ctx.roomId) return;`):

```typescript
		if (!ctx.roomId) return;
		if (data.room_address && data.room_address !== ctx.roomId) return;
```

- [ ] **Step 3: Frontend type check + build**

```bash
cd frontend && bun run tsc --noEmit
```

Expected: no errors.

- [ ] **Step 4: Commit**

```bash
git add frontend/src/hooks/socketHandlers/utils.ts \
        frontend/src/hooks/socketHandlers/geometryHandlers.ts
git commit -m "$(cat <<'EOF'
fix(frontend): gate invalidate handlers on event room_address

Closes review finding #6. ``createInvalidateHandler`` now drops events
whose ``room_address`` does not match the current room, and
``onGeometriesInvalidate`` gates the same way. Prevents stale events
from a previously-joined room mutating the active room's cache during
navigation races.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Document schema-audit limitation (Important #7)

**Files:**
- Modify: `tests/zndraw/test_event_schema_audit.py:1-19` (expand module docstring)

- [ ] **Step 1: Expand the docstring**

Replace the module docstring (lines 1-5) of `tests/zndraw/test_event_schema_audit.py` with:

```python
"""Layer 1: static schema audit for RoomScopedEvent subclasses.

Fails fast at collect time if any subclass forgets to declare both
``room_id: UUID`` and ``room_address: str``.

LIMITATION: discovery uses ``RoomScopedEvent.__subclasses__()``, which
only enumerates classes that have been imported. The audit imports
``zndraw.socket_events`` and ``zndraw_joblib.events`` directly so every
first-party event type is registered before the test runs. If a third-
party plugin defines a ``RoomScopedEvent`` subclass in a module that
the audit does not import, that subclass slips through. Plugin authors
should either expose the module via an entry-point that the audit can
load explicitly or add their own equivalent test.
"""
```

- [ ] **Step 2: Run the audit to confirm it still passes**

```bash
uv run pytest tests/zndraw/test_event_schema_audit.py -v
```

Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add tests/zndraw/test_event_schema_audit.py
git commit -m "$(cat <<'EOF'
docs(events): note __subclasses__ limit in the schema audit

Closes review finding #7. ``RoomScopedEvent.__subclasses__()`` only
returns classes that have been imported, so plugin-defined events can
slip through unless their module is imported by the audit. Document
the constraint and the suggested workarounds.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Rename MockSioServer skip_sid (Minor #8)

**Files:**
- Modify: `tests/zndraw/helpers.py:178-197`

- [ ] **Step 1: Rename the parameter**

Edit `tests/zndraw/helpers.py:178-197`. Replace the `emit` method signature:

```python
    async def emit(
        self,
        event_or_model: str | BaseModel,
        data: Any = None,
        *,
        room: str | None = None,
        skip_sid: str | None = None,
        to: str | None = None,
        **_kwargs: Any,
    ) -> None:
```

And extend the captured dict (around line 197) so callers can inspect what was skipped:

```python
        self.emitted.append(
            {"event": event, "data": data, "room": room, "to": to, "skip_sid": skip_sid}
        )
```

- [ ] **Step 2: Run every test that uses MockSioServer**

```bash
uv run pytest tests/zndraw -v
```

Expected: PASS. (Existing tests don't read `skip_sid` from captured dicts, so the additional key is harmless.)

- [ ] **Step 3: Commit**

```bash
git add tests/zndraw/helpers.py
git commit -m "$(cat <<'EOF'
test(helpers): align MockSioServer.emit signature with production

Closes review finding #8. The mock had ``_skip_sid`` (silently absorbed
by ``**_kwargs``); production uses ``skip_sid``. Rename and capture the
value in ``emitted`` for test introspection.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Final verification (no commit)

- [ ] **Step 1: Run the entire pytest suite**

```bash
uv run pytest
```

Expected: every test PASSES.

- [ ] **Step 2: Run the linter and formatter**

```bash
uvx prek --all-files
```

Expected: no failures.

- [ ] **Step 3: Frontend type check**

```bash
cd frontend && bun run tsc --noEmit
```

Expected: clean.

- [ ] **Step 4: Sanity-check the commit log**

```bash
git log --oneline 4f21d6ccc5453e6c1447ea40def1d5e6788851c9..HEAD
```

Expected: 9 new commits on top of `43b143ad` (one per task except Task 10).

---

## Self-review checklist

- **Spec coverage:** All 11 findings have an assigned task or are explicitly waived (finding #11 — `room_id` field kept for documentation).
- **No placeholders:** every code block is complete; commands have expected output.
- **Type consistency:** `also_notify_user_ids` (Task 5) flows from `prev_user_ids` in `routes/rooms.py`; `build_room_scoped_emission` (Task 3) returns `Emission`, which `emit(...)` already accepts.
- **TDD honored:** every behavioral change has a red→green pair (Tasks 2, 3, 5).
- **Frequent commits:** 9 commits, one per task, no batching.
