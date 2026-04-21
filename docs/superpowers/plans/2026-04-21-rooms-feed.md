# Rooms feed implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `room:@overview` Socket.IO broadcast channel with a shared `rooms:feed` channel every authenticated socket auto-joins on connect, then retire the `/rooms` page that was the only `@overview` consumer. This fixes two user-visible bugs: cross-room staleness of the sidebar list, and same-room staleness of the per-row `frame_count` after frame mutations.

**Architecture:** Every authenticated socket enters `rooms:feed` in `on_connect` (right after `user:{user_id}`). A new helper `broadcast_room_update(sio, session, storage, room)` in `src/zndraw/routes/rooms.py` centralizes all `RoomUpdate` fan-out: public rooms → `rooms:feed`, private rooms → per-member `user:{uid}`. The helper replaces every direct `room="room:@overview"` emit across `rooms.py`, `frames.py`, `trajectory.py`, `server_settings.py`. The frontend's `/rooms` page, the `isOverview` prop threading through `useSocketManager.ts` + `connectionHandlers.ts`, and the "Go to Room List" menu item are deleted. No frontend socket-handler code changes — `onRoomUpdate` in `roomHandlers.ts` already upserts into `useRoomsStore`; it simply starts receiving events for every visible room.

**Tech Stack:** Python 3.11+ / FastAPI / python-socketio via zndraw-socketio / SQLModel. React 18 / TypeScript / Zustand. `uv` for Python deps, `bun` for TypeScript. Real Redis via TcpFakeServer in tests.

**Branch:** `fix/919-rooms-feed` (already created).

**Spec:** `docs/superpowers/specs/2026-04-21-rooms-feed-design.md`.

---

## File structure

### Backend (modify)

- `src/zndraw/socketio.py` — `on_connect` auto-joins `rooms:feed`.
- `src/zndraw/routes/rooms.py` — add `broadcast_room_update` helper; migrate `create_room` + `update_room`.
- `src/zndraw/routes/frames.py` — migrate bulk-append + delete to helper.
- `src/zndraw/routes/trajectory.py` — migrate upload to helper.
- `src/zndraw/routes/server_settings.py` — migrate `set_default_room` + `unset_default_room` to helper.

### Frontend (delete)

- `frontend/src/pages/roomList.tsx` — the page itself.

### Frontend (modify)

- `frontend/src/App.tsx` — drop `/rooms` route + import.
- `frontend/src/components/RoomManagementMenu.tsx` — drop "Go to Room List" menu item + its handler + orphaned import.
- `frontend/src/hooks/socketHandlers/types.ts` — drop `isOverview` from `HandlerContext`.
- `frontend/src/hooks/useSocketManager.ts` — drop `isOverview` option, context field, and `@overview` leave branch.
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts` — collapse the `isOverview` / `roomId` split into a single `roomId`-only branch.

### Tests

- `tests/zndraw/test_socketio_rooms.py` — delete `test_socketio_join_system_room_overview`; add `test_rooms_feed_auto_join`, `test_same_room_frame_append_updates_sidebar`.
- `tests/zndraw/test_broadcast_room_update.py` (new) — unit tests of the helper's routing using `MockSioServer`, including the private-room exclusion invariant (no REST path exists to make a private room, so this is the only place it can be tested today).

---

## Task 1: Add `rooms:feed` auto-join on connect (TDD)

**Files:**
- Modify: `src/zndraw/socketio.py:83-108` (`on_connect` handler)
- Test: `tests/zndraw/test_socketio_rooms.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/zndraw/test_socketio_rooms.py`:

```python
@pytest.mark.asyncio
async def test_rooms_feed_auto_join(
    server: str, http_client: AsyncClient
) -> None:
    """Every authenticated socket must auto-join rooms:feed on connect,
    so room_update events from any public room are delivered without
    an explicit join.

    Regression test for #919: an idle client (in no room) receives a
    room_update when another user creates a public room.
    """
    from zndraw.socket_events import RoomUpdate

    token_a = await _get_user_token(http_client, "feed-a@example.com")
    token_b = await _get_user_token(http_client, "feed-b@example.com")

    # Client B connects but joins NO room.
    sio_b = socketio.AsyncClient()
    received: list[RoomUpdate] = []

    @sio_b.on("room_update")
    async def on_room_update(data: dict) -> None:
        received.append(RoomUpdate.model_validate(data))

    await sio_b.connect(server, auth={"token": token_b})

    # Client A creates a public room via REST.
    new_room_id = await _create_room(http_client, token_a)

    # Give the event loop a beat to deliver the broadcast.
    await asyncio.sleep(0.5)

    assert any(
        e.id == new_room_id for e in received
    ), f"Client B received no room_update for {new_room_id}; got {received}"

    await sio_b.disconnect()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/zndraw/test_socketio_rooms.py::test_rooms_feed_auto_join -v`

Expected: FAIL — the `room_update` emit still targets `room:@overview`, which nothing joins.

- [ ] **Step 3: Add the auto-join**

Edit `src/zndraw/socketio.py`:

```python
# Replace the single-line enter_room block in on_connect with:
    user_id = UUID(payload["sub"])
    await tsio.save_session(sid, {"user_id": user_id, "current_room_id": None})
    await tsio.enter_room(sid, f"user:{user_id}")
    await tsio.enter_room(sid, "rooms:feed")
    return True
```

- [ ] **Step 4: Temporarily redirect one emit site so the test can pass**

In `src/zndraw/routes/rooms.py` `create_room` (currently around L393), change:

```python
    await sio.emit(event, room="room:@overview")
```

to:

```python
    await sio.emit(event, room="rooms:feed")
```

(The other five emit sites get migrated in Task 3. This minimal change keeps Task 1 focused.)

- [ ] **Step 5: Run test to verify it passes**

Run: `uv run pytest tests/zndraw/test_socketio_rooms.py::test_rooms_feed_auto_join -v`

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/socketio.py src/zndraw/routes/rooms.py tests/zndraw/test_socketio_rooms.py
git commit -m "feat(socketio): auto-join rooms:feed for live room_update delivery"
```

---

## Task 2: `broadcast_room_update` helper + unit tests (TDD)

**Files:**
- Modify: `src/zndraw/routes/rooms.py` (add helper next to `build_room_update`)
- Test: `tests/zndraw/test_broadcast_room_update.py` (new)

- [ ] **Step 1: Write the failing tests**

Create `tests/zndraw/test_broadcast_room_update.py`:

```python
"""Unit tests for broadcast_room_update helper routing logic.

These tests cover cases that cannot be driven via REST today:
private rooms (no REST way to create them) and the exact set of
channels a broadcast targets. Real socket integration is covered by
test_socketio_rooms.py.
"""

from uuid import uuid4

import pytest
from helpers import MockSioServer
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.models import MemberRole, Room, RoomMembership
from zndraw.routes.rooms import broadcast_room_update
from zndraw.storage import FrameStorage


@pytest.mark.asyncio
async def test_broadcast_public_room_targets_feed(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Public rooms broadcast to the shared rooms:feed channel only."""
    room = Room(id="pub", is_public=True)
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = [call["room"] for call in sio.emitted]
    assert rooms_targeted == ["rooms:feed"]


@pytest.mark.asyncio
async def test_broadcast_private_room_targets_each_member(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Private rooms broadcast to each member's user:{uid} channel and
    nothing else — never rooms:feed, never non-members."""
    member_a = uuid4()
    member_b = uuid4()
    non_member = uuid4()  # noqa: F841 — intentionally unreferenced

    room = Room(id="priv", is_public=False)
    session.add(room)
    session.add(
        RoomMembership(room_id="priv", user_id=member_a, role=MemberRole.OWNER)
    )
    session.add(
        RoomMembership(room_id="priv", user_id=member_b, role=MemberRole.MEMBER)
    )
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    rooms_targeted = sorted(call["room"] for call in sio.emitted)
    assert rooms_targeted == sorted(
        [f"user:{member_a}", f"user:{member_b}"]
    )
    assert "rooms:feed" not in rooms_targeted


@pytest.mark.asyncio
async def test_broadcast_private_room_with_no_members_is_noop(
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Private room with no memberships emits nothing — no channel is
    authorized, so no broadcast occurs."""
    room = Room(id="orphan", is_public=False)
    session.add(room)
    await session.commit()

    sio = MockSioServer()
    await broadcast_room_update(sio, session, frame_storage, room)

    assert sio.emitted == []
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/zndraw/test_broadcast_room_update.py -v`

Expected: FAIL with `ImportError: cannot import name 'broadcast_room_update'`.

- [ ] **Step 3: Implement the helper**

In `src/zndraw/routes/rooms.py`, immediately after the existing `build_room_update` function (around L281), add:

```python
async def broadcast_room_update(
    sio,
    session: AsyncSession,
    storage: FrameStorage,
    room: Room,
) -> None:
    """Broadcast a full RoomUpdate to every authorized viewer.

    Public rooms go to the shared ``rooms:feed`` channel — every
    authenticated socket is a member, including any client currently
    joined to ``room:{id}``, so this single emit reaches in-room
    viewers too.

    Private rooms fan out to each member's ``user:{uid}`` channel,
    which similarly covers both in-room and out-of-room members.
    """
    event = await build_room_update(session, storage, room)
    if room.is_public:
        await sio.emit(event, room="rooms:feed")
        return
    result = await session.exec(
        select(RoomMembership.user_id).where(
            RoomMembership.room_id == room.id
        )
    )
    for uid in result.all():
        await sio.emit(event, room=f"user:{uid}")
```

Also add `RoomMembership` to the imports near the top of the file:

```python
from zndraw.models import (
    Room,
    RoomBookmark,
    RoomFigure,
    RoomGeometry,
    RoomMembership,
    SelectionGroup,
    ServerSettings,
)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/zndraw/test_broadcast_room_update.py -v`

Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/rooms.py tests/zndraw/test_broadcast_room_update.py
git commit -m "feat(rooms): broadcast_room_update helper routes by visibility"
```

---

## Task 3: Migrate every emit site to the helper

Each step migrates one file. The Task 1 temporary redirect in `create_room` is superseded by the helper call.

**Files:**
- Modify: `src/zndraw/routes/rooms.py` (two emit sites)
- Modify: `src/zndraw/routes/frames.py` (two emit sites)
- Modify: `src/zndraw/routes/trajectory.py` (one emit site)
- Modify: `src/zndraw/routes/server_settings.py` (three emit sites across two endpoints)

- [ ] **Step 1: Migrate `src/zndraw/routes/rooms.py`**

In `create_room` (currently around L391-393), replace:

```python
    # Broadcast room creation to @overview system room
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="rooms:feed")
```

with:

```python
    await broadcast_room_update(sio, session, storage, room)
```

In `update_room` (currently around L597-601), replace the whole `if changed:` branch's emit block:

```python
    if changed:
        event = await build_room_update(session, storage, room)
        log.debug("Broadcasting RoomUpdate: %s", event.model_dump())
        await sio.emit(event, room=f"room:{room.id}")
        await sio.emit(event, room="room:@overview")
```

with:

```python
    if changed:
        await broadcast_room_update(sio, session, storage, room)
```

The old block had a `log.debug(...)` line. `log` is defined at L74
(`log = logging.getLogger(__name__)`) and was only used by that one
debug call. After the replacement, delete both the `log = ...` line
(L74) and the `import logging` at the top of the file. `ruff` would
flag these as unused if left in place.

- [ ] **Step 2: Migrate `src/zndraw/routes/frames.py`**

In the bulk-append endpoint (currently L481-483), replace:

```python
    # Notify @overview of updated frame count
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")
```

with:

```python
    await broadcast_room_update(sio, session, storage, room)
```

Update the import on L44 from:

```python
from zndraw.routes.rooms import build_room_update
```

to:

```python
from zndraw.routes.rooms import broadcast_room_update
```

Apply the same emit-block replacement in the delete endpoint (currently L617-619).

Confirm `build_room_update` is no longer referenced in `frames.py` (it shouldn't be after this step):

```bash
grep -n "build_room_update" src/zndraw/routes/frames.py
```

Expected: no matches.

- [ ] **Step 3: Migrate `src/zndraw/routes/trajectory.py`**

In the upload endpoint (currently L322-324), replace:

```python
    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")
```

with:

```python
    await broadcast_room_update(sio, session, storage, room)
```

Update the import on L42 from `build_room_update` to `broadcast_room_update`, same as Task 3 Step 2.

- [ ] **Step 4: Migrate `src/zndraw/routes/server_settings.py`**

In `set_default_room` (currently L105-115), replace the whole block:

```python
    # Broadcast full snapshots for affected rooms
    if old_default_id and old_default_id != request.room_id:
        old_room = await session.get(Room, old_default_id)
        if old_room:
            event = await build_room_update(session, storage, old_room)
            await sio.emit(event, room="room:@overview")
            await sio.emit(event, room=f"room:{old_default_id}")

    event = await build_room_update(session, storage, room)
    await sio.emit(event, room="room:@overview")
    await sio.emit(event, room=f"room:{request.room_id}")
```

with:

```python
    # Broadcast full snapshots for affected rooms
    if old_default_id and old_default_id != request.room_id:
        old_room = await session.get(Room, old_default_id)
        if old_room:
            await broadcast_room_update(sio, session, storage, old_room)

    await broadcast_room_update(sio, session, storage, room)
```

In `unset_default_room` (currently L142-146), replace:

```python
        old_room = await session.get(Room, old_default_id)
        if old_room:
            event = await build_room_update(session, storage, old_room)
            await sio.emit(event, room="room:@overview")
            await sio.emit(event, room=f"room:{old_default_id}")
```

with:

```python
        old_room = await session.get(Room, old_default_id)
        if old_room:
            await broadcast_room_update(sio, session, storage, old_room)
```

Update the import on L20 from `build_room_update` to `broadcast_room_update`.

- [ ] **Step 5: Verify the `@overview` string is gone from `src/zndraw/`**

```bash
grep -rn "@overview" src/zndraw/
```

Expected: no matches.

- [ ] **Step 6: Verify `build_room_update` now has exactly one importer (rooms.py itself, internal call)**

```bash
grep -rn "build_room_update" src/zndraw/
```

Expected: only the definition and the one internal call inside `broadcast_room_update`, both in `src/zndraw/routes/rooms.py`. No imports from other route files.

- [ ] **Step 7: Run the existing backend test suite to catch regressions**

Run: `uv run pytest tests/zndraw/test_routes_rooms.py tests/zndraw/test_routes_frames.py tests/zndraw/test_trajectory.py tests/zndraw/test_socketio_rooms.py tests/zndraw/test_broadcast_room_update.py -v`

Expected: all pass. If `test_socketio_join_system_room_overview` fails because the emit no longer lands there, don't fix that test here — Task 5 removes it.

- [ ] **Step 8: Commit**

```bash
git add src/zndraw/routes/
git commit -m "refactor(rooms): migrate every RoomUpdate emit to broadcast_room_update"
```

---

## Task 4: Same-room `frame_count` live update (regression test for the screenshot bug)

**Files:**
- Test: `tests/zndraw/test_socketio_rooms.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/zndraw/test_socketio_rooms.py`:

```python
@pytest.mark.asyncio
async def test_same_room_frame_append_updates_sidebar(
    server: str, http_client: AsyncClient
) -> None:
    """Two clients viewing the same room must both receive a room_update
    with the new frame_count when one of them appends frames.

    Regression test for the screenshot in #919: two tabs on /rooms/s22
    showed 44 frames and 22 frames respectively because frame mutations
    only emitted room_update to @overview, never to room:{id}.
    """
    from zndraw.socket_events import RoomUpdate

    token_a = await _get_user_token(http_client, "same-room-a@example.com")
    token_b = await _get_user_token(http_client, "same-room-b@example.com")

    room_id = await _create_room(http_client, token_a)

    # Both clients join the same room.
    sio_a = socketio.AsyncClient()
    sio_b = socketio.AsyncClient()

    received_b: list[RoomUpdate] = []

    @sio_b.on("room_update")
    async def on_room_update(data: dict) -> None:
        received_b.append(RoomUpdate.model_validate(data))

    await sio_a.connect(server, auth={"token": token_a})
    await sio_b.connect(server, auth={"token": token_b})

    tsio_a = wrap(sio_a)
    tsio_b = wrap(sio_b)
    await tsio_a.call(
        RoomJoin(room_id=room_id), response_model=RoomJoinResponse
    )
    await tsio_b.call(
        RoomJoin(room_id=room_id), response_model=RoomJoinResponse
    )

    # Drain any room_update emissions that predate the frame append.
    await asyncio.sleep(0.3)
    received_b.clear()

    # Client A appends frames via REST. The endpoint validates required
    # rendering keys (arrays.colors, arrays.radii), so we build proper
    # frames via the serialization helper — the same path test_routes_frames
    # uses. See src/zndraw/routes/frames.py:64 for the validation.
    import ase  # local import keeps top-of-file imports untouched

    from zndraw.client import atoms_to_json_dict

    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    frame = atoms_to_json_dict(atoms)

    append_response = await http_client.post(
        f"/v1/rooms/{room_id}/frames",
        json={"frames": [frame, frame, frame]},
        headers={"Authorization": f"Bearer {token_a}"},
    )
    assert append_response.status_code == 201, append_response.text

    await asyncio.sleep(0.5)

    # Client B must see a room_update with the new total.
    assert received_b, "Client B received no room_update"
    latest = received_b[-1]
    assert latest.id == room_id
    assert latest.frame_count >= 3

    await sio_a.disconnect()
    await sio_b.disconnect()
```

- [ ] **Step 2: Run the test**

Run: `uv run pytest tests/zndraw/test_socketio_rooms.py::test_same_room_frame_append_updates_sidebar -v`

Expected: PASS (Task 3 already wired `frames.py` to the helper, which emits to `rooms:feed`; client B is in that feed via auto-join from Task 1). If it fails, the bulk-append URL or request body differs — check `src/zndraw/routes/frames.py` for the exact endpoint path and expected frame shape, then adjust the request.

- [ ] **Step 3: Commit**

```bash
git add tests/zndraw/test_socketio_rooms.py
git commit -m "test(socketio): assert same-room room_update on frame append"
```

---

## Task 5: Retire `@overview` from the test suite

**Files:**
- Modify: `tests/zndraw/test_socketio_rooms.py`

- [ ] **Step 1: Remove the obsolete test**

Delete the whole function `test_socketio_join_system_room_overview` in `tests/zndraw/test_socketio_rooms.py` (currently L166-188). The `rooms:feed` auto-join test from Task 1 covers the replacement semantics.

- [ ] **Step 2: Verify no remaining `@overview` string in tests**

```bash
grep -rn "@overview" tests/zndraw/
```

Expected: one match in `tests/zndraw/test_routes_rooms.py:78` — the REST room-ID validation test asserting `@overview` is rejected as a room_id (unrelated to socket broadcasts; this test must stay).

- [ ] **Step 3: Run the full socketio test file**

Run: `uv run pytest tests/zndraw/test_socketio_rooms.py -v`

Expected: all pass.

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw/test_socketio_rooms.py
git commit -m "test: drop obsolete @overview join test"
```

---

## Task 6: Delete `/rooms` page and route

**Files:**
- Delete: `frontend/src/pages/roomList.tsx`
- Modify: `frontend/src/App.tsx`

- [ ] **Step 1: Delete the page**

```bash
git rm frontend/src/pages/roomList.tsx
```

- [ ] **Step 2: Update `frontend/src/App.tsx`**

Remove the import on L13:

```tsx
import RoomListPage from "./pages/roomList";
```

Remove the route entry from the `createBrowserRouter` call:

```tsx
		{
			path: "/rooms",
			element: <RoomListPage />,
		},
```

Leave the other routes unchanged.

- [ ] **Step 3: Verify no remaining reference**

```bash
grep -rn "RoomListPage\|pages/roomList" frontend/src/
```

Expected: no matches.

- [ ] **Step 4: Typecheck**

Run from `frontend/`:

```bash
cd frontend && bun run tsc --noEmit && cd ..
```

Expected: no errors.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/App.tsx
git commit -m "feat(frontend): remove /rooms page; RoomsPanel fully replaces it"
```

---

## Task 7: Remove "Go to Room List" from `RoomManagementMenu`

**Files:**
- Modify: `frontend/src/components/RoomManagementMenu.tsx`

- [ ] **Step 1: Remove the handler and the menu item**

Delete the `handleGoToRoomList` function (currently L207-210):

```tsx
	const handleGoToRoomList = () => {
		navigate("/rooms");
		handleCloseMenu();
	};
```

Delete the menu item that calls it (currently L389-394):

```tsx
				<MenuItem onClick={handleGoToRoomList}>
					<ListItemIcon>
						<ListIcon />
					</ListItemIcon>
					<ListItemText>Go to Room List</ListItemText>
				</MenuItem>
```

- [ ] **Step 2: Check if `ListIcon` is still imported**

```bash
grep -n "ListIcon" frontend/src/components/RoomManagementMenu.tsx
```

If the only reference remaining is the import line (e.g., `import ListIcon from "@mui/icons-material/List";`), remove that import too. The `tsc` check in step 5 would catch this as an unused symbol, but removing it now keeps the commit clean.

- [ ] **Step 3: Verify `navigate` is still used**

If `handleGoToRoomList` was the only consumer of `navigate`, remove the `useNavigate` hook and its import as well. Check:

```bash
grep -n "navigate(" frontend/src/components/RoomManagementMenu.tsx
```

If there are other `navigate(...)` calls (e.g., `handleGoToFilesystem` at L214), keep the hook.

- [ ] **Step 4: Typecheck**

```bash
cd frontend && bun run tsc --noEmit && cd ..
```

Expected: no errors.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/components/RoomManagementMenu.tsx
git commit -m "feat(frontend): remove \"Go to Room List\" menu item"
```

---

## Task 8: Drop `isOverview` from socket-handler types

**Files:**
- Modify: `frontend/src/hooks/socketHandlers/types.ts`

- [ ] **Step 1: Inspect the current shape**

```bash
grep -n "isOverview" frontend/src/hooks/socketHandlers/types.ts
```

If present, remove the field from the `HandlerContext` interface. The field sits alongside `roomId`, `appStoreRoomId`, etc.

Example — before:

```ts
export interface HandlerContext {
	roomId?: string;
	appStoreRoomId?: string;
	isOverview: boolean;
	// ...
}
```

After:

```ts
export interface HandlerContext {
	roomId?: string;
	appStoreRoomId?: string;
	// ...
}
```

- [ ] **Step 2: Typecheck (expect errors from consumers)**

```bash
cd frontend && bun run tsc --noEmit && cd ..
```

Expected: errors in `useSocketManager.ts` and `connectionHandlers.ts` referencing `isOverview`. Task 9 and Task 10 clear them.

---

## Task 9: Remove `isOverview` from `useSocketManager.ts`

**Files:**
- Modify: `frontend/src/hooks/useSocketManager.ts`

- [ ] **Step 1: Drop the option from `SocketManagerOptions`**

Change (currently L16-19):

```ts
interface SocketManagerOptions {
	roomId?: string; // Room ID when on /rooms/:roomId page
	isOverview?: boolean; // True when on /rooms page
}
```

to:

```ts
interface SocketManagerOptions {
	roomId?: string; // Room ID when on /rooms/:roomId page
}
```

- [ ] **Step 2: Drop the destructure on L61**

Remove:

```ts
	const { isOverview = false } = options;
```

- [ ] **Step 3: Drop `effectIsOverview`, `isOverview` from the context, and the overview-leave branch**

Remove `effectIsOverview` (L66) and the `isOverview` entry from the `HandlerContext` object (around L75). Replace the whole unmount cleanup block (L165-185) that handles `isNavigatingToOverview` and `isLeavingRoom` with the simplified version:

```ts
			return () => {
				cancelled = true;
				// Note: room_leave is NOT emitted during room switching
				// Backend automatically handles leaving old room in room_join handler
				const currentState = useAppStore.getState();
				const isLeavingRoom = effectRoomId !== currentState.roomId;

				if (isLeavingRoom) {
					setSessionId(null);
					setCameraKey(null);
				}

				chatCleanup();
```

(Keep every `socket.off(...)` call below that block unchanged.)

- [ ] **Step 4: Drop `isOverview` from the `useEffect` dependency array**

Around L223-252, remove the standalone `isOverview,` line from the dependency list.

- [ ] **Step 5: Typecheck**

```bash
cd frontend && bun run tsc --noEmit && cd ..
```

Expected: the only remaining error is in `connectionHandlers.ts` (handled in Task 10).

---

## Task 10: Collapse the `onConnect` branch in `connectionHandlers.ts`

**Files:**
- Modify: `frontend/src/hooks/socketHandlers/connectionHandlers.ts`

- [ ] **Step 1: Replace the dispatch block**

The current dispatch (around L197-284) has three branches: `isOverview`, `roomId`, and the `else` fallback. Replace the whole `// Join appropriate room based on page` through the end of `onConnect` body with:

```ts
			// Join the specific room if one is set; otherwise connect is
			// sufficient — rooms:feed is auto-joined server-side.
			if (ctx.roomId) {
				socket.emit(
					"room_join",
					{ room_id: ctx.roomId, client_type: "frontend" },
					async (response: RoomJoinResponse | RoomJoinError) => {
						// Handle 404 - room doesn't exist, create it via REST API
						if ("status" in response && response.status === 404) {
							const urlCopyFrom = new URLSearchParams(
								window.location.search,
							).get("copy_from");
							try {
								await createRoom({
									room_id: ctx.roomId!,
									copy_from: urlCopyFrom ?? undefined,
								});
							} catch (error: any) {
								// 409 Conflict = room created by another client, continue
								if (error.response?.status !== 409) {
									console.error("Failed to create room:", error);
									ctx.setInitializationError({
										message: "Failed to create room",
										details:
											error instanceof Error
												? error.message
												: "Could not create room on server",
									});
									return;
								}
							}
							// Retry join after room creation
							socket.emit(
								"room_join",
								{ room_id: ctx.roomId!, client_type: "frontend" },
								(retryResponse: RoomJoinResponse | RoomJoinError) => {
									if ("status" in retryResponse) {
										console.error(
											"Failed to join room:",
											retryResponse,
										);
										ctx.setInitializationError({
											message: "Failed to join room",
											details:
												retryResponse.detail ||
												"Server rejected the connection",
										});
										return;
									}
									handleRoomJoin(retryResponse);
								},
							);
							return;
						}

						if ("status" in response) {
							console.error("Failed to join room:", response);
							ctx.setInitializationError({
								message: "Failed to join room",
								details:
									response.detail || "Server rejected the connection",
							});
							return;
						}

						handleRoomJoin(response);
					},
				);
			} else {
				ctx.setConnected(true);
			}
```

The difference from the current code is: the `isOverview` branch is gone and the `else if (ctx.roomId)` has been flattened to `if (ctx.roomId)`.

- [ ] **Step 2: Typecheck**

```bash
cd frontend && bun run tsc --noEmit && cd ..
```

Expected: no errors.

- [ ] **Step 3: Lint**

```bash
cd frontend && bun run lint && cd ..
```

Expected: clean (no unused imports, no dangling variables).

- [ ] **Step 4: Commit Tasks 8 + 9 + 10 together**

```bash
git add frontend/src/hooks/socketHandlers/types.ts \
        frontend/src/hooks/useSocketManager.ts \
        frontend/src/hooks/socketHandlers/connectionHandlers.ts
git commit -m "feat(frontend): drop isOverview plumbing; rooms:feed is auto-joined"
```

---

## Task 11: Dead-code sweep

Final verification pass. Any nonzero result here is a mistake to fix before the PR.

- [ ] **Step 1: Search for `@overview` across the whole repo**

```bash
grep -rn "@overview" src/ frontend/src/ tests/
```

Expected: exactly one match — `tests/zndraw/test_routes_rooms.py:78` (the room-ID validation parametrize case). Anything else in `src/`, `frontend/src/`, or docs is a leftover.

- [ ] **Step 2: Search for `isOverview` in the frontend**

```bash
grep -rn "isOverview" frontend/src/
```

Expected: no matches.

- [ ] **Step 3: Search for `RoomListPage` / `pages/roomList`**

```bash
grep -rn "RoomListPage\|pages/roomList" frontend/src/
```

Expected: no matches.

- [ ] **Step 4: Search for the old per-channel string literal**

```bash
grep -rn "room:@overview\|\"@overview\"\|'@overview'" src/ frontend/src/ tests/
```

Expected: only the `"@overview"` literal inside the parametrize list in `tests/zndraw/test_routes_rooms.py:78`.

- [ ] **Step 5: Search for the obsolete `handleGoToRoomList`**

```bash
grep -rn "handleGoToRoomList" frontend/src/
```

Expected: no matches.

- [ ] **Step 6: Ruff + pyright + bun lint**

```bash
uvx prek --all-files
```

```bash
cd frontend && bun run lint && bun run tsc --noEmit && cd ..
```

Expected: all clean.

- [ ] **Step 7: Full backend test run**

```bash
uv run pytest tests/zndraw/ -x
```

Expected: all pass. If any test references `@overview` in a way not covered above (e.g., comment references in unrelated files), decide inline whether to delete the comment or leave it; prefer deletion when the comment no longer describes real behavior.

- [ ] **Step 8: Manual two-tab verification (matches the bug report)**

From the repo root:

1. Start the dev server: `uv run zndraw`.
2. In a second terminal, start the frontend dev server: `cd frontend && bun run dev`.
3. Open `http://localhost:5173/rooms/s22` in two different tabs.
4. In tab A, upload a small trajectory file via the sidebar upload button.
5. In tab B, within ~1 second the sidebar row for `s22` must update to the new `frame_count`. The timeline at the bottom must also update.
6. Open tab C on `http://localhost:5173/rooms/other`. Create a new room in tab A via the `+` button — tab C's sidebar must show the new room without a refresh.

If any step fails, the broadcast path has a gap — inspect `broadcast_room_update` call sites and the `rooms:feed` auto-join.

- [ ] **Step 9: Commit any sweep fix-ups, then push**

If the sweep surfaced anything:

```bash
git add -u
git commit -m "chore: clean up residual @overview references"
```

Push the branch:

```bash
git push -u origin fix/919-rooms-feed
```

---

## Self-review notes

**Spec coverage mapping:**

- Spec §Channel topology → Task 1 (`on_connect` auto-join).
- Spec §Broadcast helper → Task 2.
- Spec §Call sites migrated → Task 3 (every site named in the spec's table).
- Spec §Frontend changes (delete) → Task 6, Task 7.
- Spec §Frontend changes (simplify) → Tasks 8–10.
- Spec §Data flow → Task 4 (same-room) + Task 1 (cross-room) regression tests.
- Spec §Private-room visibility → Task 2 (two of three unit tests).
- Spec §Tests — feed auto-join → Task 1; same-room update → Task 4; private-room isolation → Task 2; `@overview` retired → Task 5 + Task 11 sweep.
- Spec §Dead-code sweep → Task 11.

**Known limitations encoded in the plan:**

- Private-room tests are unit tests against `MockSioServer`, not integration tests. There is no REST path today to create a private room (`RoomCreate` has no `is_public` field and no admin PATCH exposes it). This is flagged in the plan's File structure section and is the single deviation from the "prefer `server_factory` over mocks" preference — it's forced by missing API surface, not preference.
