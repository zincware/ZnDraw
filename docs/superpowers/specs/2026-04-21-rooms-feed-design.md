# Rooms feed: real-time room list updates

**Issue:** [#919](https://github.com/zincware/ZnDraw/issues/919) — Retire `/rooms` page and broadcast room updates to all visible viewers.

**Branch:** `fix/919-rooms-feed`

## Context

The redesigned sidebar panel (`RoomsPanel.tsx`) is a full replacement for
the standalone `/rooms` page. However, `room_update` events are only
broadcast to `room:@overview`, a channel only the (now unused) `/rooms`
page ever joined. As a result two distinct bugs manifest today:

1. **Cross-room staleness.** Client B viewing `/rooms/foo` never
   receives a `room_update` when client A creates `/rooms/bar` — the
   sidebar list on B stays stale until a manual refresh.
2. **Same-room staleness.** Two clients viewing the *same* `/rooms/foo`
   can show different frame counts in their sidebar row. Frame
   mutations emit `room_update` only to `@overview`, not to
   `room:{id}`, so even in-room viewers' `useRoomsStore` is never
   updated. The 3D viewer timeline is refreshed via
   `frames_invalidate` → `appStore.frameCount`, but that path does not
   touch `useRoomsStore`, which is what the sidebar row reads.

A third compounding issue: the emit fan-out is ad-hoc. `PATCH` emits to
both `room:{id}` and `@overview`; frame mutations emit only to
`@overview`; `set_default_room` emits to both. Centralizing the fan-out
is part of the fix.

## Non-goals

- Adding `DELETE /v1/rooms/{id}` or a `RoomDelete` socket event. The
  sidebar row menu already shows Delete as disabled with a "not yet
  supported by the backend" tooltip; that stays for a follow-up.
- Exposing private rooms in the public listing API (today
  `list_rooms` filters by `is_public=True`). The design honors
  private-room visibility rules but does not surface them to the UI.
- Any hotfix-style re-fetch. This design is the real fix.

## Design

### Channel topology

One new Socket.IO channel, `rooms:feed`, is entered in `on_connect`
(`src/zndraw/socketio.py:107`, right after `user:{user_id}`). Every
authenticated socket is a member for its whole lifetime — no explicit
join/leave, no client-side management. `@overview` is retired.

### Broadcast helper

A single helper in `src/zndraw/routes/rooms.py` centralizes fan-out.
Every `RoomUpdate` emit site calls this helper; no site emits directly.

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
    else:
        member_ids = (
            await session.exec(
                select(RoomMembership.user_id).where(
                    RoomMembership.room_id == room.id
                )
            )
        ).all()
        for uid in member_ids:
            await sio.emit(event, room=f"user:{uid}")
```

**Why no separate `room:{id}` emit.** Both scopes already cover
in-room viewers: the public branch because every authenticated socket
is in `rooms:feed`, the private branch because every in-room viewer
is a member and therefore in some `user:{uid}`. Emitting to
`room:{id}` on top would deliver the event twice to sockets in both
scopes, wasting a handler run (`setRoom` is idempotent so state is
still correct). The same-room staleness bug is fixed by reaching
in-room viewers through the visibility scope, not by a per-room emit.

### Call sites migrated to the helper

Each of these currently emits to `room:@overview` (and sometimes
`room:{id}` ad-hoc) and is replaced by a single
`await broadcast_room_update(sio, session, storage, room)`:

| File | Function | Current L# |
|---|---|---|
| `src/zndraw/routes/rooms.py` | `create_room` | 392–393 |
| `src/zndraw/routes/rooms.py` | `update_room` | 598–601 |
| `src/zndraw/routes/frames.py` | bulk append | 482–483 |
| `src/zndraw/routes/frames.py` | delete | 618–619 |
| `src/zndraw/routes/trajectory.py` | upload | 323–324 |
| `src/zndraw/routes/server_settings.py` | `set_default_room` (old + new) | 109–115 |
| `src/zndraw/routes/server_settings.py` | `unset_default_room` (old) | 144–146 |

After migration, no file in `src/zndraw/` contains the string
`@overview`.

### Frontend changes

**Delete:**

- `frontend/src/pages/roomList.tsx`.
- `path: "/rooms"` route + `RoomListPage` import in
  `frontend/src/App.tsx`.
- `handleGoToRoomList` callback and its menu item in
  `frontend/src/components/RoomManagementMenu.tsx` (L207-210,
  L389-394). Remove the `ListIcon` import if orphaned.

**Simplify:**

- `frontend/src/hooks/useSocketManager.ts`: drop the `isOverview`
  option from `SocketManagerOptions`; drop `isOverview` from
  `HandlerContext`; delete the `room_leave("@overview")` and
  `isNavigatingToOverview` branches at L166-178.
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts`: drop the
  `if (ctx.isOverview) { ... } else if (ctx.roomId) { ... }` split
  (L197-284). One code path remains: join `ctx.roomId` if set. The
  `rooms:feed` channel is auto-joined server-side; no client action
  needed.
- `frontend/src/hooks/socketHandlers/types.ts`: remove `isOverview`
  from `HandlerContext`.

**Unchanged:**

- `onRoomUpdate` in `frontend/src/hooks/socketHandlers/roomHandlers.ts`
  (L47) already upserts into `useRoomsStore`. It just starts
  receiving events for every public room.
- `socket.on("room_update", ...)` registration at
  `useSocketManager.ts:139`.

### Data flow after fix

```
  Client A: POST /v1/rooms/foo/frames/bulk (+22 frames)
      │
      ▼
  backend: storage.extend(...)          (write through to DB + storage)
  backend: sio.emit(FramesInvalidate,   (clients in room:foo update
           room=room:foo)                 the viewer's frame cache)
  backend: broadcast_room_update(...)
           → sio.emit(RoomUpdate, room=rooms:feed)     (every authenticated
                                                        client, incl. in-room)
      │
      ▼
  Client B in room:foo:
      onFramesInvalidate: invalidate react-query cache + setFrameCount
      onRoomUpdate:       useRoomsStore.setRoom(foo, {...frame_count: 44})
      → sidebar row + timeline both show 44 frames ✓

  Client C in room:bar (no action in foo):
      onRoomUpdate:       useRoomsStore.setRoom(foo, {...frame_count: 44})
      → sidebar row for foo updates live ✓
```

### Private-room visibility

Today `list_rooms` filters by `is_public=True`, so private rooms never
appear in any sidebar list. The design still emits private-room
`RoomUpdate`s to each member's `user:{uid}` channel. This is ~5 lines
in the helper; it future-proofs the invariant stated in the issue ("A
private room change must only reach its members") so we don't have to
retrofit it when private-room listing lands.

## Tests

Real Redis via `server_factory`; no mocks (per `CLAUDE.md`).

1. **Feed auto-join.** Two sockets connect; one joins no room, the
   other creates a public room via `POST /v1/rooms`. The idle client
   receives a `room_update` event without any explicit join.
2. **Same-room update** (regression for the screenshot bug). Two
   sockets join `room:foo`; one appends frames via
   `/frames/bulk`; both clients receive a `room_update` with
   `frame_count == new_total`.
3. **Private-room isolation.** User A creates a private room with
   user B as the only member; user C (authenticated, not a member)
   receives no `room_update` for that room. User B does.
4. **`@overview` retired.** Grep-style: no source file under
   `src/zndraw/` contains `@overview`. Existing
   `test_join_overview_room` is removed or rewritten to target
   `rooms:feed` auto-join.

## Dead-code sweep (executed before PR)

A final pass verifies no leftover references remain:

1. `grep -rn "@overview" src/ frontend/src/` → zero hits.
2. `grep -rn "isOverview" frontend/src/` → zero hits.
3. `grep -rn "room:@overview" src/ frontend/src/` → zero hits.
4. `grep -rn "RoomListPage\|pages/roomList" frontend/src/` → zero hits.
5. `grep -rn "handleGoToRoomList" frontend/src/` → zero hits.
6. `uv run ruff check` and `bun run lint` clean (catches unused
   imports automatically).
7. Docs / screenshot references to `/rooms` page updated or removed
   (per issue §1).

## Acceptance

- `/rooms` route returns 404.
- Creating, updating, or renaming a public room anywhere in the app
  updates the sidebar list in real time for every authenticated
  viewer — same room, different room, or no room.
- Within a single room, frame append/delete/upload updates the
  sidebar row's frame count for all in-room viewers without a page
  reload.
- Private rooms never emit `room_update` to non-members.
- No reference to `@overview` remains in source.
