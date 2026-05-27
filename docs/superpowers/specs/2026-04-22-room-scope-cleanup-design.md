# Room scope refactor — follow-up cleanup

**Branch:** `feat/room-group-scope-refactor` (continuation)

## Context

Three loose ends remain after the room / group scope refactor.

1. `superuserLock` state in the frontend store is never written to `true`
   after `Room.locked` was removed. It is dead weight across a store slice,
   a selector, a component, and the socket-handler context.
2. `RoomJoinResponse.locked: bool` was kept as a hardcoded `False` for
   "backward compat" but nothing downstream still reads it meaningfully.
3. The frontend room-join flow auto-creates a room whenever the server
   returns 404 on `room_join`. This silently conflates two distinct
   outcomes: the room genuinely does not exist vs. the room exists but
   the caller cannot `can_read`. The second case, by design, returns 404
   to hide existence — but the auto-create retry surfaces the resulting
   bug as a user-visible "Failed to join room / Room with id X not
   found" snackbar.

A separate gap — the `/groups` frontend page is missing member
management, role changes, group deletion, share-link UI, and room
ownership transfer affordances — is acknowledged but tracked as its
own follow-up, not covered here.

## Non-goals

- New group management UI. The `/groups` page gaps are separate work.
- Any backend access-control changes. `can_read` / `can_edit` /
  `can_manage` stay as they are.
- Preserving the "paste a URL → fresh room" convenience. The feature is
  an existence-leak channel (a URL that succeeds tells the caller the
  room did not exist; a URL that fails conflates absence with denial).
  It is removed intentionally.

## Design

### Change 1 — remove `superuserLock` state

Touch points:

- `frontend/src/stores/slices/lockSlice.ts` — drop `superuserLock` field
  and `setSuperuserLock` action from the slice.
- `frontend/src/store.tsx` — `selectIsRoomReadOnly` loses its
  `superuserLock` branch; the selector still checks the superuser bypass
  and the `userLock`.
- `frontend/src/components/geometry/GeometryGrid.tsx` — drop the
  `superuserLock` selector, the corresponding branch in `canEdit`
  (line 77), and the "Room is locked" tooltip fallback in
  `editabilityTooltip` (line 141).
- `frontend/src/hooks/useSocketManager.ts` — drop the `setSuperuserLock`
  capture and pass-through (lines 46, 93, 227).
- `frontend/src/hooks/socketHandlers/types.ts` — drop `setSuperuserLock`
  from the `HandlerContext` type.
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts` — remove
  the dead write at lines 68-69 (`ctx.setSuperuserLock(response.locked ?? false)`).

### Change 2 — drop `RoomJoinResponse.locked`

Backend:

- `src/zndraw/socket_events.py` — remove `locked: bool` from the
  `RoomJoinResponse` model.
- `src/zndraw/socketio.py` — remove `locked=False` from both
  `RoomJoinResponse(...)` constructions (system-room path and normal
  path).

Frontend:

- `frontend/src/hooks/socketHandlers/connectionHandlers.ts` — remove
  `locked` from the local `RoomJoinResponse` interface declaration.

Tests:

- `tests/zndraw/test_socket_commands.py` — remove `locked=...` from the
  two `RoomJoinResponse(...)` constructions (full-field and
  required-fields tests) and drop the `resp.locked == ...` assertion.

No other producer or consumer of the field exists. The pyclient
(`src/zndraw/client/socket.py`) validates the shape through
`RoomJoinResponse.model_validate` but never reads `locked`, so dropping
it is transparent there once the model is updated.

### Change 3 — close the room-create existence leak + drop URL auto-create

Two coordinated changes. The backend change is required on its own —
the existing POST handler leaks existence directly. The frontend
change is required on its own — the auto-create retry turns any
`can_read` denial into a confusing "not found" snackbar.

**Backend — `POST /v1/rooms` always suffixes the `room_id`.**

Current behavior (`src/zndraw/routes/rooms.py:360-369`) returns 200
with the existing room's `frame_count` whenever the requested
`room_id` already exists, *regardless of whether the caller can_read
the existing row*. A caller can therefore probe any `room_id` and
learn whether it exists (201 = new, 200 = taken). This bypasses the
`404-hides-existence` contract that `GET /v1/rooms/{id}` and
`room_join` already honor, and it also leaks `frame_count` for
private rooms.

A two-tier fix (suffix only on unreadable collision, echo verbatim
otherwise) is insufficient: a caller can still compare the returned
`room_id` against the requested value to distinguish "free" (echoed
verbatim) from "taken-but-hidden" (suffixed). Any conditional
suffixing is a side-channel.

New behavior: **every successful `POST /v1/rooms` returns a new room
at `{request.room_id}-{6 random url-safe chars}`.** No idempotent
"return existing" branch, no verbatim echo, no conditionals. The
request's `room_id` becomes a prefix hint; the server chooses the
final ID. If the suffixed ID itself collides (astronomically
unlikely), regenerate; cap at 10 tries and raise 500 on exhaustion.

- Add `_resolve_available_room_id(session, prefix, max_tries=10)` in
  `rooms.py` that generates `{prefix}-{suffix}`, checks for
  collision, and retries until free. Suffix alphabet is
  `string.ascii_letters + string.digits` (consistent with the
  `^[a-zA-Z0-9\-_]+$` validator), length 6.
- Delete the `if existing is not None: return ...` early-return
  block entirely.
- `RoomCreateResponse.created: bool` is removed — every call now
  creates a new room; the field is always `True` and carries no
  information.

From the caller's perspective, `POST /v1/rooms` always returns 201
with a `room_id` that differs from the request. The network surface
is identical regardless of whether the requested prefix collides
with anything in the database. No existence leak, no `frame_count`
leak, no idempotency for replays (which no caller relied on — every
legitimate create flow generates a fresh UUID and navigates away).

**Frontend — consume the returned `room_id`; drop the socket auto-create.**

First, the socket handler. Remove the `createRoom` fallback block in
`frontend/src/hooks/socketHandlers/connectionHandlers.ts:202-259`
(the `?copy_from=` URL lookup, the `createRoom` call, the 409 catch,
and the nested retry). On any non-200 `room_join` reply, set
`initializationError` to
`{message: "Room not found or not accessible", details: \`HTTP ${status}\`}`
(the `details` slot is typed `string`, so we stringify the status
code for debuggability). The server-supplied `detail` string is
discarded — it carries the misleading "Room with id X not found"
text that would leak the spec's existence-hiding promise on the UX
surface. Drop the now-unused `createRoom` import at the top of the
file.

Second, every create-and-navigate call site must route to
`response.room_id` (the actual server-assigned ID), because the
returned ID is now *always* different from the request. This is
correctness-critical, not an edge-case audit.

- `components/DuplicateRoomDialog.tsx:59` — already uses
  `result.room_id`. No change.
- `pages/templateSelection.tsx:122` — currently navigates via state
  derived from `newRoomId` (the requested). Switch to the returned
  value.
- `panels/RoomsPanel.tsx:52` — same pattern; switch.
- `panels/roomsHeaderActions.tsx:21, 32, 50` — currently
  `navigate('/rooms/${id}')` using the locally-generated `id`.
  Switch to the returned `room_id`.
- `panels/FilesystemPanel.tsx:131-159` — `targetRoomId` is the
  requested value; replace with the returned `room_id` before the
  `leaveRoom` + navigate step, and reuse it for the subsequent
  `uploadTrajectory` calls so uploads land on the actual new room.

The request-side `room_id` remains a required field on `RoomCreate`
(callers still pass a UUID); it is now a prefix the server uses to
name the created resource. `RoomCreateResponse.created` is removed
from the TypeScript client alongside the backend field.

### Existence-hiding contract

The parent refactor spec chose 404 over 403 on `can_read` denials to
hide room existence. That choice stands. After this cleanup:

- `room_join` 404 produces a generic frontend error, never an
  auto-create side-channel.
- `POST /v1/rooms` returns the same response shape in every case:
  201 with a server-chosen `{prefix}-{suffix}` ID. A caller who
  posts `{room_id: "target-id"}` cannot distinguish the three
  underlying cases ("target-id" is free / exists-readable /
  exists-unreadable) — the response is structurally identical. No
  existence leak, no `frame_count` leak.

## Testing

- **Backend — response shape is invariant across collision states.**
  New test in `tests/zndraw/test_rooms.py`: parameterized over
  (no existing row, existing readable row, existing unreadable row
  owned by a second user). In every case, `POST /v1/rooms
  {room_id: "X"}` returns 201 with a `room_id` matching the pattern
  `^X-[A-Za-z0-9]{6}$`. Untouched rooms remain untouched.
- **Backend — double-post yields two distinct rooms.** Same user
  posts `{room_id: "X"}` twice; both calls succeed, both return
  201, and the two returned `room_id` values are different. (This
  pins the removal of the idempotent-return branch.)
- **Backend — socket_events construction.** Update
  `tests/zndraw/test_socket_commands.py` for the two
  `RoomJoinResponse` constructions: drop the `locked=...` field and
  the `resp.locked` assertion.
- **Backend — existing tests that assume verbatim `room_id`.** Sweep
  for `response.json()["room_id"] == ...` patterns in existing room
  tests (grep `tests/zndraw/test_rooms.py`,
  `tests/zndraw/test_socketio_rooms.py`, and related); update each
  to use `startswith(prefix + "-")` or capture the returned ID and
  use it for follow-up calls.
- **Frontend tests / Playwright smoke.** Any Playwright spec that
  relies on the URL after a create operation must read the ID from
  the post-create URL / response, not hard-code it. Create-and-
  navigate assertions need to accept the suffixed ID.
- **Type-check pass.** `bun run typecheck` and `uv run pyright` must
  be clean after the edits. The `superuserLock` removal cascades
  from slice → store → context; the `locked` removal cascades from
  the backend Pydantic model to the frontend interface;
  `RoomCreateResponse.created` removal cascades across the TS
  client and any reader.

## Risks

- **Hidden readers of `RoomJoinResponse.locked`.** Grep across the
  repo (including tests and the pyclient) confirms the only reader
  after Change 1 is the `connectionHandlers.ts` line being removed
  in the same patch.
- **Callers that hard-code post-create `room_id`.** The `POST
  /v1/rooms` response ID now always differs from the request.
  Anything that reuses the requested ID for follow-up navigation,
  uploads, or test assertions breaks. Mitigation: the audit in
  Change 3 covers all frontend sites; the test sweep entry in the
  Testing section covers backend tests. External pyclient / CLI
  callers must be updated to read `response.room_id` (acceptable
  per the unshipped-code compat rule).
- **URL aesthetics.** Room URLs grow from 36 to 43 chars
  (`<uuid>-<6>`). No functional impact.
