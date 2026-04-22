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

**Backend — `POST /v1/rooms` on collision.**

Current behavior (`src/zndraw/routes/rooms.py:360-369`) returns 200
with the existing room's `frame_count` whenever the requested
`room_id` already exists, *regardless of whether the caller can_read
the existing row*. A caller can therefore probe any `room_id` and
learn whether it exists (201 = new, 200 = taken). This bypasses the
`404-hides-existence` contract that `GET /v1/rooms/{id}` and
`room_join` already honor, and it also leaks `frame_count` for
private rooms.

New behavior:

- If the requested `room_id` exists and the caller `can_read` it →
  return existing (200, idempotent). No leak: the caller already sees
  the row via normal APIs.
- If the requested `room_id` exists and the caller cannot `can_read` →
  **generate a fresh `room_id` by appending `-<6 random url-safe
  chars>` to the requested value, create a new room at the suffixed
  ID, and return 201 with the actual new `room_id`.** If the suffixed
  ID also collides (astronomically unlikely), regenerate; cap at 10
  tries and raise 500 on exhaustion.
- If the requested `room_id` does not exist → create (unchanged, 201).

Helper: add `_resolve_available_room_id(session, desired, max_tries=10)`
in `rooms.py` that encapsulates the loop. The suffix alphabet is
`string.ascii_letters + string.digits` (same character class already
accepted by the `^[a-zA-Z0-9\-_]+$` validator), length 6.

From the caller's perspective, `POST /v1/rooms` now always succeeds
with 200 or 201 and the response `room_id` tells them the actual
resource location — which may differ from what they requested when
collision-with-unreadable occurs.

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

Second, audit the create-and-navigate call sites so they route to
`response.room_id` (the actual server-assigned ID) rather than the
value they passed in. Current state:

- `components/DuplicateRoomDialog.tsx:59` — already uses
  `result.room_id`. No change.
- `pages/templateSelection.tsx:~122` — verify.
- `panels/RoomsPanel.tsx:~52` — verify.
- `panels/roomsHeaderActions.tsx:21, 32, 50` — currently navigates
  with the locally-generated `id`. Change to navigate with the value
  returned from `createRoom`.
- `panels/FilesystemPanel.tsx:131-159` — `targetRoomId` is the
  requested value; switch to the returned `room_id` before the
  `leaveRoom` + navigate step.

Each site today generates a fresh `crypto.randomUUID()`, so the
suffixed branch is effectively dead code — the audit is a correctness
guarantee so that the one scenario where it fires (a
`DuplicateRoomDialog` user typing an ID that happens to collide with
an unreadable room) routes cleanly instead of stranding the UI on
the wrong URL.

### Existence-hiding contract

The parent refactor spec chose 404 over 403 on `can_read` denials to
hide room existence. That choice stands. After this cleanup:

- `room_join` 404 produces a generic frontend error, never an
  auto-create side-channel.
- `POST /v1/rooms` no longer distinguishes the two collision states
  at the network surface. A caller who probes an ID owned by an
  unreadable room gets a 201 with a different `room_id` (the
  suffixed one) — the information they learn is "my ID suggestion
  wasn't used verbatim," not "a private room exists here." That
  residual signal is acceptable: with fresh v4 UUIDs the branch is
  statistically unreachable, and the content of the private room
  remains completely opaque.

## Testing

- **Backend — existence leak.** New test in
  `tests/zndraw/test_rooms.py` (or a dedicated file): two users, user
  A creates a private room `X`, user B posts `POST /v1/rooms`
  `{room_id: "X"}`. Assert response is 201 with `room_id` *not equal
  to* `"X"` (suffixed) and `created=True`. Assert user A's room `X`
  is untouched.
- **Backend — idempotent create for authorized caller.** User A
  posts `POST /v1/rooms` `{room_id: "X"}` twice; second call returns
  200 with `room_id == "X"` and `created=False`.
- **Backend — socket_events construction.** Update
  `tests/zndraw/test_socket_commands.py` for the two
  `RoomJoinResponse` constructions: drop the `locked=...` field and
  the `resp.locked` assertion.
- **Type-check pass.** `bun run typecheck` and `uv run pyright` must
  be clean after the edits. The `superuserLock` removal cascades from
  slice → store → context; the `locked` removal cascades from the
  backend Pydantic model to the frontend interface.
- **Manual / Playwright smoke.** Navigating to a brand-new UUID URL
  should land on the initialization-error screen, not on a
  freshly-created room. Navigating to an existing accessible room
  should continue to load normally. Creating a room via the landing
  page, `roomsHeaderActions`, and `FilesystemPanel` should all
  navigate to the server-returned `room_id`.

## Risks

- **Hidden readers of `RoomJoinResponse.locked`.** Grep across the
  repo (including tests and the pyclient) confirms the only reader
  after Change 1 is the `connectionHandlers.ts` line being removed in
  the same patch.
- **Residual leak via suffixed `room_id`.** When
  `POST /v1/rooms` is called with an ID that collides with an
  unreadable room, the 201 response carries a suffixed ID — the
  caller learns their verbatim ID "was not used." With fresh v4
  UUIDs (all legitimate callers), this branch is statistically
  unreachable; only a caller who manually enters a colliding ID in
  `DuplicateRoomDialog` could trigger it. Content of the private
  room remains fully opaque. Accepted as a design tradeoff vs.
  mandating server-generated IDs across every caller.
