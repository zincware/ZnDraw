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

### Change 3 — disable URL-based room auto-creation

Frontend-only change in `frontend/src/hooks/socketHandlers/connectionHandlers.ts`.

Current behavior (lines 202-259): on a `room_join` reply with
`status === 404`, the handler reads `?copy_from=` from the URL, calls
`createRoom({room_id, copy_from})`, catches 409 to keep going on a
race, and retries `room_join`. That success path hides both "room did
not exist" (now does) and "caller cannot access the existing room"
(still 404s on retry, surfacing as the misleading snackbar).

New behavior:

- On any non-200 `room_join` reply, set `initializationError` to a
  single generic value:
  - `message`: `"Room not found or not accessible"`
  - `details`: the numeric `status` (for debuggability). The
    server-supplied `detail` string is discarded — it carries the
    misleading "Room with id X not found" message that violates the
    existence-hiding contract on the UX surface.
- The handler does not call `createRoom`, does not read `copy_from`
  from the URL, and does not retry. The fallback block and its
  `?copy_from=` URL lookup are deleted together.
- The `createRoom` import at the top of `connectionHandlers.ts`
  becomes unused and is removed.

Room creation remains available through the intentional paths, all of
which still work unchanged:

- Landing page (`pages/templateSelection.tsx`, line 122).
- `RoomsPanel` and `roomsHeaderActions`.
- `DuplicateRoomDialog`.
- `FilesystemPanel`.
- CLI `zndraw file.xyz`, which hits `POST /v1/rooms` directly and then
  redirects the browser to the new URL.

The initialization-error surface already exists and renders via the
main layout; no new component is needed. A "Go back" / "Create new
room" link to `/` is already present on the error screen.

### Existence-hiding contract

The parent refactor spec chose 404 over 403 on `can_read` denials to
hide room existence. That choice stands. With the auto-create fallback
gone, the frontend can no longer betray existence through a
"create-then-retry" side channel — a caller who hits a room they
cannot see gets the same generic "not found or not accessible"
message as a caller who mistyped a URL.

## Testing

- **Backend.** Update `tests/zndraw/test_socket_commands.py` for the
  two `RoomJoinResponse` constructions. No new backend test is needed;
  the existing access-matrix and share-link suites already cover the
  404-vs-200 semantics on `room_join`.
- **Type-check pass.** `bun run typecheck` and `uv run pyright` must
  be clean after the edits. The `superuserLock` removal cascades from
  slice → store → context; the `locked` removal cascades from the
  backend Pydantic model to the frontend interface.
- **Manual / Playwright smoke.** Navigating to a brand-new UUID URL
  should land on the initialization-error screen, not on a
  freshly-created room. Navigating to an existing accessible room
  should continue to load normally.

## Risks

- **User muscle memory on URL-based create.** Users who relied on
  pasting a UUID to spin up a fresh room lose the shortcut. The error
  screen's "Create new room" link to `/` covers the remaining flow.
- **Hidden readers of `RoomJoinResponse.locked`.** Grep across the
  repo (including tests and the pyclient) confirms the only reader
  after Change 1 is the `connectionHandlers.ts` line being removed in
  the same patch.
- **Stale `?copy_from=` URL params.** If any documentation or bookmark
  relied on entering `/rooms/new-id?copy_from=template` to
  auto-create a copy, it will now fail. The landing page and
  `DuplicateRoomDialog` still support `copy_from` explicitly.
