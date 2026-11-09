# Room namespace refactor — design

**Branch:** `feat/room-group-scope-refactor` (continuation)
**Supersedes:** `docs/superpowers/specs/2026-04-22-room-scope-cleanup-design.md`

## Context

The previous cleanup spec proposed a server-side random 6-char suffix on every `POST /v1/rooms` to close the existence-leak channel where a caller could distinguish "name free" (201) from "name taken-but-hidden" (200/200-with-frame_count). That approach honored the existence-hiding invariant but silently broke the Python client contract `ZnDraw(room="X")`: the client stores `self.room = X` and never updates it from the server response, so any auto-suffixed ID would be unreachable on reconnect.

The fundamental tension: in a flat namespace, any `POST` either succeeds (proving the name was free) or fails (proving it was taken). You cannot simultaneously preserve "user picks the name" and "no existence leak via POST" in a single namespace.

This spec resolves it by **namespacing room IDs by owner UUID**: `{owner_uuid}/{room_name}`. Cross-namespace POST is rejected at the prefix check before any DB lookup, so cross-namespace probing is impossible by construction. Within own namespace, idempotent reuse is fine (the caller can already see the room). The Python `ZnDraw(room="user-uuid/my-room")` flow works end-to-end.

This spec also folds in the two smaller cleanups from the prior spec that are still valid: dead `superuserLock` state in the frontend and the unused `RoomJoinResponse.locked` field.

## Invariants (non-negotiable, carried from parent spec)

1. **Existence hiding.** Cross-namespace probes (GET, room_join, POST) return responses indistinguishable from "absent" — no caller can determine that a room they cannot access exists.
2. **No `frame_count` leak.** 404 response bodies carry no payload.
3. **Network surface invariant.** For any caller probing a resource they don't own, all underlying states (absent / readable-by-others-only / unreadable) produce structurally identical responses.
4. **No side-channel via response shape comparison.** Verbatim-vs-composed IDs, status-code divergence, body-shape divergence all disallowed for cross-namespace probes.
5. **Access control unchanged.** `can_read` / `can_edit` / `can_manage` semantics from `src/zndraw/access.py` stay as-is.
6. **`room_join` 404 → generic frontend error.** No auto-create side-channel from URL paste.

## Non-goals

- Human-readable owner slugs. UUIDs are used as namespace prefixes throughout. A future spec may layer a slug subsystem on top; out of scope here.
- Old-link redirect table after transfer. Old composed addresses break on transfer, by design (KISS). A future spec may introduce an aliases table.
- Group management UI improvements. Separate work.
- Cross-namespace URL-paste-creates-room. Removed in this spec — it was the existence-leak channel. The **own-namespace** case (User A pastes a URL in their own namespace, room doesn't exist yet → frontend auto-creates) is preserved and is leak-free by construction.

## Architecture

### Core principle: storage UUID vs. surface address

One rule, applied everywhere:

- **Internal/storage:** `Room.id` (a surrogate UUID) is the identity. Used in FK references, frame-storage prefixes, Redis keys, pub/sub channel names. Never surfaces to humans.
- **Wire protocol & user surfaces:** the composed `{owner_id}/{room_name}` is the identity. Used in HTTP URLs, socket event payloads, share links, Python `ZnDraw(room=...)`, CLI args, response bodies.

A surface field is never the surrogate UUID; an internal key is never the composed string. If you see a surrogate UUID in a response or a composed address in a Redis key, it's a bug.

### Composed addresses

Rooms are addressed publicly as **two-segment strings**: `{owner_uuid}/{room_name}`.

- `owner_uuid` is the UUID of the owning user (`Room.owner_user_id`) or owning group (`Room.owner_group_id`). XOR constraint already enforces exactly one is non-null.
- `room_name` is user-chosen, validated against `^[a-zA-Z0-9\-_]+$`.
- The composed string is the public identity: HTTP URLs, share links, socket event payloads, Python `ZnDraw(room=...)`, CLI room args.

There is **one canonical form**: full path, always. No implicit current-user-namespace resolution (`ZnDraw(room="my-room")` is rejected as malformed). One form, like `git clone https://github.com/owner/repo`.

### Surrogate UUID identity (Option C, industry standard)

`Room.id` continues to be a server-generated UUID — the **immutable internal identity**. Foreign keys (`RoomShareLink.room_id`, frame storage prefixes, socket channel names, Redis keys) all use this surrogate. Transfer never touches them.

The composed address is **computed at API boundaries**, never stored as a column. A `Room.public_address` property returns `f"{owner_uuid}/{room_name}"`.

The mapping from composed address to surrogate UUID is enforced by a `UNIQUE INDEX(COALESCE(owner_user_id, owner_group_id), room_name)`.

### Cross-namespace probe is impossible

- `POST /v1/rooms` rejects any request where `owner_user_id` is not the caller and `owner_group_id` is not a group the caller belongs to — at the prefix-validation step, before any DB lookup. Identical 403 response regardless of what exists inside the foreign namespace.
- `GET /v1/rooms/{X}/{Y}` returns 404 with no body whether `X` is an unknown UUID, `X` is known but `Y` doesn't exist, or the row exists but caller can't `can_read`.
- `room_join` socket event same as GET.

### Public rooms still live in their owner's namespace

`PUBLIC` visibility doesn't change the address. A public user-owned room is `{owner_user_uuid}/{name}`; a public group-owned room is `{owner_group_uuid}/{name}`. Anyone can `can_read`, but the URL still encodes ownership. No separate `/public/` root.

## Data model

### `src/zndraw/models.py`

Required new imports: `Index`, `func` from sqlalchemy. (`UniqueConstraint` already imported; `Index` and `func` are not.)

Room model changes:

| Column | Type | Change |
|---|---|---|
| `id` | UUID PK (as `str`) | Unchanged — surrogate, immutable. |
| `owner_user_id` | UUID, FK→user, nullable | Unchanged. |
| `owner_group_id` | UUID, FK→group, nullable | Unchanged. |
| `room_name` | `str` | **New.** Validated by `^[a-zA-Z0-9\-_]+$`. |
| `visibility` | enum | Unchanged. |
| (other fields) | … | Unchanged (description, step, frame_selection, default_camera, created_at, created_by_id). |

New constraint added to `__table_args__`:

```python
Index(
    "ux_room_owner_name",
    func.coalesce(owner_user_id, owner_group_id),
    room_name,
    unique=True,
)
```

Existing `room_owner_exactly_one` and `room_visibility_matches_owner` CHECKs are kept verbatim.

New method on `Room`:

```python
@property
def public_address(self) -> str:
    owner = self.owner_user_id or self.owner_group_id
    return f"{owner}/{self.room_name}"
```

The old `Room.id` field-level regex validator (single-segment `^[a-zA-Z0-9\-_]+$`) is **removed** — `Room.id` is server-generated UUID, never user-supplied. The regex moves to `Room.room_name`.

### Composed-address helpers

**Backend: no parser needed.** Every wire surface carries the split form — URL routes use two path params (`{owner_id}/{room_name}`), POST bodies and all socket events use split fields. FastAPI parses path params; Pydantic validates each field independently. There is no central `parse_room_address` helper. The composed string only appears in *outbound* response fields, built by `Room.public_address`.

**Frontend `frontend/src/utils/roomAddress.ts`:** only a composer is needed.

```ts
export const composeRoomAddress = (ownerId: string, name: string) => `${ownerId}/${name}`;
```

React Router's `useParams<{ownerId, roomName}>()` already gives the split form; socket emits use split fields. The composer is only used for share-link URLs and other display contexts.

**Python client: parses once, inline.** `src/zndraw/client/core.py` accepts `ZnDraw(room="<owner_uuid>/<name>")` and splits at `__init__` (single `str.split("/", 1)` call with a validation check). No helper module needed; the parse is one line at one site.

## API contract

### URL routes (two path segments)

All room routes update from `/{room_id}` to `/{owner_id}/{room_name}`. Every consumer of `room_id: str = Path()` becomes `(owner_id: UUID = Path(), room_name: str = Path())`.

**Routes in `src/zndraw/routes/rooms.py`:** GET, PATCH, DELETE on the new two-segment path.

**Content sub-routes (15 files):**

- `src/zndraw/routes/bookmarks.py:29` — `/v1/rooms/{owner_id}/{room_name}/bookmarks`
- `src/zndraw/routes/chat.py:36` — `/v1/rooms/{owner_id}/{room_name}/chat/messages`
- `src/zndraw/routes/edit_lock.py:35` — `/v1/rooms/{owner_id}/{room_name}/edit-lock`
- `src/zndraw/routes/figures.py:31` — `/v1/rooms/{owner_id}/{room_name}/figures`
- `src/zndraw/routes/frames.py:64` — `/v1/rooms/{owner_id}/{room_name}/frames`
- `src/zndraw/routes/geometries.py:87` — `/v1/rooms/{owner_id}/{room_name}/geometries`
- `src/zndraw/routes/geometries.py:369` — secondary router; same path prefix
- `src/zndraw/routes/isosurface.py:32` — `/v1/rooms/{owner_id}/{room_name}/frames/{index}/isosurface`
- `src/zndraw/routes/presets.py:40` — `/v1/rooms/{owner_id}/{room_name}/presets`
- `src/zndraw/routes/progress.py:26` — `/v1/rooms/{owner_id}/{room_name}/progress`
- `src/zndraw/routes/screenshots.py:43` — `/v1/rooms/{owner_id}/{room_name}/screenshots`
- `src/zndraw/routes/selection_groups.py:31` — `/v1/rooms/{owner_id}/{room_name}/selection-groups`
- `src/zndraw/routes/share_links.py:28` — `/v1/rooms/{owner_id}/{room_name}/share-links`
- `src/zndraw/routes/step.py:23` — `/v1/rooms/{owner_id}/{room_name}/step`
- `src/zndraw/routes/trajectory.py:52` — `/v1/rooms/{owner_id}/{room_name}/trajectory`

A missed router silently 404s every request to that resource. Implementation plan enumerates them as a checklist.

### Access dependencies (`src/zndraw/dependencies.py:609-611`)

Replace `AccessReadDep` / `AccessEditDep` / `AccessManageDep` (and their underlying `get_readable_room` / `get_editable_room` / `get_manageable_room` at lines 546-606) with two-path-segment variants that look up `Room` by `(owner_id, room_name)` via the UNIQUE index. Same `AccessContext` (NamedTuple at lines 499-515) emitted; just sourced from two path params instead of one.

Failure paths all return 404 — never 403 — for read/list paths to preserve existence hiding.

### `POST /v1/rooms`

The existing `RoomCreate` Pydantic model (`src/zndraw/schemas.py:50-57`) currently has `{room_id: str, owner_group_id: UUID | None, ...}`. It is restructured to:

```json
{
  "owner_id": "8f3a-...",
  "name": "my-experiment",
  "description": "..." | null,
  "copy_from": "..." | null
}
```

A **single `owner_id` field** identifies the target namespace. The server resolves whether `owner_id` refers to a user or a group by lookup, then writes to the appropriate column. This follows the standard polymorphic-ownership pattern used by GitLab (`namespace_id`), HuggingFace (`{namespace}/{name}`), and similar systems — KISS for callers; the boolean type-tag (`is_group`) is redundant when the UUID is already globally unique across user/group tables.

Server validation (in order):

1. `name` matches `^[a-zA-Z0-9\-_]+$` (Pydantic field validator).
2. **Permission gate, before any room-table lookup:**
   - If `owner_id == current_user.id` (or caller is superuser) → user-owned create permitted; skip to step 4 with the user column targeted.
   - Otherwise, look up `owner_id` in the `group` table.
     - Found AND caller has ≥ MEMBER role → group-owned create permitted; skip to step 4 with the group column targeted.
     - Found but caller is not a member → 403.
     - Not found in `group` table (i.e. `owner_id` is neither caller's own UUID nor a group) → 403, identical response.
3. (Superuser only) Look up `owner_id` in `user` table when it didn't match a group; if it's a real user that isn't the caller, create on their behalf with the user column targeted.
4. Insert or fetch existing row by `(COALESCE(owner_user_id, owner_group_id), name)` via the UNIQUE index:
   - No match → INSERT → 201 `created: true`.
   - Match → return existing row → 200 `created: false` (idempotent reuse; possible only within caller's own namespace because step 2 already gated foreign-namespace POSTs).

All "you cannot create here" cases (step 2) return identical 403 bodies regardless of whether the foreign namespace exists, is a user, is a group, or contains a room with the requested name.

Two PK lookups per POST is acceptable (`user` and `group` UUID lookups are O(1)). The flow short-circuits in the common case (`owner_id == current_user.id`) with zero auxiliary lookups.

Response shape:

```json
{
  "room_id": "8f3a.../my-experiment",
  "created": true,
  "frame_count": 0
}
```

The early-return branch currently at `src/zndraw/routes/rooms.py:361-369` and its preceding format validator at lines 351-358 are removed; the unique-index + idempotent-reuse logic replaces both.

### `PATCH /v1/rooms/{owner_id}/{room_name}` — transfer (C1)

Body uses the same single-`owner_id` shape as POST. The transfer field is `new_owner_id: UUID` (replaces the two XOR fields the current PATCH uses at `rooms.py:633-665`):

```json
{
  "new_owner_id": "..." | null,
  "description": "..." | null,
  "visibility": "..." | null,
  "frame_count": 0 | null
}
```

Server logic when `new_owner_id` is set:
- Resolve kind via lookup (user table first if it equals caller, else group table — same pattern as POST step 2).
- Validate caller can transfer to that destination (group membership for group transfers; superuser for user-to-user transfers — preserve current rules).
- **Destination collision check.** If a row already exists with `(COALESCE(new_owner_user_id, new_owner_group_id), room.room_name)`, return 409 `TransferTargetInvalid`. No mutation.
- **Mutate owner columns only.** `Room.id` (surrogate UUID) and all FK references (`RoomShareLink.room_id`, frame storage prefix, `room_channel(room.id)`, all Redis keys) are untouched. The composed `public_address` changes because the underlying owner columns changed.
- **Emit `room_renamed` socket event** on the existing channel (`room_channel(room.id)` — the surrogate UUID, unchanged) with `{old_address, new_address}` so connected clients can navigate to the new URL.
- Response includes the new composed `room_id`.

Old composed addresses break. The implementation plan does not include an alias/redirect mechanism.

### Socket layer

**All four room-scoped event models in `src/zndraw/socket_events.py` split `room_id: str` into `(owner_id: UUID, room_name: str)`:**

```python
class RoomJoin(BaseModel):
    owner_id: UUID
    room_name: str
    client_type: Literal["frontend", "pyclient"] = "frontend"


class RoomLeave(BaseModel):  # was lines 28-31
    owner_id: UUID
    room_name: str


class TypingStart(BaseModel):  # was lines 38-41
    owner_id: UUID
    room_name: str


class TypingStop(BaseModel):  # was lines 44-47
    owner_id: UUID
    room_name: str
```

Server handlers for each look up the row by `(owner_id, room_name)` via the UNIQUE index, then operate on the surrogate UUID internally (e.g. `room_channel(room.id)`). Per-field Pydantic validation gives clear errors for malformed UUIDs vs malformed names.

**Outbound `RoomJoinResponse`** at `src/zndraw/socket_events.py:55-66`:
- Remove `locked: bool` (carried from parent spec).
- `room_id: str` field carries the **composed address** (`Room.public_address`), not the surrogate UUID. Surface = composed, per the core principle.
- Both construction sites in `src/zndraw/socketio.py` (system-room path at lines 244-251, normal path at lines 327-336) update accordingly.

`room_channel(room_id: str)` at `src/zndraw/dependencies.py:188-190` continues to use the surrogate UUID internally — `f"room:{room_id}"` where `room_id` is the `Room.id` column. Pub/sub channel keys are stable across transfers. The composed address never appears in Redis or pub/sub key space.

The `room_join` handler at `src/zndraw/socketio.py:178-184` resolves `(owner_id, room_name)` via the new dependency factory; 404-equivalent error response on any lookup failure or access denial.

### `room_renamed` socket event (new)

Emitted on the room's surrogate-UUID channel after a successful PATCH that changes the owner. Payload: `{old_address: str, new_address: str, room_id: str /* surrogate UUID, for client correlation */}`.

Frontend handler:
1. If `old_address` matches the currently-open room → snackbar + navigate to `/rooms/{new_address}`.
2. Let `room_join` re-run against the new address (succeeds if caller still has access, 404 otherwise).

### Response envelopes — owner display labels (frontend-dumb path)

Every room response that flows to the frontend gains:

- `room_id: str` — the composed `{owner_id}/{room_name}` (the public address).
- `owner_id: UUID` — the bare owner UUID (convenience for permission equality checks on the frontend).
- `owner_kind: Literal["user", "group"]` — determined by which owner column is non-null.
- `owner_label: str` — user's email (for user owners) or group's `name` (for group owners). Backend resolves once at serialization time. Avoids N+1 fetches on the frontend.

`RoomResponse` at `src/zndraw/schemas.py:60-72` is restructured: drop the separate `owner_user_id` / `owner_group_id` fields in favor of `owner_id` + `owner_kind` (no info loss; the frontend never needs to distinguish them as separate UUID slots). Other room response models (`RoomCreateResponse`, the room-detail response if separate, room list items) gain the same fields.

### Existence-hiding response matrix (final)

| Caller probe | Underlying state | Response |
|---|---|---|
| `GET /v1/rooms/X/Y` | owner X doesn't exist | 404, no body |
| `GET /v1/rooms/X/Y` | owner X exists, no room Y | 404, no body |
| `GET /v1/rooms/X/Y` | room exists, caller can't read | 404, no body |
| `GET /v1/rooms/X/Y` | room exists, caller can read | 200 |
| `POST /v1/rooms {owner=foreign, ...}` | (any state inside foreign namespace) | 403, identical body |
| `POST /v1/rooms {owner=self, name=Y}` | Y doesn't exist | 201, `created: true` |
| `POST /v1/rooms {owner=self, name=Y}` | Y exists in self's namespace | 200, `created: false` |
| `room_join {owner_id, room_name}` | any failure | 404-equivalent error |
| URL paste `/rooms/<self_uuid>/Y` | (any state, frontend retries via own-namespace POST) | 200/201 from POST — never reveals more than the matching POST row above |
| URL paste `/rooms/<foreign>/Y` | (any state) | generic frontend error; no POST attempt |

All cross-namespace responses are byte-identical. Within-own-namespace responses leak only the caller's own state (which they can already see).

## Frontend impact

### React Router (`frontend/src/App.tsx:38-63`)

Three route definitions updated:

- L52: `/rooms/:roomId/files` → `/rooms/:ownerId/:roomName/files`
- L56: `/rooms/:roomId` → `/rooms/:ownerId/:roomName`
- L60: `/room/:roomId` legacy alias — **removed** per no-backwards-compat policy.

### `useParams` call sites

Update all to `const { ownerId, roomName } = useParams<{ ownerId: string; roomName: string }>();` and convert to composed address at API/socket boundaries via `composeRoomAddress(ownerId, roomName)`.

- `frontend/src/App.tsx:17` (`FilesystemRedirect`)
- `frontend/src/pages/landingPage.tsx:64` (`MainPage`)
- `frontend/src/components/RoomManagementMenu.tsx:50`
- `frontend/src/components/ConnectionDialog.tsx:21`
- `frontend/src/panels/ChatPanel.tsx:184`

### URL construction sites

All updated to use the composed address from the server response, never a client-side constructed value:

| File | Line | Change |
|---|---|---|
| `App.tsx` | 19 | Use composed address from `useParams`. |
| `panels/RoomsPanel.tsx` | 52, 60, 93 | Capture response `room_id`, navigate via that. |
| `panels/roomsHeaderActions.tsx` | 21, 32, 50 | Capture response `room_id` (3 sites). |
| `panels/FilesystemPanel.tsx` | 131, 160 | Capture response `room_id`, reuse for subsequent uploads. |
| `components/RoomManagementMenu.tsx` | 207 | Use composed address. |
| `components/DuplicateRoomDialog.tsx` | 59 | Already uses `result.room_id` ✓ — no change. |
| `pages/templateSelection.tsx` | 104, 123, 126 | Capture response `room_id`. |
| `panels/roomRowMenu.tsx` | 57 | `${window.location.origin}/rooms/${room.id}` — fine, `room.id` will now be the composed string. |

### Client-side UUID generation — removed

Every `crypto.randomUUID()` call producing a room ID is removed:

- `components/DuplicateRoomDialog.tsx:51`
- `pages/templateSelection.tsx:113`
- `panels/RoomsPanel.tsx:50`
- `panels/roomsHeaderActions.tsx:18, 29, 43`

Create dialogs take a name string input from the user; POST body becomes `{owner_user_id: currentUser.id, name, ...}` (or `{owner_group_id, name, ...}` in group context). The server composes the address and returns it.

### Owner label rendering

The room-list UI uses the new `owner_label` + `owner_kind` fields from response payloads. No frontend lookup logic. `RoomManagementMenu.tsx:175-177`'s owner-equality check stays as-is (it compares UUIDs, not labels).

### `connectionHandlers.ts` rewrites

- L24-31 `RoomJoinResponse` interface: remove `locked: boolean`.
- L69 `superuserLock` write: removed.
- L204-246 `createRoom` fallback block: **rewritten, not deleted.** New logic gates on namespace ownership:
  - If `ownerId === currentUser.id` (URL points at the caller's own user namespace): POST `/v1/rooms {owner_user_id: currentUser.id, name: roomName, copy_from?: <from URL param>}`. On success, retry `room_join` at the new room. The `?copy_from=` URL param (L205-207) is **kept** for this flow — it remains a useful URL-template entry point. The 409 catch-and-retry stays (handles the idempotent-reuse race).
  - Otherwise (cross-namespace, including group URLs the caller may not be a member of): set `initializationError = {message: "Room not found or not accessible", details: \`HTTP ${status}\`}` and stop. No POST attempt. The server's `detail` string is discarded.
  - The `createRoom` import at the top of the file is **kept** (still used by the gated path).
- `RoomJoin` socket emit: send `{owner_id, room_name, client_type}` (split form).

Note: scoping auto-create to the caller's own *user* namespace (not their groups) is intentional KISS — group rooms are created via the explicit group-room UI, not URL paste. A user who wants a group room can use the rooms panel from the group view. This avoids the frontend needing to enumerate the caller's group memberships before deciding whether to POST.

### Store: `selectIsRoomReadOnly` (`frontend/src/store.tsx:59-68`)

Remove the `superuserLock` branch at L62-66. Selector reduces to the `userLock` check only.

### `superuserLock` removal (full sweep, carried from parent spec)

- `frontend/src/stores/slices/lockSlice.ts`: drop `superuserLock` field and `setSuperuserLock` action.
- `frontend/src/store.tsx`: `selectIsRoomReadOnly` simplification (above).
- `frontend/src/components/geometry/GeometryGrid.tsx`: drop `superuserLock` selector, the corresponding branch in `canEdit` (L77), and the "Room is locked" tooltip fallback (L141).
- `frontend/src/hooks/useSocketManager.ts`: drop `setSuperuserLock` capture and pass-through (L46, 93, 227).
- `frontend/src/hooks/socketHandlers/types.ts`: drop `setSuperuserLock` from `HandlerContext`.
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts:69`: drop the dead write.

## Python client (`src/zndraw/client/`)

### Constructor parses once

`src/zndraw/client/core.py:193`: `self.room` becomes the full composed address (validated against `^[0-9a-f-]+/[a-zA-Z0-9\-_]+$`). Parsing into `(owner_id, room_name)` happens once at `__init__`. Auto-UUID for the no-arg case becomes auto-generation of `{current_user.id}/{uuid4()}`.

`ZnDraw(room="my-room")` (single-segment) is rejected with a clear error pointing at the new format.

### `APIManager.create_room` (`api.py:132-159`)

Payload changes from `{room_id, copy_from?, description?}` to:

```python
{
    "owner_user_id": <parsed owner UUID> | None,
    "owner_group_id": None | <parsed owner UUID>,
    "name": <parsed name>,
    "copy_from": ...,
    "description": ...,
}
```

`copy_from` itself remains a composed address (the source room).

### `core.py:231-239` GET-then-POST flow

Stays as today's structure (GET to detect existence, POST to create if missing) but with two-segment URLs and the split POST body. The Python flow remains `ZnDraw(room="user-uuid/foo")` → reconnects to existing room, or creates it under the caller's own namespace if free.

### Server reconnect on transfer

The pyclient subscribes to `room_renamed` events. On rename, updates `self.room` to the new composed address and continues. (Symmetric to the frontend handler.)

## Other in-repo callers (CLI, joblib, executor, MCP, docs)

Audit findings — every site that constructs room URLs or instantiates `ZnDraw` with a room argument needs updating in lockstep with the core refactor.

### CLI (`src/zndraw/cli.py` and `src/zndraw/cli_agent/`)

| File | Line | Site | Change |
|---|---|---|---|
| `src/zndraw/cli.py` | 184 | `upload_file()` passes a positional room string to `ZnDraw(...)` | Accept the composed `<owner_uuid>/<name>` form; reject single-segment input. |
| `src/zndraw/cli_agent/rooms.py` | 47 | `list_rooms()` calls `GET /v1/rooms` | Response shape gains `owner_id`/`owner_kind`/`owner_label`/`room_id` (composed). Update display. |
| `src/zndraw/cli_agent/rooms.py` | 71 | `create_room()` POSTs `RoomCreate.model_dump()` | Build `{owner_id, name, ...}` body. Default `owner_id` to the authenticated caller's UUID. |
| `src/zndraw/cli_agent/rooms.py` | 118 | `open_room()` builds `room_url = f"{settings.url}/rooms/{room}"` | `{room}` is now a composed string; URL becomes `/rooms/{owner_id}/{room_name}`. |
| `src/zndraw/cli_agent/connection.py` | 214 | `get_current_step()` calls `/v1/rooms/{room}/step` | Pass the composed address as a two-segment path. |
| `src/zndraw/cli_agent/connection.py` | 309 | `get_zndraw()` passes `room=` to `ZnDraw(...)` | `room` must be composed form by this point. |
| `src/zndraw/cli_agent/mount.py` | 49 | `url: f"{conn.base_url}/room/{vis.room}"` | Two issues: legacy `/room/` (single) is removed; URL becomes `/rooms/{vis.room}` where `vis.room` is composed. |

### joblib (`src/zndraw_joblib/client.py`)

| Line | Site | Change |
|---|---|---|
| 102-108 | `ClaimedTask.__init__(room_id: str)` | The `room_id` field carries the composed address from the API layer. Type stays `str`; semantics change. |
| 424 | `GET /v1/joblib/rooms/{room_id}/jobs` | Use composed two-segment form (route on the backend updates correspondingly). |
| 575 | `POST /v1/joblib/rooms/{room}/tasks/{job_name}` | Same. |
| 806 | `POST /v1/rooms/{room_id}/chat/messages` | Same. |

Backend joblib routes under `/v1/joblib/rooms/` also need their path params split — flag for the implementation plan to enumerate alongside the 15 content sub-routes.

### Executor (`src/zndraw/executor.py:63-67`)

`_run()` instantiates `ZnDraw(..., room=room_id)`. `room_id` arriving here originates from a `ClaimedTask` (joblib). With the joblib chain updated, the value is already composed; no logic change here, just verify the type carries through.

### CLI scripts integration test surface

E2E specs use CLI invocations like `rooms create --room ${ROOM}` (`frontend/e2e/socket-sync.spec.ts`). The CLI's `--room` arg now requires the composed form. Rewriting specs to capture the composed address from a real `rooms create` response (rather than hardcoding) is the cleaner approach.

### Documentation

| File | Lines | Current | Change |
|---|---|---|---|
| `README.md` | 60, 70 | `ZnDraw(url="...", room="my-room")` | Use a composed example, e.g. `ZnDraw(url="...", room=f"{user_uuid}/my-room")`. Add a note that `room=` accepts the composed form returned by the server. |
| `docs/source/python-api.rst` | 22 | Same | Same. |
| `docs/source/python-api.rst` | 1033, 1040 | `vis.register_job(..., room="my-room-id")` / `room="@global"` | Update examples; verify whether `@global` magic still applies (likely retained as a sigil for the system room; flag for implementation). |

### Tests (Python)

The audit identified ~144 `ZnDraw(...)` instantiations across 14 test files. Patterns to rewrite:

- `ZnDraw(url=server, room="literal-string")` → construct the room via the model layer (or via the CLI/API), capture the composed address, use that.
- `ZnDraw(url=server, room=uuid.uuid4().hex)` → no longer a valid room ID by itself; prepend the test user's UUID: `room=f"{test_user.id}/{uuid.uuid4().hex}"`.

The test fixtures should expose a helper (e.g., `make_room_address(user_id, name)`) — `tests/zndraw/conftest.py` is the natural home.

## Removals summary

Backend:
- `src/zndraw/routes/rooms.py:351-358` — old `room_id` user-input format validator.
- `src/zndraw/routes/rooms.py:361-369` — "return existing" 200 early-return branch.
- `Room.id` field-level regex validator (single-segment).
- `RoomJoinResponse.locked: bool` field + both construction sites in `src/zndraw/socketio.py` (lines 244-251, 327-336).
- Single-segment `RoomDep`-shaped path-param dependencies (functions at `dependencies.py:546-606`, dep types at lines 609-611) — replaced by two-segment variants.
- `src/zndraw/schemas.py:50-57` `RoomCreate`: the `{room_id: str, owner_group_id: UUID | None}` shape replaced by `{owner_id: UUID, name: str, ...}`.
- `src/zndraw/schemas.py:60-72` `RoomResponse`: separate `owner_user_id` / `owner_group_id` fields replaced by `owner_id` + `owner_kind` + `owner_label`; gains composed `room_id`.
- `src/zndraw/socket_events.py:28-31, 38-41, 44-47` — single `room_id: str` on `RoomLeave`, `TypingStart`, `TypingStop` replaced by split `(owner_id, room_name)`.

Frontend:
- `frontend/src/App.tsx:60` — `/room/:roomId` legacy alias.
- `connectionHandlers.ts:69` — `setSuperuserLock` write.
- `connectionHandlers.ts:204-246` `createRoom` fallback block: **rewritten** (not removed) — gated on own-user-namespace; cross-namespace 404 → generic error. See the connectionHandlers.ts rewrites section above.
- `connectionHandlers.ts:24-31` — `locked: boolean` from `RoomJoinResponse` interface.
- All `crypto.randomUUID()` client-side room ID generation: `DuplicateRoomDialog.tsx:51`, `templateSelection.tsx:113`, `RoomsPanel.tsx:50`, `roomsHeaderActions.tsx:18, 29, 43`.
- `superuserLock` field/action/selector/branch across `lockSlice.ts`, `store.tsx`, `GeometryGrid.tsx`, `useSocketManager.ts`, `socketHandlers/types.ts`.

Tests:
- `tests/zndraw/test_socket_commands.py`: `locked=…` in `RoomJoinResponse` constructions, `resp.locked` assertions, single-segment `RoomJoin` constructions.
- 10 backend test files with hardcoded single-segment IDs (`test_access_deps.py`, `test_socket_commands.py`, `test_share_token_dep.py`, `test_redis.py`, `test_socketio_rooms.py`, `test_access_matrix.py`, `test_gif.py`, `test_routes_bookmarks.py`, `test_routes_figures.py`, `test_routes_geometries.py`) — rewritten to construct rooms via the model layer and use `Room.public_address`.
- 5 Playwright specs (`registration.spec.ts`, `socket-sync.spec.ts`, `chat-features.spec.ts`, `constraint-visualization.spec.ts`, `ui-panels-chat.spec.ts`) with hardcoded `ROOM = "test-..."` constants — rewritten to capture the response address from a real create flow.

## Migration

All databases are purged on cutover (project policy: no backwards compat, no shipped releases). No migration code. The next `create_db_and_tables()` produces the new schema; first room created lives at `{owner_uuid}/{name}`.

External pyclient/CLI users on a previous build must update their `ZnDraw(room=...)` arguments and any URL bookmarks. Acceptable per the unshipped-code policy.

## Test plan

### Existence-leak invariants (centerpiece)

**Parametric cross-namespace POST test** in `tests/zndraw/test_rooms.py`:

States to parameterize:
1. Name free in caller's own namespace → 201, `created: true`.
2. Name taken in caller's own namespace → 200, `created: false`, references same surrogate UUID.
3. POST with `owner_user_id` set to another user (caller is not superuser) → 403.
4. POST with `owner_user_id` set to another user who happens to have a room with that name → 403.

Assert: responses (3) and (4) are byte-identical.

**Parametric GET probe test:**

`GET /v1/rooms/{other_user_id}/X` across:
- owner UUID is unknown
- owner exists, no room X
- room exists, caller cannot read
- (positive control) room exists, caller can read

First three return 404 with empty body, byte-identical. Fourth returns 200.

**Idempotent reuse:**

Same user POSTs same `(owner_user_id, name)` twice. Both succeed; second returns `created: false`; both responses reference the same surrogate UUID.

**Group POST gating:**

Non-member POSTs `{owner_group_id: G, name: "Y"}` → 403, identical shape to user-namespace-foreign 403.

### Transfer semantics

- Transfer collision: PATCH `{owner_user_id: B}` where `(B, room_name)` already exists → 409 `TransferTargetInvalid`, no mutation.
- Transfer happy path: PATCH succeeds; `Room.id` (surrogate UUID) unchanged; `RoomShareLink.room_id` and frame storage prefix unchanged; response carries new composed `room_id`.
- `room_renamed` event emission on transfer: payload includes correct `old_address` and `new_address`.

### Socket layer

- `tests/zndraw/test_socket_commands.py` updated: split `RoomJoin` constructions, no `locked=…`, no `resp.locked` assertions.
- New test: `room_join` against a cross-namespace room returns 404-equivalent error, no body distinguishing absent/hidden.

### Backend test data sweep

Rewrite the 10 listed backend test files to construct rooms via the model layer (`Room(owner_user_id=test_user.id, room_name="r1", ...)`) and use `Room.public_address` for URL construction.

### Frontend E2E

Rewrite the 5 listed Playwright specs to capture the response address from a real create flow (no hardcoded room name constants).

### Static checks

- `uv run pyright` clean.
- `bun run typecheck` clean.

## Risks

- **15 sub-route surface.** Every content router needs path schema and access-dep update in lockstep. A missed router silently 404s. Implementation plan enumerates them.
- **Public URL vs internal channel divergence.** Public URLs and POST/PATCH paths use composed `{owner_id}/{room_name}`; pub/sub channels and Redis keys use surrogate UUID `Room.id`. Intentional (transfer doesn't churn pub/sub). Reviewers may assume they're the same — called out here and in implementation-plan section headers.
- **Owner display label sourcing.** Backend serializers gain `owner_label` (email for users, name for groups). Need to confirm no PII concern; emails are visible to anyone granted access today, so no new exposure. Document in implementation plan.
- **Python `ZnDraw(room="X")` breaking change.** Single-segment IDs now error out. Acceptable per unshipped-code policy; clearly documented in error message.
- **5 E2E specs break on day one.** Implementation plan must update them in the same PR.
- **CLI room args.** `rooms create --room ${ROOM}` and similar in tests need the new two-segment form. Implementation plan must audit the CLI surface — concrete sites enumerated in the "Other in-repo callers" section.
- **joblib routes (`/v1/joblib/rooms/{room_id}/...`).** Backend joblib routes parallel the content sub-routes and need the same two-segment path-param treatment. Locations: under `src/zndraw_joblib/`. Implementation plan should enumerate them as a separate checklist (the 15 content sub-routes do *not* include these).
- **Polymorphic owner_id resolution requires user/group UUID uniqueness.** Both tables use UUID v4 PKs; collision probability is astronomical (2^122 namespace, single shared global UUID space). No additional disambiguator is needed. If the design ever moves to non-UUID identifiers per table, the polymorphic `owner_id` field must be revisited.
- **URL length.** Public URLs grow from ~36 chars to ~80 chars. Cosmetic; accepted as the price of KISS (no slug subsystem).
