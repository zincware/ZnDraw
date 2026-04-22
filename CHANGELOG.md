# Changelog

## [Unreleased]

### BREAKING: Room scope refactor (2026-04)

**Database recreation required.** No Alembic migration is provided — drop
the old SQLite/Postgres schema before upgrading. Existing rooms,
memberships, and ACLs are **not preserved**.

#### Data model

- Replaced `Room.is_public` and `Room.locked` with:
  - `Room.visibility`: `private` | `group` | `public` (three-value scope)
  - `Room.owner_user_id` XOR `Room.owner_group_id` (polymorphic owner)
  - Two `CHECK` constraints enforce exactly-one-owner and
    visibility/owner consistency at the DB level.
- Removed `RoomMembership` and `MemberRole` entirely. Per-room ACLs are
  replaced by per-group membership.
- Added `Group` (UUID id, globally-unique name, creator), `GroupMembership`
  (group ↔ user, role `viewer` | `member` | `admin`), and `RoomShareLink`
  (bearer-token capability scoped to one room, `view` or `edit`, optional
  expiry + revocation).
- Added `User.is_guest: bool = False`. Set by `POST /v1/auth/guest`. The
  old `@guest.user` email-suffix identification is gone.

#### Access control

- Single source of truth: `src/zndraw/access.py` exposes pure predicates
  `can_read`, `can_edit`, `can_manage` over `(user, room, share,
  group_role)`.
- FastAPI composites: `AccessReadDep`, `AccessEditDep`, `AccessManageDep`
  bundle room + share + group-role resolution and enforce the predicates.
- Unauthorized reads of `PRIVATE` / `GROUP` rooms return `404 Not Found`
  (not `403`) to avoid leaking existence.
- The per-room admin lock is deleted. The Redis per-geometry edit lock
  stays as pure serialization (no authorization).

#### Auth

- Every REST route that previously accepted `OptionalUserDep` now
  requires a JWT (`CurrentUserDep`). Guest tokens via
  `POST /v1/auth/guest` remain the canonical anonymous path.
- One exception: `GET /v1/rooms/{id}/trajectory` retains a local
  `OptionalUserDep` alias to support the pre-existing `zndraw-cli
  download` token-only flow.
- `get_local_token_or_admin` no longer depends on `OptionalUserDep`;
  JWT decode is inline.

#### New routes

```
POST   /v1/groups                              create group (creator → ADMIN)
GET    /v1/groups                              list groups the caller belongs to
GET    /v1/groups/{id}                         group details + my_role
PATCH  /v1/groups/{id}                         rename / describe (admin only)
DELETE /v1/groups/{id}                         delete (admin only; rooms must be reassigned first)

POST   /v1/groups/{id}/members                 add member (admin only)
DELETE /v1/groups/{id}/members/{user_id}       remove member (admin only, or self)
PATCH  /v1/groups/{id}/members/{user_id}       change role (admin only)

POST   /v1/rooms/{id}/share-links              create share link (manager only)
GET    /v1/rooms/{id}/share-links              list active links (manager only)
DELETE /v1/rooms/{id}/share-links/{link_id}    revoke link (manager only)
```

#### Modified routes

- `GET /v1/rooms` — list filter is now union of public ∪ owner-user ∪
  group-member (previously `is_public=True` only).
- `GET /v1/rooms/{id}` / `presence` / `sessions` — gated on `can_read`.
- `PATCH /v1/rooms/{id}` — gained `owner_user_id`, `owner_group_id`,
  `visibility`, `description` as updatable fields. Ownership transfer is
  a PATCH that sets exactly one owner field; server enforces invariants
  and the caller-must-be-in-target-group rule for group transfers. The
  `locked` field is gone.
- Content routes (geometries, figures, frames, trajectory, bookmarks,
  selection groups, presets, step) gate on `can_read` or `can_edit` via
  the new composites.

#### Socket.IO

- `on_connect` stashes a `share_token` field from the auth payload onto
  the session (alongside `user_id` and `current_room_id`).
- `room_join` replaces the `RoomMembership` lookup with `can_read`.
  Share tokens from the connect-time auth payload participate. Denial
  raises `RoomNotFound` (404) to avoid existence leaks.
- Broadcast routing (`broadcast_room_update`) splits by visibility:
  - `PUBLIC` → `rooms:feed`
  - `GROUP` → every group member's `user:{uid}` channel
  - `PRIVATE` → the owner's `user:{owner_user_id}` channel

#### Problem types

Added (RFC 9457):
- `GroupNotFound` (404) — group missing or not visible to caller
- `GroupNameTaken` (409) — duplicate group name
- `NotGroupMember` (403) — caller must be any-role member
- `NotGroupAdmin` (403) — admin role required
- `LastGroupAdmin` (409) — cannot demote/remove sole admin
- `GroupHasRooms` (409) — group owns rooms; delete or reassign first
- `TransferTargetInvalid` (409) — caller not in the target group
- `ShareLinkNotFound` (404) — revoke/list target missing or revoked
- `ShareLinkInvalid` (401) — `X-Room-Share-Token` unknown / revoked /
  expired / wrong room (documentation-only today; resolver returns
  `None` instead of raising)

Removed:
- `NotRoomMember`
- `AlreadyRoomMember`

#### Frontend

- New types in `frontend/src/myapi/client.ts`: `Visibility`, `GroupRole`,
  `ShareAccess`, `Group`, `GroupMember`, `ShareLink`.
- Room interfaces drop `locked`; add `visibility`, `owner_user_id`,
  `owner_group_id`.
- New API clients: `listGroups`, `createGroup`, `getGroup`, `updateGroup`,
  `deleteGroup`, `listGroupMembers`, `addGroupMember`, `removeGroupMember`,
  `updateGroupMemberRole`, `listShareLinks`, `createShareLink`,
  `revokeShareLink`.
- New components: `VisibilitySelector`, `ShareDialog`.
- New page: `/groups` (`GroupsPage` — list my groups + create form).
- `?share={token}` URL param is parsed on room-view mount; stored
  per-tab in memory (not localStorage); attached to REST as
  `X-Room-Share-Token` via an axios interceptor, and to socket.io as
  `share_token` in the connect auth payload.
- The lock/unlock menu items in `RoomManagementMenu`, `RoomsPanel`, and
  `roomRowMenu` are gone. Visibility is now controlled via the
  `VisibilitySelector` inside `RoomManagementMenu`.

#### Settings

- Added `ZNDRAW_SERVER_DEFAULT_ROOM_VISIBILITY` (default `public`) —
  controls the visibility assigned when `POST /v1/rooms` omits it.

#### Non-goals

- External auth (OAuth / OIDC / SSO)
- Automated email verification (superuser can still flip
  `is_verified` manually via `PATCH /v1/auth/users/{id}`)
- Per-room individual ACLs (`RoomCollaborator`) — share links + groups
  cover current use cases
- Cross-group room sharing (a room has exactly one owner)
