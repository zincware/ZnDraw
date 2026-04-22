# Room / User / Group scope refactor

**Branch:** `feat/room-group-scope-refactor`

## Context

Today's room access model has three collapsed concerns:

- `Room.is_public: bool` — binary visibility.
- `Room.locked: bool` — admin lock that blocks non-owner edits.
- `RoomMembership(room_id, user_id, role)` — per-room ACL with a
  `MemberRole` enum (`MEMBER | MODERATOR | OWNER`) that is mostly
  inert (no code path enforces MODERATOR/OWNER distinctions today).

Visibility is binary and per-room ACLs are the only way to grant
narrower access, which does not scale to "a team of 5 people who share
20 rooms." There is also no concept of user grouping on the principal
side: `User.is_superuser` is the sole role, and groups/orgs/teams do
not exist (the only "group" in the repo is `SelectionGroup`, unrelated
to access).

Anonymous access is implemented through a mixture of `OptionalUserDep`
(14 call sites) and a real guest-user row created by
`POST /v1/auth/guest` (`src/zndraw/routes/auth.py:25-43`). The guest
row is indistinguishable from a registered user except for an email
suffix (`{8-hex}@guest.user`). Every client — CLI, frontend — already
obtains a guest JWT before making API calls, so `user is None` paths
are dead weight that also forces every route to carry a "maybe auth"
branch.

The refactor introduces **users + groups + polymorphic room owners +
visibility scope + share links**, replaces the lock with pure
ownership semantics, and eliminates the anonymous code path.

## Non-goals

- External auth providers (OAuth, OIDC, institutional SSO). The
  design anchors every access edge on `User.id`; provider integration
  is a future additive change.
- Email verification flow. `is_verified` stays a column on User; a
  new superuser-only route lets admins flip it manually. Automated
  email verification is out of scope.
- Per-room individual ACLs (`RoomCollaborator`). Share links cover
  ad-hoc sharing; groups cover durable teams. If the gap between them
  becomes a real problem, `RoomCollaborator` is an additive change.
- Per-user group-external resource sync (SAML group claims, GitHub
  orgs). Out of scope until external auth lands.
- Cross-group room sharing (one room shared with multiple groups).
  A room has exactly one owner; groups do not overlap.

## Design

### Data model

```python
class Visibility(str, Enum):
    PRIVATE = "private"   # only the owning user (user-owned rooms)
    GROUP   = "group"     # only members of the owning group (group-owned rooms)
    PUBLIC  = "public"    # anyone (including guests) can view


class GroupRole(str, Enum):
    VIEWER = "viewer"   # read-only for group rooms
    MEMBER = "member"   # read + edit for group rooms
    ADMIN  = "admin"    # member + manage membership + manage group rooms


class ShareAccess(str, Enum):
    VIEW = "view"
    EDIT = "edit"


class Room(SQLModel, table=True):
    __table_args__ = (
        CheckConstraint(
            "(owner_user_id IS NOT NULL) <> (owner_group_id IS NOT NULL)",
            name="room_owner_exactly_one",
        ),
        CheckConstraint(
            "(visibility = 'private' AND owner_user_id IS NOT NULL) OR "
            "(visibility = 'group'   AND owner_group_id IS NOT NULL) OR "
            "(visibility = 'public')",
            name="room_visibility_matches_owner",
        ),
    )
    id: str = Field(default_factory=lambda: str(uuid4()), primary_key=True)
    description: str | None = None
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime())
    created_by_id: UUID | None = Field(default=None, index=True)  # audit, immutable
    owner_user_id:  UUID | None = Field(default=None, foreign_key="user.id",  index=True)
    owner_group_id: UUID | None = Field(default=None, foreign_key="group.id", index=True)
    visibility: Visibility = Field(default=Visibility.PUBLIC)
    step: int = Field(default=0)
    frame_selection: str | None = Field(default=None)
    default_camera: str | None = Field(default=None)


class Group(SQLModel, table=True):
    id: UUID = Field(default_factory=uuid4, primary_key=True)
    name: str = Field(unique=True, index=True)
    description: str | None = None
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime())
    created_by_id: UUID = Field(foreign_key="user.id", index=True)


class GroupMembership(SQLModel, table=True):
    __table_args__ = (UniqueConstraint("group_id", "user_id"),)
    id: int | None = Field(default=None, primary_key=True)
    group_id: UUID = Field(foreign_key="group.id", index=True)
    user_id:  UUID = Field(foreign_key="user.id",  index=True)
    role: GroupRole = Field(default=GroupRole.VIEWER)
    joined_at: datetime = Field(default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime())


class RoomShareLink(SQLModel, table=True):
    id: UUID = Field(default_factory=uuid4, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    token: str = Field(unique=True, index=True)  # url-safe, ~32 bytes
    access: ShareAccess = Field(default=ShareAccess.VIEW)
    created_by_id: UUID = Field(foreign_key="user.id")
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC), sa_type=UTCDateTime())
    expires_at: datetime | None = Field(default=None, sa_type=UTCDateTime())
    revoked_at: datetime | None = Field(default=None, sa_type=UTCDateTime())
```

**Removals:** `Room.is_public`, `Room.locked`, `RoomMembership`,
`MemberRole`.

**zndraw-auth additions:** `User.is_guest: bool = False`. Set to
`True` in `create_guest_session`. Used instead of email-suffix
matching; survives future external auth changes cleanly.

### Access control

All permission checks derive from ownership and group role. There is
no lock, no capability override except share tokens (bearer tokens
scoped to one room).

The `can_read` predicate is the single source of truth; REST and
socketio both use it. Unauthorized access to `PRIVATE` or `GROUP`
rooms returns `404 Not Found` (not `403`) to avoid leaking room
existence.

```python
def can_read(user: User, room: Room, share: ShareContext | None) -> bool:
    # `share`, when not None, is guaranteed valid + room-matched by resolver
    if user.is_superuser:
        return True
    if room.visibility == Visibility.PUBLIC:
        return True
    if room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None:
        if user_in_group(user.id, room.owner_group_id):
            return True
    if share is not None:
        return True  # any non-None ShareContext implies at least VIEW access
    return False


def can_edit(user: User, room: Room, share: ShareContext | None) -> bool:
    if user.is_superuser:
        return True
    if room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None:
        role = group_role(user.id, room.owner_group_id)
        if role in (GroupRole.MEMBER, GroupRole.ADMIN):
            return True
    if room.owner_user_id is not None and room.visibility == Visibility.PUBLIC:
        return True  # chaotic-edit default for user-owned public rooms
    if share is not None:
        return share.access == ShareAccess.EDIT
    return False


def can_manage(user: User, room: Room) -> bool:
    """Delete, transfer, change visibility. Share tokens never suffice."""
    if user.is_superuser:
        return True
    if room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None:
        return group_role(user.id, room.owner_group_id) == GroupRole.ADMIN
    return False
```

Read/edit matrix:

| Room                       | VIEWER | MEMBER | ADMIN | Owner (user) | Non-member (authed) | Guest | Superuser |
|----------------------------|--------|--------|-------|--------------|---------------------|-------|-----------|
| user-owned, private        | —      | —      | —     | R/E          | —                   | —     | all       |
| user-owned, public         | —      | —      | —     | R/E          | R/E                 | R/E   | all       |
| group-owned, group-visible | R      | R/E    | R/E/M | —            | —                   | —     | all       |
| group-owned, public        | R      | R/E    | R/E/M | —            | R                   | R     | all       |

(R = read, E = edit, M = manage.) Share-link bearers get R or R/E
according to the link's `access`, but never M.

### Listing queries

"Rooms visible to me" is a single SQL predicate:

```python
my_group_ids = select(GroupMembership.group_id).where(GroupMembership.user_id == me.id)

select(Room).where(
    or_(
        Room.visibility == Visibility.PUBLIC,
        Room.owner_user_id == me.id,
        Room.owner_group_id.in_(my_group_ids),
    )
)
```

Existing indexes on `owner_user_id` and `owner_group_id` carry the
weight. Group-member cache (`MyGroupIdsDep`) avoids the subquery per
request.

### Share links

Bearer-token capability grants that compose with visibility. URL
carries `?share={token}`. Frontend:

1. Parses `share` from the URL on navigation.
2. Holds it in memory for the tab (not localStorage — ephemeral).
3. Attaches `X-Room-Share-Token: {token}` to every REST call for that
   room's namespace.
4. Passes it in the socketio `join` auth payload.

Server middleware (`ShareTokenDep`) resolves the token: looks up the
row by `token`, verifies `revoked_at IS NULL` and
`expires_at IS NULL OR expires_at > now()`, returns a `ShareContext`
(room_id, access) — or `None` if the token is missing, revoked,
expired, or belongs to a different room. The permission checks above
consume it; absence reduces to "no share context" without raising.

Share tokens grant view/edit, never management. Revoking a link is
setting `revoked_at` (preferred over row delete, keeps audit).

### Concurrency coordination (orthogonal to permissions)

The existing per-geometry edit lock
(`dependencies.py:293-330`) is runtime coordination, not permission,
and stays. The rule going forward:

- Lock key: `lock:{room_id}:{resource}` in Redis, TTL ~30s, refreshed
  on activity.
- Holder: `user_id + session_id`, stored as the value.
- Auto-expire if the TTL runs out (holder went silent).
- The permission check and the lock acquisition are separate calls.
  The lock never authorizes; it only serializes among already-
  authorized writers.

**Deferred (follow-up, not this refactor):** explicit force-release
endpoint so any caller with `can_edit` on the same resource can steal
an abandoned lock without waiting for TTL. TTL-based expiry covers
the common case; force-release is a UX enhancement for the next
iteration.

### API surface

All routes follow REST resource conventions — nouns, not verbs.
State transitions are `PATCH`es on the resource; collections are
plural nouns; writes return the updated resource.

**New routes:**

```
POST   /v1/groups                              create group (creator → ADMIN)
GET    /v1/groups                              list groups the caller belongs to
GET    /v1/groups/{id}                         group details + member list (any group member, regardless of role)
PATCH  /v1/groups/{id}                         rename / describe (admin only)
DELETE /v1/groups/{id}                         delete (admin only; rooms must be reassigned first)

POST   /v1/groups/{id}/members                 add member (admin only)
DELETE /v1/groups/{id}/members/{user_id}       remove member (admin only, or self)
PATCH  /v1/groups/{id}/members/{user_id}       change role (admin only)

POST   /v1/rooms/{id}/share-links              create share link
GET    /v1/rooms/{id}/share-links              list active links (manager only)
DELETE /v1/rooms/{id}/share-links/{link_id}    revoke link (manager only)
```

**Modified routes (state transitions as PATCH on the resource):**

- `PATCH /v1/rooms/{id}` gains `owner_user_id`, `owner_group_id`,
  `visibility`, `description` as updatable fields. Ownership transfer
  is just "set the owner field; server validates target exists + the
  caller is allowed + invariants hold." Body is a partial — any
  subset of the above. Server enforces:
  - exactly-one-owner invariant (XOR via CHECK);
  - caller satisfies `can_manage(room)`;
  - if `owner_group_id` is set, caller is a MEMBER or ADMIN of that
    group (you can transfer into a group you're in);
  - if `owner_user_id` is changed to another user, transfer is
    one-sided (effective immediately, no acceptance step). A future
    "recipient must accept" flow is a deliberate non-goal of this
    refactor.
- `PATCH /v1/users/{id}` — this endpoint already exists via
  fastapi-users. `is_verified` is editable only by superusers; other
  fields follow fastapi-users' existing authorization rules. No new
  route; new server-side policy.
- Every `OptionalUserDep` call site → `CurrentUserDep`. Site list in
  the Refactor sites section below.
- `GET /v1/rooms` list filter: union of public / owner-user /
  group-member, replacing today's `is_public=True` filter.
- `GET /v1/rooms/{id}` detail: uses `can_read` (supports share
  tokens).
- Content routes (geometry, figures, trajectory, frames) gate on
  `can_edit`.
- Delete gates on `can_manage`.

**Socketio:**

- `room_join` (`socketio.py:171-323`): replaces the
  `RoomMembership`-based check with `can_read`. Share token passed in
  the `auth` payload at connect, stored on the session, attached to
  permission checks.
- `room_leave` unchanged.
- Existing emits (`SessionJoined`, `GeometryInvalidate`, etc.) are
  unchanged; the room-channel pubsub (`room:{room_id}`) still scopes
  fan-out.

### Error types (RFC 9457 Problem Details)

All new errors follow the existing `ProblemType` pattern in
`src/zndraw/exceptions.py` (kebab-case problem IDs under
`/v1/problems/`, `ProblemDetail` responses, `.exception()` raisers
wired through the global handler). Reuse existing types where they
fit; add new ones below.

**Reused:**

- `Forbidden` (403) — caller lacks required role for the operation.
  Used when `can_edit` / `can_manage` fails but the caller can still
  see the room.
- `RoomNotFound` (404) — room ID does not exist, OR caller cannot
  `can_read` a PRIVATE / GROUP room (404 over 403 to avoid leaking
  existence).
- `NotAuthenticated` (401) — raised by fastapi-users' JWT strategy
  on any route that uses `CurrentUserDep` when the token is missing,
  malformed, or the user is inactive. Also raised by
  `get_local_token_or_admin` when neither auth path succeeds.
- `InvalidPayload` (422) — transfer requests with both owner fields
  set, unknown enum values, etc.

**New (to add in `exceptions.py`):**

| Type | Status | When |
|---|---|---|
| `GroupNotFound` | 404 | group ID unknown, or caller cannot see it |
| `GroupNameTaken` | 409 | `POST /v1/groups` with a name already in use |
| `NotGroupMember` | 403 | caller must be a member (any role) of the group to perform the op |
| `NotGroupAdmin` | 403 | caller must be ADMIN of the group |
| `LastGroupAdmin` | 409 | attempt to demote / remove the sole ADMIN of a group |
| `GroupHasRooms` | 409 | `DELETE /v1/groups/{id}` when one or more rooms are still owned by the group |
| `TransferTargetInvalid` | 409 | `PATCH /v1/rooms/{id}` changes owner_group_id to a group the caller is not in |
| `ShareLinkNotFound` | 404 | revoke / list target does not exist, expired, or already revoked |
| `ShareLinkInvalid` | 401 | provided `X-Room-Share-Token` is unknown, revoked, expired, or targets a different room |

Each new type inherits from `ProblemType`, defines `title` + `status`,
and provides `raise_for_client` for symmetric client-side mapping
matching existing conventions (e.g. `raise PermissionError` for 403,
`raise ValueError` for 409). Route decorators use `problem_responses(
...)` to register them in OpenAPI.

### Defaults and group-management rules

- Room creation: `visibility = Visibility.PUBLIC`,
  `owner_user_id = creator.id`, `owner_group_id = None`. Preserves
  the current `zndraw file.xyz` CLI workflow (chaotic-edit public
  demos). Configurable via `Settings.default_room_visibility`
  (Pydantic field).
- Group creation: any active user may create. Creator is sole
  `ADMIN`.
- Group names are **globally unique** — `Group.name` has a UNIQUE
  index. Collisions return `GroupNameTaken` (409).
- Group join: new members default to `VIEWER` (least privilege).
  Admins explicitly promote to `MEMBER` / `ADMIN`.
- Group deletion: **forbidden while any room is still owned by the
  group.** `DELETE /v1/groups/{id}` with group-owned rooms returns
  `GroupHasRooms` (409). Admins must reassign each room's owner
  (`PATCH /v1/rooms/{id}`) before deleting.
- **Last-admin protection.** The sole ADMIN of a group cannot
  self-remove, demote themselves, or otherwise leave the group admin-
  less. Returns `LastGroupAdmin` (409). To leave, they must first
  promote another member to ADMIN (or delete the group after
  reassigning rooms).
- Share links: `access = VIEW`, `expires_at = None`.

### Schema change (no data migration)

**There is no alembic migration and none is wanted.** The project is
pre-v1.0.0; schema breakage is acceptable. On the release that ships
this refactor, **users must recreate their database** — drop the old
SQLite/Postgres schema, let the app re-run `create_all` on startup
against the new models, and resume from a blank state. Existing
rooms, memberships, and ACLs are not preserved. Document this
prominently in the release notes / CHANGELOG.

This keeps the refactor surface clean: no backfill scripts, no
dual-read shim, no compatibility window.

### Refactor sites (code)

All call-site changes needed to realize the new model. Each of these
gets covered by the implementation plan.

- `OptionalUserDep` → `CurrentUserDep`:
  - `routes/rooms.py:430` (`list_rooms`)
  - `routes/geometries.py:97, 129, 149, 383`
  - `routes/figures.py:41, 59`
  - `routes/trajectory.py:93`
- `get_local_token_or_admin` (`dependencies.py:47-76`): stops
  depending on `OptionalUserDep`; chain its two paths with a
  dedicated local helper. Remove the `OptionalUserDep` export.
- `socketio.py:171-323` `room_join`: replace `RoomMembership` lookup
  with `can_read`, accept share token in auth payload.
- Geometry edit-lock logic (`dependencies.py:293-330`): permission
  gate becomes `can_edit`; the lock itself stays.
- Guest creation (`routes/auth.py:32`): set `is_guest=True` on the
  `UserCreate` payload. Remove any existing `@guest.user` suffix
  checks elsewhere.
- Room model: drop `is_public`, `locked`. Add `owner_user_id`,
  `owner_group_id`, `visibility`, the two CHECK constraints.
- `RoomMembership` + `MemberRole`: delete. Remove
  `NotRoomMember` / `AlreadyRoomMember` problem types.
- Frontend: update room-creation UI to include the three-value
  visibility selector (replacing the public checkbox), add the Share
  dialog with link management, add a Groups section, wire the
  `?share=` URL param → `X-Room-Share-Token` header.

### Capability and guest policy

- Guests (`is_guest=True`) are full users at the schema level. They
  own rooms, create groups, get invited. The only current capability
  gate is informational (`is_verified`); future email-verification
  gates can consult `is_guest` or `is_verified` as needed.
- Admins (`is_superuser=True`) can manually verify any user by
  `PATCH /v1/users/{id}` with `is_verified: true`. The superuser-
  only field-level check lives in the user-update handler (the route
  itself is from fastapi-users).

### Testing strategy

Integration tests against `server_factory` + real Redis (per
`CLAUDE.md`: "prefer testing against server_factory instead of
mockups"; no mock Redis).

Key scenarios:

- **Visibility matrix.** For each combination of (owner kind,
  visibility, caller relation, role), assert `can_read`, `can_edit`,
  `can_manage` return the expected values. Parameterized table test.
- **CHECK constraints.** Inserting a room with both owner FKs set,
  or visibility mismatched to owner, raises `IntegrityError`.
- **Ownership transfer via PATCH.** user→user, user→group, group→
  user, and combined visibility-change cases. Edge cases: target
  group caller is not in (expect `TransferTargetInvalid`);
  non-manager caller attempts transfer (expect `Forbidden`);
  invariant-violating payload (both owner fields) rejected as
  `InvalidPayload`.
- **Share links.** Create (view/edit), use by guest, expire, revoke,
  unknown token → `ShareLinkInvalid`, token for a different room →
  `ShareLinkInvalid`. Include a socketio join test with a share
  token in the auth payload.
- **Group membership.** Add/remove, role change, last-admin
  protection (`LastGroupAdmin` when demoting or removing the sole
  admin, including self-removal), default-role-on-add is VIEWER,
  duplicate-name creation → `GroupNameTaken`.
- **Group deletion.** Empty group deletable by admin; group with
  owned rooms → `GroupHasRooms`; non-admin → `Forbidden`.
- **Auth gate.** Every route that previously used `OptionalUserDep`
  now returns 401 on missing auth (previously returned 200 with
  partial data).
- **Problem-type OpenAPI coverage.** Each new `ProblemType` is
  referenced by at least one route's `responses=problem_responses(
  ...)` so the OpenAPI schema documents the error envelope.

## Interaction with existing features

- **Rooms feed** (`2026-04-21-rooms-feed-design.md`): broadcasts
  `room_update` to the `rooms:feed` channel. Under the new model,
  feed delivery must be scoped — clients only receive updates for
  rooms they `can_read`. The channel stays one; delivery is filtered
  server-side before emit (per-subscriber filter) or clients
  subscribe to targeted channels (`rooms:feed:user:{id}`,
  `rooms:feed:group:{id}`, `rooms:feed:public`). The second approach
  scales better; defer the decision to implementation.
- **Default room / `@empty` preset**: seeded as user-owned public (by
  `DEFAULT_ADMIN`), matches today.
- **Joblib workers** (`/v1/joblib/rooms/@global/...`): the `@global`
  scope is server-side infrastructure, not a Room row; unaffected.
- **RoomGeometry, RoomBookmark, RoomFigure, RoomPreset** tables:
  unaffected by this refactor; they key on `room_id` and their
  existing access paths go through the permission layer above.
