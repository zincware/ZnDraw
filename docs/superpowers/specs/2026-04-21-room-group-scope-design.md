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
    created_at: datetime = Field(default_factory=..., sa_type=UTCDateTime())
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
    created_at: datetime = Field(default_factory=..., sa_type=UTCDateTime())
    created_by_id: UUID = Field(foreign_key="user.id", index=True)


class GroupMembership(SQLModel, table=True):
    __table_args__ = (UniqueConstraint("group_id", "user_id"),)
    id: int | None = Field(default=None, primary_key=True)
    group_id: UUID = Field(foreign_key="group.id", index=True)
    user_id:  UUID = Field(foreign_key="user.id",  index=True)
    role: GroupRole = Field(default=GroupRole.VIEWER)
    joined_at: datetime = Field(default_factory=..., sa_type=UTCDateTime())


class RoomShareLink(SQLModel, table=True):
    id: UUID = Field(default_factory=uuid4, primary_key=True)
    room_id: str = Field(foreign_key="room.id", index=True)
    token: str = Field(unique=True, index=True)  # url-safe, ~32 bytes
    access: ShareAccess = Field(default=ShareAccess.VIEW)
    created_by_id: UUID = Field(foreign_key="user.id")
    created_at: datetime = Field(default_factory=..., sa_type=UTCDateTime())
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
    if user.is_superuser:
        return True
    if room.visibility == Visibility.PUBLIC:
        return True
    if room.owner_user_id == user.id:
        return True
    if room.owner_group_id is not None:
        if user_in_group(user.id, room.owner_group_id):
            return True
    if share is not None and share.room_id == room.id and share.valid:
        return share.access in (ShareAccess.VIEW, ShareAccess.EDIT)
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
    if share is not None and share.room_id == room.id and share.valid:
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
and stays. Formalize the rule:

- Lock key: `lock:{room_id}:{resource}` in Redis, TTL ~30s, refreshed
  on activity.
- Holder: `user_id + session_id`, stored as the value.
- Auto-expire if the TTL runs out (holder went silent).
- Force-release allowed: any caller whose `can_edit` check passes for
  the same resource may steal a stale lock. The UI shows "X is
  editing — take over?" after a threshold.

The permission decision and the coordination decision are separate
calls. The lock never authorizes; it only serializes among already-
authorized writers.

### API surface

**New routes:**

```
POST   /v1/groups                         create group (creator → ADMIN)
GET    /v1/groups                         list groups the caller belongs to
GET    /v1/groups/{id}                    group details + member list (members+admins only)
PATCH  /v1/groups/{id}                    rename / describe (admin only)
DELETE /v1/groups/{id}                    delete (admin only; rooms must be transferred first)

POST   /v1/groups/{id}/members            add member (admin only)
DELETE /v1/groups/{id}/members/{user_id}  remove member (admin only, or self)
PATCH  /v1/groups/{id}/members/{user_id}  change role (admin only)

POST   /v1/rooms/{id}/transfer            transfer ownership
                                          body: {to_user_id?, to_group_id?, new_visibility?}
POST   /v1/rooms/{id}/share               create share link; returns {token, url, access, expires_at?}
GET    /v1/rooms/{id}/share               list active links (manager only)
DELETE /v1/rooms/{id}/share/{link_id}     revoke link (manager only)

PATCH  /v1/users/{id}/verify              flip is_verified (superuser only)
```

**Modified routes:**

- Every `OptionalUserDep` call site → `CurrentUserDep`. List of sites
  in Migration below.
- `GET /v1/rooms` list filter: union of public / owner-user / group-member,
  replacing today's `is_public=True` filter.
- `GET /v1/rooms/{id}` detail: uses `can_read` (supports share
  tokens).
- Content routes (geometry, figures, trajectory, frames) gate on
  `can_edit`.
- Delete / visibility change / transfer gate on `can_manage`.

**Socketio:**

- `room_join` (`socketio.py:171-323`): replaces the
  `RoomMembership`-based check with `can_read`. Share token passed in
  the `auth` payload at connect, stored on the session, attached to
  permission checks.
- `room_leave` unchanged.
- Existing emits (`SessionJoined`, `GeometryInvalidate`, etc.) are
  unchanged; the room-channel pubsub (`room:{room_id}`) still scopes
  fan-out.

### Defaults

- Room creation: `visibility = Visibility.PUBLIC`,
  `owner_user_id = creator.id`, `owner_group_id = None`. Preserves the
  current `zndraw file.xyz` CLI workflow (chaotic-edit public demos).
  Configurable via `Settings.default_room_visibility` (Pydantic
  field).
- Group creation: any active user may create. Creator is sole
  `ADMIN`.
- Group join: new members default to `VIEWER` (least privilege).
  Admins explicitly promote to `MEMBER` / `ADMIN`.
- Share links: `access = VIEW`, `expires_at = None`.

### Migration

Schema migration is additive then subtractive; backfill happens
between the two phases.

**Phase 1 — additive (Alembic revision 1):**

1. Create `Group`, `GroupMembership`, `RoomShareLink` tables.
2. Add `Room.owner_user_id`, `Room.owner_group_id`, `Room.visibility`
   (all nullable).
3. Add `User.is_guest` to zndraw-auth (default False).

**Phase 2 — backfill (data migration script, runs between the two
alembic revisions):**

1. `UPDATE room SET owner_user_id = created_by_id`
   for rows where `created_by_id IS NOT NULL`.
2. `UPDATE room SET visibility = CASE WHEN is_public THEN 'public'
   ELSE 'private' END`.
3. `UPDATE user SET is_guest = TRUE WHERE email LIKE '%@guest.user'`.
4. For any room where `created_by_id IS NULL` (orphaned): assign to
   a bootstrap admin user. The implementation step MUST verify the
   bootstrap user exists (env `DEFAULT_ADMIN_EMAIL` → lookup); if
   none is configured, the migration fails loudly rather than
   fabricating ownership.
5. **`RoomMembership` data loss.** The existing `RoomMembership`
   table holds the only per-user ACL for private rooms. Under the
   new model there is no per-room ACL. Two choices:
   - **Drop silently.** Acceptable only if `RoomMembership` is
     known-empty or known-redundant in production (the role enum is
     underutilized today — no MODERATOR/OWNER distinction is
     enforced anywhere, and ownership is already tracked in
     `created_by_id`). Document in release notes.
   - **Reconstitute as groups.** For each room with >1 distinct
     `user_id` in `RoomMembership` (beyond the creator), auto-create
     a `Group` named `room-{room_id[:8]}` with the members as MEMBER
     role (creator as ADMIN), set `owner_group_id = group.id`,
     `visibility = GROUP`. More work; preserves intent.

   Recommendation: inspect the prod DB row count before the
   implementation plan locks this in. If there are zero meaningful
   non-creator memberships, drop. Otherwise reconstitute.

**Phase 3 — subtractive (Alembic revision 2):**

1. Add `CheckConstraint`s on `Room` (owner XOR, visibility matches owner).
2. `ALTER COLUMN Room.visibility SET NOT NULL`.
3. Drop `Room.is_public`, `Room.locked`.
4. Drop `RoomMembership` table.
5. Drop `MemberRole` enum.

**Code migration:**

- `OptionalUserDep` call sites — rewrite to `CurrentUserDep`:
  - `routes/rooms.py:430` (`list_rooms`)
  - `routes/geometries.py:97, 129, 149, 383`
  - `routes/figures.py:41, 59`
  - `routes/trajectory.py:93`
- `get_local_token_or_admin` (`dependencies.py:47-76`): stops depending
  on `OptionalUserDep`; chain its two paths with a dedicated local
  helper. Remove the `OptionalUserDep` export.
- `socketio.py:171-323` `room_join`: replace `RoomMembership` lookup
  with `can_read`, accept share token in auth payload.
- Geometry edit-lock logic (`dependencies.py:293-330`): permission
  gate becomes `can_edit`; the lock itself stays. Add force-release
  endpoint.

### Capability and guest policy

- Guests (`is_guest=True`) are full users at the schema level. They
  own rooms, create groups, get invited. The only current capability
  gate is informational (`is_verified`); future email-verification
  gates can consult `is_guest` or `is_verified` as needed.
- Admins (`is_superuser=True`) can manually verify any user via
  `PATCH /v1/users/{id}/verify`.

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
- **Transfer flows.** user→user, user→group, group→user, all edge
  cases including non-member target group rejection and
  already-ADMIN preservation.
- **Share links.** Create (view/edit), use by guest, expire, revoke,
  wrong-token 404, token for a different room 404. Include a socketio
  join test with a share token in the auth payload.
- **Group membership.** Add/remove, role change, last-admin
  protection (cannot remove or demote the sole admin),
  default-role-on-add is VIEWER.
- **Guest promotion.** Guest creates room → registers with email →
  room ownership preserved (future-proofing assertion on linking
  behavior, even though external auth is not in this scope).
- **Migration determinism.** A snapshot of a pre-migration DB (with
  `is_public`, `locked`, `RoomMembership` rows) runs through Phase 1
  + backfill + Phase 3 and lands on identical permission decisions as
  the old code for the same caller/room pairs.
- **`OptionalUserDep` removal.** Every route that used it now returns
  401 on missing auth (previously returned 200 with partial data).

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

## Open questions

- **Group deletion with rooms present.** Forbid (require transfer
  first), or cascade orphan rooms to the ADMIN who deletes the group?
  Recommendation: forbid. Matches GitHub org-delete semantics.
- **Group name uniqueness.** Global vs per-creator namespace? Global
  simpler, but users may compete for common names ("research",
  "team-a"). Recommendation: global with a friendly slug collision
  error. Worth revisiting if friction emerges.
- **Last-admin protection.** If a group has one ADMIN and they leave,
  what happens? Recommendation: block self-removal when sole admin;
  require either promoting another member or deleting the group.
