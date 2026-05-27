# Coolname display names Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace UUID-shaped owner segments in URLs and chat payloads with auto-generated coolname display names so users see friendly identifiers like `happy-blue-rabbit/my-room` instead of `a1b2c3d4-…/my-room`.

**Architecture:** Add a unique `display_name` column on `User` (auto-filled with `coolname.generate_slug(3)`), translate the two-segment route's first path param from `UUID` to a regex-validated display-name string, resolve that string back to a user/group UUID at the dependency boundary, and surface the display name (instead of email/UUID) in REST + socket payloads. No migration, no backwards compat — break as you go.

**Tech Stack:** FastAPI, SQLModel, fastapi-users, coolname, Socket.IO, React + Zustand + TypeScript.

---

## Spec reference

`docs/superpowers/specs/2026-05-27-coolname-display-names-design.html` (issue #931).

## File map

### Backend — new
- `src/zndraw_auth/display_names.py` — coolname slug generator, regex + reserved-word validator.
- `src/zndraw/routes/users.py` — `GET /v1/users/available-display-name` suggestion endpoint.
- `tests/zndraw_auth/test_display_names.py` — unit tests for slug generator + validator.
- `tests/zndraw_auth/test_user_create.py` — registration with explicit / omitted / duplicate / malformed / reserved display name.

### Backend — modified
- `src/zndraw_auth/db.py` — `User.display_name` column.
- `src/zndraw_auth/schemas.py` — `UserRead.display_name`, `UserCreate.display_name`.
- `src/zndraw_auth/users.py` — `UserManager.create` override fills/validates display_name.
- `src/zndraw_auth/__init__.py` — re-export `display_names` helpers + `UserNotFound`-adjacent exception (existing in `zndraw.exceptions`; auth keeps its own neutral import).
- `src/zndraw/routes/auth.py` — guest session returns `GuestSessionResponse` (typed) with display_name.
- `src/zndraw/dependencies.py` — `owner: str = Path(pattern=...)` everywhere; new `get_owner_uuid_from_segment` dep; `resolve_owner` returns display_name for users; `verify_room` and `get_writable_room_id` parse display-name composed addresses.
- `src/zndraw/routes/rooms.py` — `RoomCreate.owner` (display name), `RoomPatchRequest.new_owner` (display name); `_resolve_room_owner` returns display-name label; address resolution at room CRUD entry points.
- `src/zndraw/routes/chat.py` — message author label changes from `email` → `display_name`; path params widened.
- `src/zndraw/routes/bookmarks.py`, `edit_lock.py`, `figures.py`, `frames.py`, `geometries.py`, `isosurface.py`, `presets.py`, `progress.py`, `screenshots.py`, `selection_groups.py`, `share_links.py`, `step.py`, `trajectory.py` — flip `owner_id: UUID = Path()` → `owner: str = Path(pattern=...)` so FastAPI's regex gate fires before any DB lookup.
- `src/zndraw/schemas.py` — `RoomCreate.owner_id` → `RoomCreate.owner: str`; `RoomPatchRequest.new_owner_id` → `RoomPatchRequest.new_owner: str`; `RoomResponse.owner_id` → `RoomResponse.owner: str`; `MessageResponse.email` → `MessageResponse.display_name`; `PresenceSessionResponse.email` → `display_name`; `SessionItem.email` → `display_name`; `GroupMemberResponse.email` → `display_name`.
- `src/zndraw/socket_events.py` — `RoomJoin.owner_id` → `RoomJoin.owner`; `RoomLeave.owner_id` → `RoomLeave.owner`; `TypingStart.owner_id` → `TypingStart.owner`; `TypingStop.owner_id` → `TypingStop.owner`; `UserGetResponse.email` → `UserGetResponse.display_name`; `MessageNew.email` → `MessageNew.display_name`; `SessionJoined.email` → `SessionJoined.display_name`; `Typing.email` → `Typing.display_name`.
- `src/zndraw/socketio.py` — emits/reads display_name in `user_get`, `room_join` (camera key + entry), typing handlers.
- `src/zndraw/models.py` — `Room.public_address` stays UUID-prefixed by primary key fallback; add module-level helper `async def build_public_address(session, room) -> str` returning `{display_name}/{room_name}` for owner_user, `{group_name}/{room_name}` for owner_group.
- `src/zndraw/cli.py` — `_resolve_owner_id` → `_resolve_owner_name` returns display name; `_validate_room_arg` checks display-name regex; help text references "display name", not "UUID".

### Frontend — modified
- `frontend/src/utils/auth.ts` — `UserInfo.display_name: string`; guest POST reads `display_name` from response.
- `frontend/src/stores/slices/connectionSlice.ts` — no shape change (UserInfo already imported).
- `frontend/src/store.tsx` — `userEmail` selector renamed to `userDisplayName`; ownership comparison switches to display_name.
- `frontend/src/components/RegisterDialog.tsx` — adds display-name input pre-filled via `/v1/users/available-display-name`, regenerate button, client-side regex validation, 409 error mapping.
- `frontend/src/components/UserProfileDialog.tsx` — show display_name (+ email below as secondary).
- `frontend/src/panels/ChatPanel.tsx` — author label + ownership check on `display_name`; typing emits send `owner` instead of `owner_id`.
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts` — `room_join` payload key renamed `owner_id` → `owner`.
- `frontend/src/myapi/client.ts` — `CreateRoomRequest.owner_id` → `CreateRoomRequest.owner`; `RoomInfo.owner_id` → `RoomInfo.owner`; `Room.owner_id` → `Room.owner`; `RoomUpdateRequest.new_owner_id` → `RoomUpdateRequest.new_owner`; `GroupMember.email` → `GroupMember.display_name`.
- `frontend/src/types/chat.ts` — `ChatMessage.email` → `display_name`; `MessageNewEvent.email` → `display_name`.
- `frontend/src/panels/RoomsPanel.tsx`, `panels/roomsHeaderActions.tsx`, `panels/FilesystemPanel.tsx`, `components/DuplicateRoomDialog.tsx`, `pages/templateSelection.tsx` — `owner_id: currentUser.id` → `owner: currentUser.display_name`.
- `frontend/src/roomsStore.tsx` — `updates.owner_id` → `updates.owner` field rename.

### Tests — new / extend
- New: `tests/zndraw_auth/test_display_names.py`, `tests/zndraw_auth/test_user_create.py`.
- Extend: `tests/zndraw/test_auth_endpoints.py`, `tests/zndraw/test_e2e_guest_auth.py`, `tests/zndraw/test_routes_rooms.py` (and any `test_routes_*.py` using owner segments), `tests/zndraw/test_owner_resolution.py`, `tests/zndraw/test_socket_commands.py`, `tests/zndraw/test_cli.py`, `frontend/e2e/`.

## Pre-flight (read once, before Task 1)

1. Confirm worktree on branch `worktree-coolname-display-names` (created via EnterWorktree).
2. Run `grep -rn "@pytest.mark.protected" tests/` — current matches: `tests/zndraw/test_event_schema_audit.py`, `tests/zndraw/test_broadcast_contract.py`. Both check `RoomScopedEvent` structural invariants (room_id/room_address) and per-route broadcast contracts. **Neither touches owner UUIDs or the email payload field** — they remain unchanged. If a later task surfaces a protected test that does depend on the old UUID-owner contract, **stop and ask the user**; do not edit it.
3. Install dependency: `uv add coolname`.
4. Baseline: `uv run pytest -x -q` and `cd frontend && bun run build` (or `bun run typecheck`). If the baseline is red, stop and report before changing anything.

---

### Task 1: Add coolname dependency + helper module

**Files:**
- Modify: `pyproject.toml`
- Modify: `tests/zndraw_auth/conftest.py` (add `session` fixture)
- Create: `src/zndraw_auth/display_names.py`
- Create: `tests/zndraw_auth/test_display_names.py`

- [ ] **Step 1: Add coolname to project dependencies**

Run: `uv add coolname`

Expected: `pyproject.toml` shows `coolname = "*"` under `[project.dependencies]` (or equivalent), `uv.lock` updated.

- [ ] **Step 1b: Add a reusable `session` fixture**

Append to `tests/zndraw_auth/conftest.py` (after the existing `client` fixture, before any helper-only blocks):

```python
@pytest.fixture
async def session(app: FastAPI) -> AsyncGenerator[AsyncSession, None]:
    """Yield an AsyncSession bound to the test app's engine."""
    session_maker = app.state.session_maker
    async with session_maker() as s:
        yield s
```

(`AsyncSession`, `AsyncGenerator`, `FastAPI`, and `pytest` are already imported at the top of the file.)

- [ ] **Step 2: Write the failing test file**

Create `tests/zndraw_auth/test_display_names.py`:

```python
"""Unit tests for display-name generation and validation."""

from __future__ import annotations

import pytest
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import ProblemError
from zndraw_auth.db import User
from zndraw_auth.display_names import (
    DISPLAY_NAME_PATTERN,
    RESERVED_DISPLAY_NAMES,
    generate_unique_display_name,
    validate_display_name,
)


def test_pattern_accepts_valid_names() -> None:
    for name in ("happy-blue-rabbit", "abc", "a1b-2c", "z" * 64):
        assert DISPLAY_NAME_PATTERN.fullmatch(name), name


@pytest.mark.parametrize(
    "name",
    ["", "AB", "AbCdEf", "1abc", "-abc", "ab", "ab!cd", "z" * 65, "_abc"],
)
def test_pattern_rejects_invalid_names(name: str) -> None:
    assert DISPLAY_NAME_PATTERN.fullmatch(name) is None


def test_reserved_set_is_lowercase_only() -> None:
    for token in RESERVED_DISPLAY_NAMES:
        assert token == token.lower()


@pytest.mark.parametrize("bad", ["me", "admin", "internal", "overview", "global"])
def test_validate_display_name_rejects_reserved(bad: str) -> None:
    with pytest.raises(ProblemError):
        validate_display_name(bad)


@pytest.mark.parametrize("bad", ["Happy", "1ab", "ab!cd", ""])
def test_validate_display_name_rejects_malformed(bad: str) -> None:
    with pytest.raises(ProblemError):
        validate_display_name(bad)


def test_validate_display_name_accepts_good() -> None:
    validate_display_name("happy-blue-rabbit")


@pytest.mark.asyncio
async def test_generate_unique_returns_regex_valid(session: AsyncSession) -> None:
    name = await generate_unique_display_name(session)
    assert DISPLAY_NAME_PATTERN.fullmatch(name), name
    assert name not in RESERVED_DISPLAY_NAMES


@pytest.mark.asyncio
async def test_generate_unique_avoids_existing(
    session: AsyncSession, monkeypatch: pytest.MonkeyPatch
) -> None:
    """If the first slug clashes with a user, the generator retries."""
    seq = iter(["happy-blue-rabbit", "merry-red-otter"])
    monkeypatch.setattr(
        "zndraw_auth.display_names.coolname.generate_slug",
        lambda n: next(seq),
    )
    session.add(
        User(
            email="taken@example.com",
            hashed_password="x",
            display_name="happy-blue-rabbit",
        )
    )
    await session.commit()
    name = await generate_unique_display_name(session)
    assert name == "merry-red-otter"


@pytest.mark.asyncio
async def test_generate_unique_falls_back_with_hex_suffix(
    session: AsyncSession, monkeypatch: pytest.MonkeyPatch
) -> None:
    """After max_attempts exhaustion the generator appends a hex suffix."""
    monkeypatch.setattr(
        "zndraw_auth.display_names.coolname.generate_slug",
        lambda n: "clashing-slug-name",
    )
    session.add(
        User(
            email="taken@example.com",
            hashed_password="x",
            display_name="clashing-slug-name",
        )
    )
    await session.commit()
    name = await generate_unique_display_name(session, max_attempts=3)
    assert name.startswith("clashing-slug-name-")
    assert DISPLAY_NAME_PATTERN.fullmatch(name)
```

- [ ] **Step 3: Run the tests and verify they fail**

Run: `uv run pytest tests/zndraw_auth/test_display_names.py -v`

Expected: ImportError or collection failure (module `zndraw_auth.display_names` doesn't exist yet).

- [ ] **Step 4: Implement `src/zndraw_auth/display_names.py`**

Create `src/zndraw_auth/display_names.py`:

```python
"""Display-name generation and validation."""

from __future__ import annotations

import re
import secrets
from typing import Final

import coolname
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth.db import User

DISPLAY_NAME_PATTERN: Final[re.Pattern[str]] = re.compile(r"^[a-z][a-z0-9-]{2,63}$")

RESERVED_DISPLAY_NAMES: Final[frozenset[str]] = frozenset(
    {"me", "admin", "internal", "overview", "global", "system", "guest"}
)


def validate_display_name(name: str) -> None:
    """Raise UnprocessableContent if ``name`` is malformed or reserved."""
    from zndraw.exceptions import UnprocessableContent

    if not DISPLAY_NAME_PATTERN.fullmatch(name) or name in RESERVED_DISPLAY_NAMES:
        raise UnprocessableContent.exception(f"Display name '{name}' is not allowed")


async def generate_unique_display_name(
    session: AsyncSession, *, max_attempts: int = 8
) -> str:
    """Generate a coolname slug unique against ``User.display_name``.

    Retries on reserved-word collision and DB collision. After
    ``max_attempts`` exhaust without progress, appends a 4-char hex
    suffix to guarantee forward progress.
    """
    for _ in range(max_attempts):
        slug = coolname.generate_slug(3)
        if slug in RESERVED_DISPLAY_NAMES:
            continue
        if not DISPLAY_NAME_PATTERN.fullmatch(slug):
            continue
        exists = await session.scalar(
            select(User.id).where(User.display_name == slug).limit(1)
        )
        if exists is None:
            return slug
    return f"{coolname.generate_slug(3)}-{secrets.token_hex(2)}"
```

- [ ] **Step 5: Run the tests and verify they pass**

Run: `uv run pytest tests/zndraw_auth/test_display_names.py -v`

Expected: PASS (10+ tests).

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml uv.lock src/zndraw_auth/display_names.py tests/zndraw_auth/test_display_names.py
git commit -m "feat(auth): coolname-based display-name generator + validator"
```

---

### Task 2: Add `User.display_name` column + UserCreate/UserRead schemas

**Files:**
- Modify: `src/zndraw_auth/db.py`
- Modify: `src/zndraw_auth/schemas.py`
- Modify: `src/zndraw_auth/__init__.py`
- Test: `tests/zndraw_auth/test_display_names.py` (already created)

- [ ] **Step 1: Add the column to `User`**

In `src/zndraw_auth/db.py`, inside `class User(SQLAlchemyBaseUserTableUUID, Base):` body, immediately after the `is_guest` block (around line 55), add:

```python
    if TYPE_CHECKING:  # pragma: no cover
        display_name: str
    else:
        display_name: Mapped[str] = mapped_column(
            sa.String(64),
            unique=True,
            index=True,
            nullable=False,
        )
```

(Mirror the `is_guest` TYPE_CHECKING pattern. Keep `nullable=False` — every row must have one.)

- [ ] **Step 2: Update `ensure_default_admin` to fill display_name**

In `src/zndraw_auth/db.py`, inside `ensure_default_admin`, where the `admin = User(...)` block is constructed (around line 182), add `display_name=…`:

```python
        from zndraw_auth.display_names import generate_unique_display_name

        display_name = await generate_unique_display_name(session)
        admin = User(
            email=settings.default_admin_email,
            hashed_password=hashed,
            is_active=True,
            is_superuser=True,
            is_verified=True,
            display_name=display_name,
        )
```

(Local import avoids a top-level cycle with the new module.)

- [ ] **Step 3: Update schemas**

Edit `src/zndraw_auth/schemas.py`:

```python
class UserRead(schemas.BaseUser[uuid.UUID]):
    """Schema for reading user data (responses).

    Overrides ``email`` to plain ``str`` so that addresses with reserved
    TLDs (``.local``, ``.test``) stored by ``ensure_default_admin`` can
    be serialized without a pydantic ``EmailStr`` validation error.
    """

    email: str  # type: ignore[assignment]
    is_guest: bool = False
    display_name: str


class UserCreate(schemas.BaseUserCreate):
    """Schema for creating a new user."""

    is_guest: bool = False
    display_name: str | None = None  # server fills when omitted
```

(`UserUpdate` stays unchanged — no rename endpoint in this task.)

- [ ] **Step 4: Export the helpers from `zndraw_auth`**

In `src/zndraw_auth/__init__.py`, after the `from zndraw_auth.db import (…)` block, add:

```python
from zndraw_auth.display_names import (
    DISPLAY_NAME_PATTERN,
    RESERVED_DISPLAY_NAMES,
    generate_unique_display_name,
    validate_display_name,
)
```

…and add `"DISPLAY_NAME_PATTERN"`, `"RESERVED_DISPLAY_NAMES"`, `"generate_unique_display_name"`, `"validate_display_name"` to `__all__` (alphabetical).

- [ ] **Step 5: Run the display-name tests again to confirm `User.display_name` works**

Run: `uv run pytest tests/zndraw_auth/test_display_names.py -v`

Expected: still PASS (User row construction with `display_name=` now works against the real column).

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_auth/db.py src/zndraw_auth/schemas.py src/zndraw_auth/__init__.py
git commit -m "feat(auth): User.display_name column + schema fields"
```

---

### Task 3: UserManager fills/validates display_name on create

**Files:**
- Modify: `src/zndraw_auth/users.py`
- Create: `tests/zndraw_auth/test_user_create.py`

- [ ] **Step 1: Write failing tests**

Create `tests/zndraw_auth/test_user_create.py`:

```python
"""End-to-end UserManager.create coverage for display_name handling."""

from __future__ import annotations

import pytest
from httpx import AsyncClient
from sqlmodel import select
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth.db import User
from zndraw_auth.display_names import DISPLAY_NAME_PATTERN


@pytest.mark.asyncio
async def test_register_with_explicit_display_name(
    client: AsyncClient, session: AsyncSession
) -> None:
    resp = await client.post(
        "/auth/register",
        json={
            "email": "alice@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "alice-the-explorer",
        },
    )
    assert resp.status_code == 201, resp.text
    body = resp.json()
    assert body["display_name"] == "alice-the-explorer"


@pytest.mark.asyncio
async def test_register_without_display_name_fills_one(
    client: AsyncClient, session: AsyncSession
) -> None:
    resp = await client.post(
        "/auth/register",
        json={"email": "bob@example.com", "password": "very-strong-passw0rd"},
    )
    assert resp.status_code == 201, resp.text
    body = resp.json()
    name = body["display_name"]
    assert DISPLAY_NAME_PATTERN.fullmatch(name), name

    row = (
        await session.exec(select(User).where(User.email == "bob@example.com"))
    ).one()
    assert row.display_name == name


@pytest.mark.asyncio
async def test_register_duplicate_display_name_returns_409(
    client: AsyncClient, session: AsyncSession
) -> None:
    await client.post(
        "/auth/register",
        json={
            "email": "first@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "shared-name-here",
        },
    )
    resp = await client.post(
        "/auth/register",
        json={
            "email": "second@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "shared-name-here",
        },
    )
    assert resp.status_code == 409, resp.text
    body = resp.json()
    assert body["type"].endswith("/username-exists")


@pytest.mark.asyncio
async def test_register_malformed_display_name_returns_422(
    client: AsyncClient,
) -> None:
    resp = await client.post(
        "/auth/register",
        json={
            "email": "carol@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "BadName!",
        },
    )
    assert resp.status_code == 422, resp.text


@pytest.mark.asyncio
async def test_register_reserved_display_name_returns_422(
    client: AsyncClient,
) -> None:
    resp = await client.post(
        "/auth/register",
        json={
            "email": "dave@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "admin",
        },
    )
    assert resp.status_code == 422, resp.text
```

- [ ] **Step 2: Run tests to confirm they fail**

Run: `uv run pytest tests/zndraw_auth/test_user_create.py -v`

Expected: all five fail — UserManager hasn't been wired yet, so display_name is rejected by `User.display_name NOT NULL`.

- [ ] **Step 3: Override `UserManager.create`**

In `src/zndraw_auth/users.py`, inside `class UserManager(...):` add the override after `is_dev_mode = False`:

```python
    async def create(
        self,
        user_create,  # type: ignore[no-untyped-def]
        safe: bool = False,
        request=None,
    ):
        from sqlalchemy.exc import IntegrityError

        from zndraw.exceptions import UsernameExists
        from zndraw_auth.display_names import (
            generate_unique_display_name,
            validate_display_name,
        )

        if user_create.display_name is None:
            user_create.display_name = await generate_unique_display_name(
                self.user_db.session  # type: ignore[attr-defined]
            )
        else:
            validate_display_name(user_create.display_name)

        try:
            return await super().create(user_create, safe=safe, request=request)
        except IntegrityError as exc:
            raise UsernameExists.exception(
                f"Display name '{user_create.display_name}' is already taken"
            ) from exc
```

(`UserManager.user_db` is the `SQLAlchemyUserDatabase` bound at construction; its `.session` attribute is the active AsyncSession — type-ignored to avoid private-API friction. Local imports avoid cycles.)

- [ ] **Step 4: Run tests to confirm they pass**

Run: `uv run pytest tests/zndraw_auth/test_user_create.py -v`

Expected: all five PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_auth/users.py tests/zndraw_auth/test_user_create.py
git commit -m "feat(auth): UserManager.create fills+validates display_name"
```

---

### Task 4: `/v1/users/available-display-name` suggestion endpoint

**Files:**
- Create: `src/zndraw/routes/users.py`
- Modify: `src/zndraw/app.py`
- Test: extend `tests/zndraw/test_auth_endpoints.py`

- [ ] **Step 1: Write a failing test**

Append to `tests/zndraw/test_auth_endpoints.py`:

```python
@pytest.mark.asyncio
async def test_available_display_name_returns_valid_slug(
    client: AsyncClient,
) -> None:
    resp = await client.get("/v1/users/available-display-name")
    assert resp.status_code == 200, resp.text
    body = resp.json()
    from zndraw_auth.display_names import DISPLAY_NAME_PATTERN
    assert DISPLAY_NAME_PATTERN.fullmatch(body["display_name"])


@pytest.mark.asyncio
async def test_available_display_name_never_collides_with_existing(
    client: AsyncClient, session: AsyncSession
) -> None:
    # Register one user, then ask for a suggestion — must differ.
    await client.post(
        "/v1/auth/register",
        json={
            "email": "eve@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "eve-the-curious",
        },
    )
    resp = await client.get("/v1/users/available-display-name")
    assert resp.json()["display_name"] != "eve-the-curious"
```

- [ ] **Step 2: Run tests to confirm 404**

Run: `uv run pytest tests/zndraw/test_auth_endpoints.py -v -k available_display_name`

Expected: both fail with 404 (route missing).

- [ ] **Step 3: Create the router**

Create `src/zndraw/routes/users.py`:

```python
"""User-facing utility endpoints (display-name suggestions)."""

from __future__ import annotations

from fastapi import APIRouter
from pydantic import BaseModel

from zndraw.exceptions import problem_responses
from zndraw_auth import SessionDep, generate_unique_display_name

router = APIRouter(prefix="/v1/users", tags=["users"])


class DisplayNameSuggestion(BaseModel):
    """Server-suggested display name for the registration form."""

    display_name: str


@router.get("/available-display-name", responses=problem_responses())
async def get_available_display_name(
    session: SessionDep,
) -> DisplayNameSuggestion:
    return DisplayNameSuggestion(
        display_name=await generate_unique_display_name(session)
    )
```

- [ ] **Step 4: Wire it into the app**

Edit `src/zndraw/app.py`. Add the import next to the other route imports (after `from zndraw.routes.trajectory import router as trajectory_router`):

```python
from zndraw.routes.users import router as users_router
```

And register it next to the other `app.include_router(...)` calls (after `utility_router`):

```python
app.include_router(users_router)
```

- [ ] **Step 5: Run tests to confirm pass**

Run: `uv run pytest tests/zndraw/test_auth_endpoints.py -v -k available_display_name`

Expected: both PASS.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/routes/users.py src/zndraw/app.py tests/zndraw/test_auth_endpoints.py
git commit -m "feat(api): GET /v1/users/available-display-name"
```

---

### Task 5: Guest flow returns `GuestSessionResponse` with display_name

**Files:**
- Modify: `src/zndraw/routes/auth.py`
- Test: extend `tests/zndraw/test_e2e_guest_auth.py`

- [ ] **Step 1: Write failing test**

In `tests/zndraw/test_e2e_guest_auth.py`, add (or extend an existing fixture/test):

```python
@pytest.mark.asyncio
async def test_guest_session_includes_display_name(client: AsyncClient) -> None:
    resp = await client.post("/v1/auth/guest")
    assert resp.status_code == 200, resp.text
    body = resp.json()
    from zndraw_auth.display_names import DISPLAY_NAME_PATTERN
    assert DISPLAY_NAME_PATTERN.fullmatch(body["display_name"])
    assert body["email"].endswith("@guest.user")
    assert body["token_type"] == "bearer"
    assert isinstance(body["access_token"], str) and body["access_token"]
```

- [ ] **Step 2: Confirm it fails**

Run: `uv run pytest tests/zndraw/test_e2e_guest_auth.py::test_guest_session_includes_display_name -v`

Expected: KeyError / KeyError-style failure ("display_name" not in response).

- [ ] **Step 3: Implement the typed response**

In `src/zndraw/routes/auth.py`, add a Pydantic model and replace the `-> dict` return type:

```python
from typing import Literal

from pydantic import BaseModel
from zndraw_auth import generate_unique_display_name


class GuestSessionResponse(BaseModel):
    access_token: str
    token_type: Literal["bearer"] = "bearer"
    email: str
    display_name: str


@router.post("/guest")
async def create_guest_session(
    auth_settings: AuthSettingsDep,
    user_manager: Annotated[UserManager, Depends(get_user_manager)],
    settings: Annotated[Settings, Depends(get_zndraw_settings)],
    session: SessionDep,
) -> GuestSessionResponse:
    """Create anonymous guest user (is_guest=True) and return JWT token."""
    email = f"{uuid4().hex[:8]}@guest.user"
    password = settings.guest_password.get_secret_value()
    display_name = await generate_unique_display_name(session)

    user = await user_manager.create(
        UserCreate(
            email=email,
            password=password,
            is_guest=True,
            display_name=display_name,
        )
    )

    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    token = await strategy.write_token(user)

    return GuestSessionResponse(
        access_token=token, email=email, display_name=user.display_name
    )
```

Also add `from zndraw_auth import SessionDep` to the existing imports if not already present.

- [ ] **Step 4: Confirm it passes**

Run: `uv run pytest tests/zndraw/test_e2e_guest_auth.py -v`

Expected: new test PASS; existing tests still PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/auth.py tests/zndraw/test_e2e_guest_auth.py
git commit -m "feat(auth): GuestSessionResponse carries display_name"
```

---

### Task 6: `resolve_owner` returns display_name; new path-segment resolver

**Files:**
- Modify: `src/zndraw/dependencies.py`
- Modify: `src/zndraw/models.py`
- Test: extend `tests/zndraw/test_owner_resolution.py`

- [ ] **Step 1: Write failing tests**

In `tests/zndraw/test_owner_resolution.py` (extend the existing file; if a function-name list exists, append):

```python
@pytest.mark.asyncio
async def test_resolve_owner_returns_display_name_for_user(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import OwnerKind, resolve_owner
    from zndraw_auth.db import User

    user = User(
        email="ada@example.com",
        hashed_password="x",
        display_name="ada-lovelace-coder",
    )
    session.add(user)
    await session.commit()
    await session.refresh(user)

    kind, label = await resolve_owner(session, user.id)
    assert kind == OwnerKind.USER
    assert label == "ada-lovelace-coder"


@pytest.mark.asyncio
async def test_resolve_owner_returns_group_name_unchanged(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import OwnerKind, resolve_owner
    from zndraw.models import Group
    from zndraw_auth.db import User
    import uuid

    creator = User(
        email="creator@example.com",
        hashed_password="x",
        display_name="creator-of-groups",
    )
    session.add(creator)
    await session.commit()
    await session.refresh(creator)

    group = Group(name="my-group", created_by_id=creator.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    kind, label = await resolve_owner(session, group.id)
    assert kind == OwnerKind.GROUP
    assert label == "my-group"


@pytest.mark.asyncio
async def test_get_owner_uuid_from_segment_resolves_user(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import get_owner_uuid_from_segment
    from zndraw_auth.db import User

    user = User(
        email="seg@example.com",
        hashed_password="x",
        display_name="seg-test-user",
    )
    session.add(user)
    await session.commit()
    await session.refresh(user)

    assert (
        await get_owner_uuid_from_segment(session, "seg-test-user")
    ) == user.id


@pytest.mark.asyncio
async def test_get_owner_uuid_from_segment_resolves_group(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import get_owner_uuid_from_segment
    from zndraw.models import Group
    from zndraw_auth.db import User

    creator = User(
        email="gcreator@example.com",
        hashed_password="x",
        display_name="g-creator-display",
    )
    session.add(creator)
    await session.commit()
    await session.refresh(creator)

    group = Group(name="visible-group", created_by_id=creator.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    assert (
        await get_owner_uuid_from_segment(session, "visible-group")
    ) == group.id


@pytest.mark.asyncio
async def test_get_owner_uuid_from_segment_unknown_raises_user_not_found(
    session: AsyncSession,
) -> None:
    from zndraw.dependencies import get_owner_uuid_from_segment
    from zndraw.exceptions import ProblemError

    with pytest.raises(ProblemError) as excinfo:
        await get_owner_uuid_from_segment(session, "no-such-owner-here")
    assert excinfo.value.problem.status == 404
```

- [ ] **Step 2: Run tests to confirm they fail**

Run: `uv run pytest tests/zndraw/test_owner_resolution.py -v`

Expected: the two existing tests pass (or whatever the baseline is), the five new tests fail at import or assertion.

- [ ] **Step 3: Update `resolve_owner` to return display_name**

In `src/zndraw/dependencies.py`, replace the existing `resolve_owner`:

```python
async def resolve_owner(
    session: AsyncSession, owner_id: UUID
) -> tuple[OwnerKind, str] | None:
    """Look up ``owner_id`` as a user (returns display_name) or group (returns name)."""
    user = await session.get(User, owner_id)
    if user is not None:
        return OwnerKind.USER, user.display_name
    group = await session.get(Group, owner_id)
    if group is not None:
        return OwnerKind.GROUP, group.name
    return None
```

- [ ] **Step 4: Add `get_owner_uuid_from_segment`**

In `src/zndraw/dependencies.py` (somewhere after `resolve_owner`), add:

```python
async def get_owner_uuid_from_segment(
    session: AsyncSession, owner: str
) -> UUID:
    """Resolve a path display-name segment to a user UUID, or a group UUID by name.

    Tries ``User.display_name`` first, then ``Group.name`` (lowercase groups
    that happen to satisfy the path regex still work). Raises ``UserNotFound``
    when neither match — the path regex already gated malformed input.
    """
    from zndraw.exceptions import UserNotFound

    user_id = await session.scalar(
        select(User.id).where(User.display_name == owner).limit(1)
    )
    if user_id is not None:
        return user_id
    group_id = await session.scalar(
        select(Group.id).where(Group.name == owner).limit(1)
    )
    if group_id is not None:
        return group_id
    raise UserNotFound.exception(f"Owner '{owner}' not found")
```

(The existing `from zndraw.models import Group, GroupMembership, Room, RoomGeometry, RoomShareLink` import already imports `Group`. Add `from zndraw_auth import User` if it's not already in scope — it is, via `from zndraw_auth import ( SessionDep, User, ... )`.)

- [ ] **Step 5: Run tests to confirm pass**

Run: `uv run pytest tests/zndraw/test_owner_resolution.py -v`

Expected: all owner_resolution tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/dependencies.py tests/zndraw/test_owner_resolution.py
git commit -m "feat(api): resolve_owner returns display_name; add path-segment resolver"
```

---

### Task 7: Flip two-segment path params from UUID → display-name string

**Files:**
- Modify: `src/zndraw/dependencies.py`
- Modify: `src/zndraw/routes/rooms.py`
- Modify: `src/zndraw/routes/bookmarks.py`
- Modify: `src/zndraw/routes/chat.py`
- Modify: `src/zndraw/routes/edit_lock.py`
- Modify: `src/zndraw/routes/figures.py`
- Modify: `src/zndraw/routes/frames.py`
- Modify: `src/zndraw/routes/geometries.py`
- Modify: `src/zndraw/routes/isosurface.py`
- Modify: `src/zndraw/routes/presets.py`
- Modify: `src/zndraw/routes/progress.py`
- Modify: `src/zndraw/routes/screenshots.py`
- Modify: `src/zndraw/routes/selection_groups.py`
- Modify: `src/zndraw/routes/share_links.py`
- Modify: `src/zndraw/routes/step.py`
- Modify: `src/zndraw/routes/trajectory.py`

This is the largest mechanical change. Each file follows the same pattern:
- Path param `owner_id: UUID = Path()` → `owner: Annotated[str, Path(pattern=r"^[a-z][a-z0-9-]{2,63}$")]`.
- Lookups via `_load_room_by_address(session, owner_id, room_name)` now go through a new helper that first resolves display_name → UUID.

- [ ] **Step 1: Add a `_load_room_by_segment` helper to `dependencies.py`**

Right above the existing `_load_room_by_address`, add:

```python
async def _load_room_by_segment(
    session: AsyncSession, owner: str, room_name: str
) -> Room | None:
    """Resolve owner display_name → UUID, then look up the room."""
    owner_id = await get_owner_uuid_from_segment(session, owner)
    return await _load_room_by_address(session, owner_id, room_name)
```

- [ ] **Step 2: Update the central access deps**

Inside `src/zndraw/dependencies.py`, replace the existing definitions. Path params use the `Annotated[str, Path(pattern=...)]` form (no Python default; FastAPI binds them by position):

```python
async def get_share_context_two_segment(
    session: SessionDep,
    owner: Annotated[str, Path(pattern=r"^[a-z][a-z0-9-]{2,63}$")],
    room_name: str = Path(),
    x_room_share_token: str | None = Header(default=None, alias="X-Room-Share-Token"),
) -> ShareContext | None:
    """Resolve the share-token header for the two-segment path."""
    if x_room_share_token is None:
        return None
    room = await _load_room_by_segment(session, owner, room_name)
    if room is None:
        return None
    return await resolve_share_token(session, x_room_share_token, room.id)


async def _load_access_context(
    session: AsyncSession,
    owner: str,
    room_name: str,
    current_user: User,
    share: ShareContext | None,
) -> AccessContext:
    room = await _load_room_by_segment(session, owner, room_name)
    if room is None:
        raise RoomNotFound.exception(f"Room {owner}/{room_name} not found")
    group_role: GroupRole | None = None
    if room.owner_group_id is not None:
        group_role = await fetch_group_role(
            session, current_user.id, room.owner_group_id
        )
    return AccessContext(room=room, share=share, group_role=group_role)


async def get_readable_room(
    session: SessionDep,
    current_user: CurrentUserDep,
    share: TwoSegmentShareTokenDep,
    owner: Annotated[str, Path(pattern=r"^[a-z][a-z0-9-]{2,63}$")],
    room_name: str = Path(),
) -> AccessContext:
    ctx = await _load_access_context(session, owner, room_name, current_user, share)
    if not can_read(current_user, ctx.room, ctx.share, group_role=ctx.group_role):
        raise RoomNotFound.exception(f"Room {owner}/{room_name} not found")
    return ctx
```

- [ ] **Step 3: Update `get_verified_session_id` and `get_active_session_cam_id`**

Replace `owner_id: UUID = Path()` with `owner: Annotated[str, Path(pattern=r"^[a-z][a-z0-9-]{2,63}$")]`, and replace `_load_room_by_address(session, owner_id, room_name)` with `_load_room_by_segment(session, owner, room_name)`. Update the error message to use `{owner}/{room_name}` (no longer a UUID).

- [ ] **Step 4: Update `verify_room` composed-address parser**

`src/zndraw/dependencies.py:verify_room` currently parses `"<owner_uuid>/<name>"` and calls `UUID(owner_str)`. Replace with:

```python
async def verify_room(session: AsyncSession, room_id: str) -> Room:
    """Verify room exists and return it, or raise RoomNotFound.

    Accepts either the surrogate UUID primary key or a composed
    ``<owner_display>/<room_name>`` address.
    """
    if "/" in room_id:
        owner_str, _, name_part = room_id.partition("/")
        room = await _load_room_by_segment(session, owner_str, name_part)
    else:
        room = await session.get(Room, room_id)
    if room is None:
        raise RoomNotFound.exception(f"Room with id {room_id} not found")
    return room
```

- [ ] **Step 5: Update `get_writable_room_id`**

In `src/zndraw/dependencies.py:get_writable_room_id`, replace the UUID parse:

```python
async def get_writable_room_id(
    request: Request,
    session: SessionDep,
    current_user: CurrentUserDep,
    redis: RedisDep,
    room_id: str = Path(),
    x_room_share_token: str | None = Header(default=None, alias="X-Room-Share-Token"),
) -> str:
    """Verify a room is writable and return the surrogate UUID string."""
    validate_room_id(room_id)
    if room_id in ("@global", "@internal"):
        return room_id
    owner_part, _, name_part = room_id.partition("/")
    room = await _load_room_by_segment(session, owner_part, name_part)
    if room is None:
        raise RoomNotFound.exception(f"Room {room_id} not found")
    ...
```

- [ ] **Step 6: Update each two-segment route file**

For every file in the list above, look for path declarations using `{owner_id}` and signatures using `owner_id: UUID = Path()`. Replace the noqa comments and decorator paths from `/{owner_id}/{room_name}` to `/{owner}/{room_name}`, and update parameter signatures to `owner: Annotated[str, Path(pattern=r"^[a-z][a-z0-9-]{2,63}$")]`. Where the handler directly calls `_load_room_by_address`, replace with `_load_room_by_segment`. Most handlers only use `AccessReadDep`/`AccessEditDep`/`AccessManageDep`, which are now display-name-aware after Step 2.

Concrete edits in `src/zndraw/routes/rooms.py`:

- `_load_room_by_address(session, request.owner_id, name)` (in `create_room`, around the existence check) — see Task 8 for the schema change.
- The path `{owner_id}/{room_name}` (router decorators around lines 542, 569, 600, 636) becomes `{owner}/{room_name}`.

Concrete edits in `src/zndraw/routes/chat.py`:

- Prefix update: `prefix="/v1/rooms/{owner}/{room_name}/chat/messages"`.

Concrete edits in `src/zndraw/routes/share_links.py` — search for `/{owner_id}/` and replace each occurrence.

For every other file: same pattern — search-and-replace `{owner_id}` → `{owner}` in path strings and `owner_id: UUID = Path()` → `owner: Annotated[str, Path(pattern=…)]` in signatures.

- [ ] **Step 7: Type-check Python**

Run: `uv run python -c "import zndraw.app"`

Expected: imports without error.

- [ ] **Step 8: Run the access-dep test suite**

Run: `uv run pytest tests/zndraw/test_owner_resolution.py tests/zndraw/test_routes_rooms.py -v`

Expected: known failures while specs (room_id wire contract, owner_id body field) are still UUID-shaped — these are addressed in Tasks 8–10. Continue.

- [ ] **Step 9: Commit (work-in-progress)**

```bash
git add src/zndraw/dependencies.py src/zndraw/routes/
git commit -m "refactor(api): two-segment path owner is now a display-name string"
```

---

### Task 8: Rename body fields `owner_id` → `owner`, `new_owner_id` → `new_owner`

**Files:**
- Modify: `src/zndraw/schemas.py`
- Modify: `src/zndraw/routes/rooms.py`
- Test: `tests/zndraw/test_routes_rooms.py` (extend / adjust fixtures — Task 12)

- [ ] **Step 1: Update body schemas**

In `src/zndraw/schemas.py`:

```python
class RoomCreate(BaseModel):
    """Request body for POST /v1/rooms."""

    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    name: str = Field(pattern=r"^[a-zA-Z0-9\-_]+$", min_length=1, max_length=128)
    description: str | None = None
    copy_from: str | None = None
    visibility: Visibility | None = None


class RoomResponse(BaseModel):
    """Response body for room details — matches frontend Room interface."""

    room_id: str  # composed: {owner}/{room_name}
    id: str
    description: str | None = None
    frame_count: int = 0
    visibility: Visibility = Visibility.PUBLIC
    owner: str  # display_name (user) OR group name
    owner_kind: Literal["user", "group"]
    owner_label: str
    is_default: bool = False
    metadata: dict[str, str] | None = None

    model_config = ConfigDict(from_attributes=True)


class RoomPatchRequest(BaseModel):
    """Request body for PATCH /v1/rooms/{owner}/{room_name}."""

    description: str | None = None
    frame_count: int | None = Field(None, ge=0)
    visibility: Visibility | None = None
    new_owner: str | None = Field(default=None, pattern=r"^[a-z][a-z0-9-]{2,63}$")
```

(`owner_label` stays — for backwards-readable response data; same value as `owner` for users.)

- [ ] **Step 2: Update `create_room` to resolve the new body field**

In `src/zndraw/routes/rooms.py:create_room`, replace every reference to `request.owner_id` with the resolution:

```python
from zndraw.dependencies import get_owner_uuid_from_segment
...
owner_id = await get_owner_uuid_from_segment(session, request.owner)
```

Then use `owner_id` for the existing `resolve_owner` / `_load_room_by_address` calls instead of `request.owner_id`. Replace `if request.owner_id == current_user.id or ...` with `if owner_id == current_user.id or ...`.

- [ ] **Step 3: Update `update_room`**

In the same file, `update_room`:

```python
if updates.new_owner is not None:
    new_owner_id = await get_owner_uuid_from_segment(session, updates.new_owner)
    resolved = await resolve_owner(session, new_owner_id)
    ...
```

Then everywhere `updates.new_owner_id` was used, substitute `new_owner_id`.

- [ ] **Step 4: Update `_resolve_room_owner` return label**

`_resolve_room_owner` already returns `owner_label` from `resolve_owner`, which now returns display_name. Update the response composition where the field was `owner_id=…`: change to `owner=label_string` (the display_name or group name). Drop the UUID. **Note**: many downstream callers (frontend, tests) currently read `owner_id`; once we change to `owner`, those need updating in Tasks 11–12.

```python
async def _resolve_room_owner(
    session: AsyncSession, room: Room
) -> tuple[str, Literal["user", "group"], str]:
    """Return ``(owner_label, owner_kind, owner_label)``.

    The label and the new ``owner`` field are now identical (display_name for users,
    group name for groups). They remain two separate fields for response stability.
    """
    owner_uuid = room.owner_user_id or room.owner_group_id
    if owner_uuid is None:
        raise RuntimeError(...)  # same as before
    resolved = await resolve_owner(session, owner_uuid)
    if resolved is None:
        return "", "user", ""
    kind, label = resolved
    return label, "group" if kind is OwnerKind.GROUP else "user", label
```

Update every callsite to unpack `(owner_label, owner_kind, _label)` or just `(owner, kind, label)` and pass `owner=owner`, `owner_kind=kind`, `owner_label=label` into `RoomResponse` / `RoomUpdate` / `RoomCreateResponse`. **`RoomCreateResponse.room_id`** still uses `room.public_address` (changes in Task 10).

- [ ] **Step 5: Skip tests for now; revisit in Task 12**

Body schema is mid-flight. The test-suite update happens after Task 10 / 11.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/schemas.py src/zndraw/routes/rooms.py
git commit -m "feat(api): RoomCreate.owner / RoomPatchRequest.new_owner (display-name strings)"
```

---

### Task 9: Email → display_name in REST response payloads

**Files:**
- Modify: `src/zndraw/schemas.py`
- Modify: `src/zndraw/routes/chat.py`
- Modify: `src/zndraw/routes/rooms.py` (presence + sessions)
- Modify: `src/zndraw/routes/groups.py` (member listing)

- [ ] **Step 1: Replace `email` with `display_name` in response schemas**

In `src/zndraw/schemas.py`:

```python
class MessageResponse(BaseModel):
    id: int
    room_id: str
    user_id: UUID
    content: str
    created_at: datetime
    updated_at: datetime | None = None
    display_name: str | None = None

    model_config = ConfigDict(from_attributes=True)


class PresenceSessionResponse(BaseModel):
    sid: str
    user_id: UUID
    display_name: str | None


class SessionItem(BaseModel):
    sid: str
    display_name: str
    camera_key: str


class GroupMemberResponse(BaseModel):
    user_id: UUID
    display_name: str | None
    role: GroupRole
    joined_at: datetime

    model_config = ConfigDict(from_attributes=True)
```

- [ ] **Step 2: Update chat route**

In `src/zndraw/routes/chat.py`:
- Rename `_message_to_response(msg, email=None)` to `_message_to_response(msg, display_name=None)` and pass `display_name=display_name` into `MessageResponse`.
- In `list_messages`, replace `email_map: dict[str, str | None] = {}` lookups with `display_name_map`; pull `user.display_name`.
- In `create_message`, replace `email = current_user.email` with `display_name = current_user.display_name`; pass `display_name=display_name` into both `MessageNew.for_room(...)` and `_message_to_response(...)`.
- In `edit_message`, return `_message_to_response(msg, current_user.display_name)`.

- [ ] **Step 3: Update room presence + sessions**

In `src/zndraw/routes/rooms.py:get_room_presence`, the Redis-stored hash currently caches `email`. After Task 13 the hash will store `display_name`. For now:

```python
sessions_list.append(
    PresenceSessionResponse(
        sid=entry["sid"],
        user_id=_UUID(camera.owner),
        display_name=entry.get("display_name") or entry.get("email"),  # transition guard
    )
)
```

The `or entry.get("email")` fallback is dropped at the end of Task 13 once the Redis payload writes `display_name`.

In `list_sessions`, replace `email` filter param with `display_name`; pivot `entry_email` → `entry_display_name = entry.get("display_name", entry.get("email", ""))`; pass `display_name=` into `SessionItem`.

- [ ] **Step 4: Update group member listing**

In `src/zndraw/routes/groups.py` (search for `email=`), wherever a `GroupMemberResponse` is constructed, replace `email=user.email` with `display_name=user.display_name`.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/schemas.py src/zndraw/routes/chat.py src/zndraw/routes/rooms.py src/zndraw/routes/groups.py
git commit -m "feat(api): REST response payloads carry display_name (no more email)"
```

---

### Task 10: `Room.public_address` returns display-name composed string

**Files:**
- Modify: `src/zndraw/models.py`
- Modify: `src/zndraw/routes/rooms.py`
- Modify: `src/zndraw/socket_events.py`
- Modify: `src/zndraw/socketio.py`
- Modify: `src/zndraw/broadcast.py`
- Modify: anywhere `room.public_address` is read

- [ ] **Step 1: Add an async helper**

In `src/zndraw/models.py`, after the `Room` class definition, add a module-level helper:

```python
async def build_public_address(session, room: Room) -> str:
    """Return ``{display_name}/{room_name}`` for the room's owner.

    For users → ``User.display_name``; for groups → ``Group.name``.
    Falls back to the owner UUID if the row is missing (orphaned FK).
    """
    from sqlmodel import select
    from zndraw_auth.db import User

    if room.owner_user_id is not None:
        label = await session.scalar(
            select(User.display_name).where(User.id == room.owner_user_id).limit(1)
        )
        return f"{label or room.owner_user_id}/{room.room_name}"
    if room.owner_group_id is not None:
        label = await session.scalar(
            select(Group.name).where(Group.id == room.owner_group_id).limit(1)
        )
        return f"{label or room.owner_group_id}/{room.room_name}"
    return room.room_name
```

Leave `Room.public_address` (the existing synchronous UUID-form property) in place — it remains useful for `verify_room` internal paths and for tests asserting the UUID form. The helper is the canonical "what the user sees" address.

- [ ] **Step 2: Replace every `room.public_address` user-facing emission**

For each spot that constructs a `RoomResponse`, `RoomUpdate`, `RoomCreateResponse`, `RoomPatchResponse`, `RoomLeaveResponse`, or `RoomRenamed` — and for the `RoomScopedEvent.for_room(room, ...)` helper — switch to using the new helper.

Cleanest approach: bake the lookup into the event helper.

In `src/zndraw/socket_events.py`, change `RoomScopedEvent.for_room` to take an explicit address:

```python
class RoomScopedEvent(BaseModel):
    room_id: UUID
    room_address: str

    @classmethod
    def for_room(cls, room: Room, /, *, room_address: str | None = None, **kwargs) -> Self:
        return cls(
            room_id=UUID(room.id),
            room_address=room_address or room.public_address,
            **kwargs,
        )
```

In every callsite (e.g., `chat.py:create_message`, `rooms.py:broadcast_room_update`, `socketio.py:room_join`/`typing_*`), compute `room_address = await build_public_address(session, room)` ahead of the `for_room(...)` call and pass it via `room_address=`.

`broadcast_to_room` also uses `room.public_address` for room-feed dispatching; switch the same way.

In `RoomResponse`/`RoomCreateResponse`/`RoomPatchResponse`/`RoomUpdate` construction (in `routes/rooms.py`):

```python
public_address = await build_public_address(session, room)
return RoomResponse(
    room_id=public_address,
    id=room.id,
    ...
    owner=owner_label,   # already set in Task 8
    ...
)
```

- [ ] **Step 3: Confirm protected schema test still passes**

Run: `uv run pytest tests/zndraw/test_event_schema_audit.py -v`

Expected: PASS (the audit asserts presence of `room_id`+`room_address`; not the values).

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/models.py src/zndraw/socket_events.py src/zndraw/socketio.py src/zndraw/routes/ src/zndraw/broadcast.py
git commit -m "feat(api): public_address renders display_name owner segment"
```

---

### Task 11: Socket event schemas use `owner` + `display_name`

**Files:**
- Modify: `src/zndraw/socket_events.py`
- Modify: `src/zndraw/socketio.py`

- [ ] **Step 1: Schema renames**

In `src/zndraw/socket_events.py`:

```python
class RoomJoin(BaseModel):
    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    room_name: str
    client_type: Literal["frontend", "pyclient"] = "frontend"


class RoomLeave(BaseModel):
    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    room_name: str


class TypingStart(BaseModel):
    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    room_name: str


class TypingStop(BaseModel):
    owner: str = Field(pattern=r"^[a-z][a-z0-9-]{2,63}$")
    room_name: str


class UserGetResponse(BaseModel):
    id: UUID
    display_name: str
    is_superuser: bool


class SessionJoined(RoomScopedEvent):
    user_id: UUID
    sid: str
    display_name: str | None = None


class MessageNew(RoomScopedEvent):
    id: int
    user_id: UUID
    content: str
    created_at: datetime
    updated_at: datetime | None = None
    display_name: str | None = None


class Typing(RoomScopedEvent):
    user_id: UUID
    display_name: str | None = None
    is_typing: bool
```

(Field added: `from pydantic import BaseModel, Field, model_validator` — append `Field` to the existing import.)

- [ ] **Step 2: Update handlers in `src/zndraw/socketio.py`**

Replace `data.owner_id` with `data.owner` everywhere. After lookup, resolve UUID first:

```python
from zndraw.dependencies import get_owner_uuid_from_segment, _load_room_by_segment

# in room_join:
room = await _load_room_by_segment(session, data.owner, data.room_name)
if room is None:
    raise RoomNotFound.exception(f"Room {data.owner}/{data.room_name} not found")
```

In `user_get`, change `UserGetResponse(id=user.id, email=user.email, ...)` → `UserGetResponse(id=user.id, display_name=user.display_name, ...)`.

In `room_join`, replace:

```python
email = user.email
...
camera_key = f"cam:{email}:{sid[:8]}"
camera_value = json.dumps(
    {"sid": sid, "email": email, "data": camera.model_dump()}
)
```

with:

```python
display_name = user.display_name
...
camera_key = f"cam:{display_name}:{sid[:8]}"
camera_value = json.dumps(
    {"sid": sid, "display_name": display_name, "data": camera.model_dump()}
)
```

And in `SessionJoined.for_room(...)` pass `display_name=display_name` instead of `email=email`.

In `_handle_typing`, replace `email = user.email if user else None` with `display_name = user.display_name if user else None`; pass `display_name=` into `Typing.for_room`.

- [ ] **Step 3: Drop the email fallback in presence/sessions**

Now that the camera hash stores `display_name`, remove the `or entry.get("email")` transition guards added in Task 9. Update the entry-read code in `routes/rooms.py:get_room_presence` and `list_sessions` to read `entry["display_name"]`.

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/socket_events.py src/zndraw/socketio.py src/zndraw/routes/rooms.py
git commit -m "feat(socket): owner segment + display_name in socket payloads"
```

---

### Task 12: Extend backend tests for the new contract

**Files:**
- Modify: `tests/zndraw/test_routes_rooms.py`
- Modify: any other `tests/zndraw/test_routes_*.py` that constructs owner URLs
- Modify: `tests/zndraw/test_socket_commands.py`
- Modify: `tests/zndraw/test_cli.py` (covered in Task 14, but record here)

- [ ] **Step 1: Update room-route tests**

Search the file for `"{owner_id}/{room_name}"` formatting, `owner_id=user.id`, and `request.owner_id` patterns. Replace:

- POST `/v1/rooms` body: `{"owner_id": str(user.id), …}` → `{"owner": user.display_name, …}`.
- Path construction: `f"/v1/rooms/{user.id}/{name}"` → `f"/v1/rooms/{user.display_name}/{name}"`.
- Response asserts on `body["owner_id"]` → `body["owner"]` (equals the display_name).

- [ ] **Step 2: Add new tests for the regex gate**

In `tests/zndraw/test_routes_rooms.py`, add:

```python
@pytest.mark.asyncio
async def test_owner_segment_uuid_returns_422(
    auth_client: AsyncClient, user_with_display_name
) -> None:
    fake_uuid = "00000000-0000-0000-0000-000000000000"
    resp = await auth_client.get(f"/v1/rooms/{fake_uuid}/any-room")
    assert resp.status_code == 422


@pytest.mark.asyncio
async def test_owner_segment_unknown_display_name_returns_404(
    auth_client: AsyncClient,
) -> None:
    resp = await auth_client.get("/v1/rooms/no-such-display-name/any-room")
    assert resp.status_code == 404
```

- [ ] **Step 3: Update socket-command tests**

In `tests/zndraw/test_socket_commands.py`, find `RoomJoin(owner_id=…)`, `MessageNew(email=…)`, `Typing(email=…)`, `SessionJoined(email=…)`, `UserGetResponse(email=…)`. Replace `owner_id` → `owner` (string), and `email` → `display_name`. Update payload asserts the same way.

- [ ] **Step 4: Run the full test suite**

Run: `uv run pytest -x -q tests/zndraw/`

Expected: green. **If a protected test fails, stop immediately and surface it to the user — do not silently edit a `@pytest.mark.protected` test.**

- [ ] **Step 5: Commit**

```bash
git add tests/zndraw/
git commit -m "test: switch room-route and socket fixtures to display-name owners"
```

---

### Task 13: CLI uses display_name throughout

**Files:**
- Modify: `src/zndraw/cli.py`
- Modify: `src/zndraw/cli_agent/auth.py` (if it parses owner segments — verify)
- Modify: `tests/zndraw/test_cli.py`

- [ ] **Step 1: Write a failing test**

In `tests/zndraw/test_cli.py`, replace existing UUID validation tests and add:

```python
def test_validate_room_arg_accepts_display_name() -> None:
    from zndraw.cli import _validate_room_arg

    _validate_room_arg("happy-blue-rabbit/my-room", owner_name="happy-blue-rabbit")


def test_validate_room_arg_rejects_uuid_owner() -> None:
    import typer
    from zndraw.cli import _validate_room_arg

    with pytest.raises(typer.BadParameter):
        _validate_room_arg(
            "00000000-0000-0000-0000-000000000000/x",
            owner_name="happy-blue-rabbit",
        )


def test_validate_room_arg_rejects_unprefixed() -> None:
    import typer
    from zndraw.cli import _validate_room_arg

    with pytest.raises(typer.BadParameter):
        _validate_room_arg("just-a-name", owner_name="happy-blue-rabbit")
```

- [ ] **Step 2: Confirm they fail**

Run: `uv run pytest tests/zndraw/test_cli.py -v -k validate_room_arg`

Expected: import/signature mismatch.

- [ ] **Step 3: Rename + rewrite the helpers in `src/zndraw/cli.py`**

Replace `_resolve_owner_id` and `_validate_room_arg`:

```python
def _resolve_owner_name(server_url: str, token: str) -> str:
    """Fetch the authenticated user's display_name from /v1/auth/users/me."""
    with httpx.Client(base_url=server_url, timeout=30.0) as client:
        resp = client.get(
            "/v1/auth/users/me",
            headers={"Authorization": f"Bearer {token}"},
        )
        resp.raise_for_status()
        return resp.json()["display_name"]


_DISPLAY_NAME_RE = re.compile(r"^[a-z][a-z0-9-]{2,63}$")


def _validate_room_arg(value: str, owner_name: str) -> None:
    """Validate that --room is in '<display_name>/<name>' form."""
    if "/" not in value:
        raise typer.BadParameter(
            f"--room must be '<display-name>/<name>'. Got '{value}'.\n"
            f"Your display name is {owner_name}. Try: --room {owner_name}/{value}"
        )
    owner_part, _, name_part = value.partition("/")
    if not _DISPLAY_NAME_RE.fullmatch(owner_part):
        raise typer.BadParameter(
            f"--room owner '{owner_part}' is not a valid display name. "
            f"Your display name is {owner_name}."
        )
    if not re.fullmatch(r"[a-zA-Z0-9\-_]+", name_part):
        raise typer.BadParameter(
            f"--room name '{name_part}' contains invalid characters. "
            f"Allowed: letters, digits, '-', '_'."
        )
```

In the `main()` body, replace `owner_id = _resolve_owner_id(url, token)` with `owner_name = _resolve_owner_name(url, token)`, then update downstream:

```python
if room is not None:
    _validate_room_arg(room, owner_name)

if room is not None:
    room_names = [room] * len(path or [])
else:
    simple_names = get_room_names(path or [], room, append)
    room_names = [f"{owner_name}/{n}" for n in simple_names]

first_room = (
    room_names[0]
    if room_names
    else f"{owner_name}/workspace-{uuid.uuid4().hex[:8]}"
)
```

(`uuid` import stays — used for room-name suffix randomness.)

Remove the `from uuid import UUID` import and any leftover `UUID` references in this file.

- [ ] **Step 4: Confirm CLI tests pass**

Run: `uv run pytest tests/zndraw/test_cli.py -v`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/cli.py tests/zndraw/test_cli.py
git commit -m "feat(cli): use display_name in room addresses"
```

---

### Task 14: Frontend — UserInfo type, store, auth flow

**Files:**
- Modify: `frontend/src/utils/auth.ts`
- Modify: `frontend/src/store.tsx`

- [ ] **Step 1: Extend the UserInfo type**

In `frontend/src/utils/auth.ts`:

```typescript
export interface UserInfo {
    id: string;
    email: string;
    display_name: string;
    is_active: boolean;
    is_superuser: boolean;
    is_verified: boolean;
}
```

The `login()` function uses `fetchMeWithToken` which already returns whatever `/v1/auth/users/me` produces — backend Task 2 added `display_name`, so the field flows through unchanged.

For the guest path inside `login()`, prefer the new field:

```typescript
if (!password) {
    const response = await fetch("/v1/auth/guest", { method: "POST" });
    ...
    const data = await response.json();
    // data is GuestSessionResponse: {access_token, token_type, email, display_name}
    const user = await fetchMeWithToken(data.access_token);
    localStorage.setItem(TOKEN_KEY, data.access_token);
    return { token: data.access_token, user };
}
```

`fetchMeWithToken` returns the full UserInfo including display_name, so no further plumbing is needed.

- [ ] **Step 2: Rename `userEmail` selector to `userDisplayName`**

In `frontend/src/store.tsx`:

```typescript
export const selectIsRoomReadOnly = (state: AppState): boolean => {
    const isSuperuser = state.user?.is_superuser ?? false;
    if (isSuperuser) return false;
    if (state.userLock) {
        const userDisplayName = state.user?.display_name ?? null;
        return state.userLock !== userDisplayName;
    }
    return false;
};
```

- [ ] **Step 3: Type-check**

Run: `cd frontend && bun run typecheck`

Expected: errors will fan out across consumers (ChatPanel, UserProfileDialog, RegisterDialog) — fixed in the next tasks. No new errors in `store.tsx` itself.

- [ ] **Step 4: Commit**

```bash
git add frontend/src/utils/auth.ts frontend/src/store.tsx
git commit -m "feat(frontend): UserInfo.display_name + selector rename"
```

---

### Task 15: Frontend — Registration dialog with coolname suggestion

**Files:**
- Modify: `frontend/src/components/RegisterDialog.tsx`

- [ ] **Step 1: Fetch a suggestion when the dialog opens**

Replace the body of `RegisterDialog.tsx` so it:
- Adds a `displayName` state, populated when `open` flips true via `GET /v1/users/available-display-name`.
- Validates it client-side against the regex.
- Includes it in the `registerUser` call.
- Adds a "Regenerate" icon button that re-fetches a suggestion.

```typescript
import RefreshIcon from "@mui/icons-material/Refresh";
import {
    Alert,
    Box,
    Button,
    Dialog,
    DialogActions,
    DialogContent,
    DialogTitle,
    IconButton,
    InputAdornment,
    TextField,
    Tooltip,
    Typography,
} from "@mui/material";
import type React from "react";
import { useEffect, useState } from "react";
import { connectWithAuth } from "../socket";
import { useAppStore } from "../store";
import { registerUser } from "../utils/auth";

const DISPLAY_NAME_RE = /^[a-z][a-z0-9-]{2,63}$/;

async function fetchSuggestion(): Promise<string> {
    const resp = await fetch("/v1/users/available-display-name");
    if (!resp.ok) throw new Error("Failed to suggest a display name");
    const data = await resp.json();
    return data.display_name as string;
}

interface RegisterDialogProps {
    open: boolean;
    onClose: () => void;
}

export default function RegisterDialog({ open, onClose }: RegisterDialogProps) {
    const [email, setEmail] = useState("");
    const [displayName, setDisplayName] = useState("");
    const [password, setPassword] = useState("");
    const [passwordConfirm, setPasswordConfirm] = useState("");
    const [error, setError] = useState<string | null>(null);
    const [loading, setLoading] = useState(false);

    const setUser = useAppStore((state) => state.setUser);
    const showSnackbar = useAppStore((state) => state.showSnackbar);
    const userDisplayName = useAppStore((state) => state.user?.display_name ?? null);

    useEffect(() => {
        if (!open) return;
        fetchSuggestion().then(setDisplayName).catch(() => {});
    }, [open]);

    const regenerate = async () => {
        try {
            setDisplayName(await fetchSuggestion());
        } catch (err) {
            // Silent: keep the user's typed value.
        }
    };

    const handleRegister = async () => {
        setError(null);
        if (!email.trim()) return setError("Email is required");
        if (!DISPLAY_NAME_RE.test(displayName))
            return setError(
                "Display name must be 3–64 chars, lowercase, digits and hyphens, starting with a letter",
            );
        if (!password) return setError("Password is required");
        if (password !== passwordConfirm)
            return setError("Passwords do not match");

        setLoading(true);
        try {
            await registerUser(email, password, displayName);
            const { user } = await connectWithAuth();
            setUser(user);
            showSnackbar(`Registered as ${user.display_name}`, "success");
            onClose();
            setEmail("");
            setDisplayName("");
            setPassword("");
            setPasswordConfirm("");
        } catch (err) {
            const message = err instanceof Error ? err.message : "Registration failed";
            setError(message);
        } finally {
            setLoading(false);
        }
    };

    const handleClose = () => {
        if (!loading) {
            setError(null);
            setEmail("");
            setDisplayName("");
            setPassword("");
            setPasswordConfirm("");
            onClose();
        }
    };

    const handleKeyDown = (event: React.KeyboardEvent) => {
        if (event.key === "Enter" && !loading) {
            if (email && displayName && password && passwordConfirm) handleRegister();
        }
    };

    return (
        <Dialog open={open} onClose={handleClose} maxWidth="xs" fullWidth>
            <DialogTitle>Register Account</DialogTitle>
            <DialogContent>
                <Box sx={{ pt: 1, display: "flex", flexDirection: "column", gap: 2 }}>
                    <Typography variant="body2" color="text.secondary">
                        Current temporary name: <strong>{userDisplayName}</strong>
                    </Typography>
                    <Typography variant="body2" color="text.secondary">
                        Pick a display name (or accept the suggestion), then enter an
                        email and password.
                    </Typography>

                    {error && (
                        <Alert severity="error" onClose={() => setError(null)}>
                            {error}
                        </Alert>
                    )}

                    <TextField
                        label="Display name"
                        value={displayName}
                        onChange={(e) => setDisplayName(e.target.value)}
                        onKeyDown={handleKeyDown}
                        disabled={loading}
                        fullWidth
                        autoComplete="off"
                        InputProps={{
                            endAdornment: (
                                <InputAdornment position="end">
                                    <Tooltip title="Regenerate suggestion">
                                        <IconButton
                                            size="small"
                                            onClick={regenerate}
                                            disabled={loading}
                                        >
                                            <RefreshIcon fontSize="small" />
                                        </IconButton>
                                    </Tooltip>
                                </InputAdornment>
                            ),
                        }}
                    />

                    <TextField
                        label="Email"
                        type="email"
                        value={email}
                        onChange={(e) => setEmail(e.target.value)}
                        onKeyDown={handleKeyDown}
                        disabled={loading}
                        fullWidth
                        autoComplete="email"
                    />

                    <TextField
                        label="Password"
                        type="password"
                        value={password}
                        onChange={(e) => setPassword(e.target.value)}
                        onKeyDown={handleKeyDown}
                        disabled={loading}
                        fullWidth
                        autoComplete="new-password"
                    />

                    <TextField
                        label="Confirm Password"
                        type="password"
                        value={passwordConfirm}
                        onChange={(e) => setPasswordConfirm(e.target.value)}
                        onKeyDown={handleKeyDown}
                        disabled={loading}
                        fullWidth
                        autoComplete="new-password"
                    />
                </Box>
            </DialogContent>
            <DialogActions sx={{ px: 3, pb: 2 }}>
                <Button onClick={handleClose} disabled={loading}>
                    Cancel
                </Button>
                <Button
                    onClick={handleRegister}
                    disabled={
                        loading || !email || !displayName || !password || !passwordConfirm
                    }
                    variant="contained"
                >
                    {loading ? "Registering..." : "Register"}
                </Button>
            </DialogActions>
        </Dialog>
    );
}
```

- [ ] **Step 2: Extend `registerUser` to accept display_name**

In `frontend/src/utils/auth.ts`:

```typescript
export async function registerUser(
    email: string,
    password: string,
    display_name: string,
): Promise<AuthResult> {
    const response = await fetch("/v1/auth/register", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ email, password, display_name }),
    });
    if (!response.ok) {
        const errorData = await response
            .json()
            .catch(() => ({ detail: response.statusText }));
        // 409 → display-name taken; surface to caller as-is
        throw new Error(
            errorData.detail || `Registration failed: ${response.statusText}`,
        );
    }
    return login(email, password);
}
```

- [ ] **Step 3: Commit**

```bash
git add frontend/src/components/RegisterDialog.tsx frontend/src/utils/auth.ts
git commit -m "feat(frontend): RegisterDialog drives display_name suggestion + submit"
```

---

### Task 16: Frontend — UserProfileDialog + ChatPanel

**Files:**
- Modify: `frontend/src/components/UserProfileDialog.tsx`
- Modify: `frontend/src/panels/ChatPanel.tsx`
- Modify: `frontend/src/types/chat.ts`

- [ ] **Step 1: UserProfileDialog**

Replace the selector at `frontend/src/components/UserProfileDialog.tsx:35`:

```typescript
const userDisplayName = useAppStore((state) => state.user?.display_name ?? null);
const userEmail = useAppStore((state) => state.user?.email ?? null);
```

Render `userDisplayName` as the primary identity line and `userEmail` as the secondary line (where the existing layout currently uses `userName`). Mechanical update — preserve existing styling.

- [ ] **Step 2: ChatPanel selector + ownership check + typing emits**

In `frontend/src/panels/ChatPanel.tsx`:

```typescript
const userDisplayName = useAppStore((state) => state.user?.display_name ?? null);
...
const isOwnMessage = message.display_name === userDisplayName;
...
<Typography variant="caption" color="text.secondary">
    {message.display_name} – {format(new Date(message.created_at), "HH:mm")}
</Typography>
```

Replace `socket.emit("typing_start", { owner_id: ownerId, room_name: roomName })` with `socket.emit("typing_start", { owner: ownerId, room_name: roomName })` and the same for `typing_stop`.

(`ownerId` from `useParams` is the display-name segment after Task 7's URL flip — verify the path-routing config uses the same name; URLs may be `/room/:ownerId/:roomName`. If yes, no router config change needed — only the variable contents change.)

- [ ] **Step 3: Update ChatMessage type**

`frontend/src/types/chat.ts`:

```typescript
export interface ChatMessage {
    id: number;
    room_id: string;
    user_id: string;
    content: string;
    created_at: string;
    updated_at: string | null;
    display_name: string | null;
}

export interface MessageNewEvent {
    id: number;
    room_id: string;
    user_id: string;
    content: string;
    created_at: string;
    updated_at: string | null;
    display_name: string | null;
}
```

- [ ] **Step 4: Type-check + smoke run**

Run: `cd frontend && bun run typecheck`

Expected: errors remaining are limited to `myapi/client.ts` (Task 17) and the room-create callers (Task 18).

- [ ] **Step 5: Commit**

```bash
git add frontend/src/components/UserProfileDialog.tsx frontend/src/panels/ChatPanel.tsx frontend/src/types/chat.ts
git commit -m "feat(frontend): chat + profile display the display_name"
```

---

### Task 17: Frontend — client.ts type renames

**Files:**
- Modify: `frontend/src/myapi/client.ts`

- [ ] **Step 1: Rename fields**

Edit `frontend/src/myapi/client.ts` around the line ranges listed in the spec (425–823):

- `CreateRoomRequest`: `owner_id: string` → `owner: string` (and update the inline `room_id` comment "composed: `{owner_id}/{room_name}`" → `{owner}/{room_name}`).
- `RoomInfo`: `owner_id: string` → `owner: string`.
- `Room` (interface around line 805): `owner_id: string` → `owner: string`; comment update on `room_id`.
- `RoomUpdateRequest`: `new_owner_id?: string | null` → `new_owner?: string | null`.
- `GroupMember`: `email: string | null` → `display_name: string | null`.

- [ ] **Step 2: Type-check**

Run: `cd frontend && bun run typecheck`

Expected: errors now point at room-create callers (next task).

- [ ] **Step 3: Commit**

```bash
git add frontend/src/myapi/client.ts
git commit -m "refactor(frontend): rename owner_id → owner in API types"
```

---

### Task 18: Frontend — switch callers to `currentUser.display_name`

**Files:**
- Modify: `frontend/src/panels/RoomsPanel.tsx`
- Modify: `frontend/src/panels/roomsHeaderActions.tsx`
- Modify: `frontend/src/panels/FilesystemPanel.tsx`
- Modify: `frontend/src/components/DuplicateRoomDialog.tsx`
- Modify: `frontend/src/pages/templateSelection.tsx`
- Modify: `frontend/src/hooks/socketHandlers/connectionHandlers.ts`
- Modify: `frontend/src/roomsStore.tsx`

- [ ] **Step 1: RoomsPanel**

In `frontend/src/panels/RoomsPanel.tsx:57` (inside a `createRoom` call):

```typescript
createRoom({ owner: currentUser.display_name, name, ... });
```

- [ ] **Step 2: roomsHeaderActions**

In `frontend/src/panels/roomsHeaderActions.tsx`, lines 25/43/65 — same pattern.

- [ ] **Step 3: FilesystemPanel**

In `frontend/src/panels/FilesystemPanel.tsx:138` — same pattern.

- [ ] **Step 4: DuplicateRoomDialog**

In `frontend/src/components/DuplicateRoomDialog.tsx:60` — same pattern.

- [ ] **Step 5: templateSelection**

In `frontend/src/pages/templateSelection.tsx:117,125`:

```typescript
const owner = user.display_name;
...
console.log("[Startup] Creating default room for user:", owner);
await createRoom({ owner, name, ... });
```

- [ ] **Step 6: connectionHandlers**

In `frontend/src/hooks/socketHandlers/connectionHandlers.ts:201,219,239`:

```typescript
socket.emit("room_join", { owner: ownerId, room_name: roomName, client_type: "frontend" });
...
socket.emit("typing_start", { owner: ownerId, room_name: roomName });
...
socket.emit("typing_stop", { owner: ownerId, room_name: roomName });
```

(`ownerId` here comes from `useParams` and is the URL display-name segment.)

- [ ] **Step 7: roomsStore**

In `frontend/src/roomsStore.tsx:68,80`:

```typescript
if (
    updates.owner === undefined ||
    ...
) {
    ...
}
...
const payload = {
    ...
    owner: updates.owner,
};
```

- [ ] **Step 8: Type-check + build**

Run: `cd frontend && bun run typecheck && bun run build`

Expected: clean.

- [ ] **Step 9: Commit**

```bash
git add frontend/src/panels/ frontend/src/components/DuplicateRoomDialog.tsx frontend/src/pages/templateSelection.tsx frontend/src/hooks/socketHandlers/connectionHandlers.ts frontend/src/roomsStore.tsx
git commit -m "feat(frontend): callers send display_name as owner"
```

---

### Task 19: Frontend — E2E tests for registration + room URLs

**Files:**
- Add or modify: `frontend/e2e/` test files (exact filenames depend on existing layout)

- [ ] **Step 1: Identify existing E2E suite**

Run: `ls frontend/e2e/`

Expected: existing Playwright-style tests under `frontend/e2e/`. (If none, fall back to component-level tests with React Testing Library where appropriate.)

- [ ] **Step 2: Registration coverage**

Add a Playwright test asserting:
- Opening the register dialog issues `GET /v1/users/available-display-name` and pre-fills the input.
- Clicking the regenerate icon issues another `GET` and updates the input.
- Submitting with a valid email + password + display_name lands on a fresh room URL using the display_name as the first segment.

- [ ] **Step 3: Room URL coverage**

In an existing room-URL fixture, swap UUID segments for display-name segments. Confirm the page renders normally.

- [ ] **Step 4: Run E2E (if a runner is configured)**

Run: `cd frontend && bun run test:e2e` (or the project-configured equivalent).

Expected: green.

- [ ] **Step 5: Commit**

```bash
git add frontend/e2e/
git commit -m "test(frontend): E2E coverage for display-name registration + room URLs"
```

---

### Task 20: Final sweep — full test suite + manual smoke

**Files:** none (verification only).

- [ ] **Step 1: Backend full suite**

Run: `uv run pytest -x -q`

Expected: green. **Stop and surface to the user if a protected test fails.**

- [ ] **Step 2: Frontend type-check + build**

Run: `cd frontend && bun run typecheck && bun run build`

Expected: clean.

- [ ] **Step 3: Lint / format pre-commit**

Run: `uvx prek --all-files`

Expected: green (or only auto-fixable noise — stage the fixes, re-run).

- [ ] **Step 4: Manual smoke** (per project rule: "For UI or frontend changes, start the dev server and use the feature in a browser before reporting the task as complete")

Start the stack (`uv run zndraw` + the configured frontend dev server) and verify:
1. Visiting the root URL as a fresh guest produces a friendly URL like `/rooms/<display-name>/<workspace>`.
2. The register dialog suggests a display name; the regenerate icon swaps it; submitting persists.
3. After logging in, the URL for the user's default room uses the display name.
4. Chat: own messages show your display name; another browser tab as a second user shows their display name.
5. Console: no 422/404 surprises, no `owner_id is not defined` warnings.

Record what was checked in the PR description; mention any items that could not be verified (e.g., second-user flow without a partner browser).

- [ ] **Step 5: Final commit if there are formatting/lint changes**

```bash
git add -A
git commit -m "chore: post-rollout sweep — lint/format fixes"
```

---

## After execution

When all tasks pass:
- Run `git log --oneline main..HEAD` to review the commit graph.
- Open a PR via `gh pr create` referencing issue #931. Reference the spec (`docs/superpowers/specs/2026-05-27-coolname-display-names-design.html`). Include the manual-smoke checklist.

## Out of scope (explicit)

- PATCH route for user-initiated display-name rename + UI.
- Backfill migration script for existing DBs.
- Backwards-compat acceptance of UUID-shaped owner segments.
- Visual redesign of the profile/rooms screens beyond replacing rendered identifiers.
