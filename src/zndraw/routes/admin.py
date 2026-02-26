"""Admin REST API endpoints for user management and server control.

All endpoints require admin privileges via AdminUserDep.
"""

import asyncio
import os
import signal
from typing import Annotated
from uuid import UUID

from fastapi import APIRouter, Query
from pydantic import BaseModel
from sqlmodel import func, select
from zndraw_auth import User

from zndraw.dependencies import AdminUserDep, SessionDep
from zndraw.exceptions import Forbidden, UserNotFound, problem_responses
from zndraw.schemas import OffsetPage, StatusResponse

router = APIRouter(prefix="/v1/admin", tags=["admin"])


# =============================================================================
# Response Schemas
# =============================================================================


class AdminUserResponse(BaseModel):
    """Admin view of a user with full details."""

    id: UUID
    email: str
    is_superuser: bool

    model_config = {"from_attributes": True}


class UserUpdateRequest(BaseModel):
    """Request to update a user."""

    is_superuser: bool | None = None


class ShutdownResponse(BaseModel):
    """Response for shutdown request."""

    status: str = "shutting_down"
    message: str = "Server shutdown initiated"


# =============================================================================
# User Management Endpoints
# =============================================================================


@router.get(
    "/users",
    responses=problem_responses(Forbidden),
)
async def list_users(
    session: SessionDep,
    _admin: AdminUserDep,
    offset: Annotated[int, Query(ge=0)] = 0,
    limit: Annotated[int, Query(ge=1, le=100)] = 50,
) -> OffsetPage[AdminUserResponse]:
    """List all users with pagination.

    Requires admin privileges.
    """
    # Get total count
    count_result = await session.execute(select(func.count()).select_from(User))
    total = count_result.scalar_one()

    # Get paginated users
    result = await session.execute(select(User).offset(offset).limit(limit))
    users = list(result.scalars().all())

    return OffsetPage(
        items=[AdminUserResponse.model_validate(u) for u in users],
        total=total,
        offset=offset,
        limit=limit,
    )


@router.get(
    "/users/{user_id}",
    response_model=AdminUserResponse,
    responses=problem_responses(Forbidden, UserNotFound),
)
async def get_user(
    session: SessionDep,
    _admin: AdminUserDep,
    user_id: UUID,
) -> AdminUserResponse:
    """Get a specific user by ID.

    Requires admin privileges.
    """
    user = await session.get(User, user_id)
    if user is None:
        raise UserNotFound.exception(f"User with id {user_id} not found")
    return AdminUserResponse.model_validate(user)


@router.patch(
    "/users/{user_id}",
    response_model=AdminUserResponse,
    responses=problem_responses(Forbidden, UserNotFound),
)
async def update_user(
    session: SessionDep,
    admin: AdminUserDep,
    user_id: UUID,
    request: UserUpdateRequest,
) -> AdminUserResponse:
    """Update a user (promote/demote superuser).

    Requires admin privileges.

    Note: Admins cannot demote themselves to prevent lockout.
    """
    user = await session.get(User, user_id)
    if user is None:
        raise UserNotFound.exception(f"User with id {user_id} not found")

    # Prevent self-demotion
    if (
        request.is_superuser is not None
        and user.id == admin.id
        and not request.is_superuser
    ):
        raise Forbidden.exception("Cannot demote yourself")

    # Apply updates
    if request.is_superuser is not None:
        user.is_superuser = request.is_superuser

    await session.commit()
    await session.refresh(user)

    return AdminUserResponse.model_validate(user)


@router.delete(
    "/users/{user_id}",
    responses=problem_responses(Forbidden, UserNotFound),
)
async def delete_user(
    session: SessionDep,
    admin: AdminUserDep,
    user_id: UUID,
) -> StatusResponse:
    """Delete a user.

    Requires admin privileges.

    Note: Admins cannot delete themselves.
    """
    user = await session.get(User, user_id)
    if user is None:
        raise UserNotFound.exception(f"User with id {user_id} not found")

    # Prevent self-deletion
    if user.id == admin.id:
        raise Forbidden.exception("Cannot delete yourself")

    await session.delete(user)
    await session.commit()

    return StatusResponse()


# =============================================================================
# Server Control Endpoints
# =============================================================================


@router.post(
    "/shutdown",
    response_model=ShutdownResponse,
    responses=problem_responses(Forbidden),
)
async def shutdown_server(
    _admin: AdminUserDep,
) -> ShutdownResponse:
    """Initiate graceful server shutdown.

    Requires admin privileges.

    The server will complete current requests before shutting down.
    """

    # Schedule shutdown after response is sent
    async def delayed_shutdown() -> None:
        await asyncio.sleep(0.5)  # Small delay to allow response to be sent
        os.kill(os.getpid(), signal.SIGTERM)

    asyncio.create_task(delayed_shutdown())

    return ShutdownResponse()
