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
