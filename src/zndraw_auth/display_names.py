"""Display-name generation and validation."""

from __future__ import annotations

import re
import secrets
from typing import TYPE_CHECKING, Final

import coolname
from sqlmodel import select

from zndraw_auth.db import User

if TYPE_CHECKING:
    from sqlmodel.ext.asyncio.session import AsyncSession

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
    # Fallback: append a random hex suffix and validate the same way the
    # looped slugs were — reserved-set + regex + DB uniqueness.
    for _ in range(max_attempts):
        candidate = f"{coolname.generate_slug(3)}-{secrets.token_hex(2)}"
        if candidate in RESERVED_DISPLAY_NAMES:
            continue
        if not DISPLAY_NAME_PATTERN.fullmatch(candidate):
            continue
        exists = await session.scalar(
            select(User.id).where(User.display_name == candidate).limit(1)
        )
        if exists is None:
            return candidate
    raise RuntimeError(
        "Could not generate a unique display name after fallback attempts"
    )
