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
