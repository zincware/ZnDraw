"""End-to-end UserManager.create coverage for display_name handling."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
from sqlmodel import select

from zndraw_auth.db import User
from zndraw_auth.display_names import DISPLAY_NAME_PATTERN

if TYPE_CHECKING:
    from httpx import AsyncClient
    from sqlmodel.ext.asyncio.session import AsyncSession


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
