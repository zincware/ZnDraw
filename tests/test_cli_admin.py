"""Tests for CLI admin commands (users list, users login)."""

from __future__ import annotations

import json
from datetime import UTC, datetime
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from zndraw.cli_agent import app
from zndraw.server_manager import TokenEntry, TokenStore

runner = CliRunner()


@pytest.fixture
def token_store(tmp_path):
    return TokenStore(directory=tmp_path)


@pytest.fixture
def admin_token_store(tmp_path):
    """TokenStore with an admin token pre-stored."""
    store = TokenStore(directory=tmp_path)
    store.set(
        "http://localhost:8000",
        TokenEntry(
            access_token="admin.jwt.token",
            email="admin@example.com",
            stored_at=datetime(2026, 3, 1, tzinfo=UTC),
        ),
    )
    return store


def _mock_client(**responses):
    """Create a mock httpx.Client context manager with configurable responses."""
    mock = MagicMock()
    mock.__enter__ = MagicMock(return_value=mock)
    mock.__exit__ = MagicMock(return_value=False)
    for method, response in responses.items():
        getattr(mock, method).return_value = response
    return mock


# ── admin users list ─────────────────────────────────────────────────


def test_admin_users_list(admin_token_store):
    """admin users list should return user list."""
    validate_resp = MagicMock(status_code=200)

    list_resp = MagicMock()
    list_resp.status_code = 200
    list_resp.raise_for_status = MagicMock()
    list_resp.json.return_value = {
        "items": [
            {"id": "u1", "email": "user1@example.com", "is_superuser": False},
            {"id": "u2", "email": "admin@example.com", "is_superuser": True},
        ],
        "total": 2,
        "limit": 100,
        "offset": 0,
    }

    mock = _mock_client(get=validate_resp)
    admin_mock = _mock_client(get=list_resp)

    with (
        patch(
            "zndraw.cli_agent.connection.get_token_store",
            return_value=admin_token_store,
        ),
        patch(
            "zndraw.cli_agent.admin.resolve_url", return_value="http://localhost:8000"
        ),
        patch("zndraw.cli_agent.connection.httpx.Client", return_value=mock),
        patch("zndraw.cli_agent.admin.httpx.Client", return_value=admin_mock),
    ):
        result = runner.invoke(app, ["admin", "users", "list"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert "items" in data
    assert len(data["items"]) == 2


# ── admin users login ────────────────────────────────────────────────


def test_admin_users_login(admin_token_store):
    """admin users login should mint token and save to store."""
    validate_resp = MagicMock(status_code=200)

    mint_resp = MagicMock()
    mint_resp.status_code = 200
    mint_resp.raise_for_status = MagicMock()
    mint_resp.json.return_value = {
        "access_token": "impersonated.jwt.token",
        "token_type": "bearer",
    }

    me_resp = MagicMock()
    me_resp.status_code = 200
    me_resp.json.return_value = {
        "id": "target-id",
        "email": "target@example.com",
        "is_superuser": False,
    }

    mock = _mock_client(get=validate_resp)

    admin_mock = _mock_client()
    admin_mock.post.return_value = mint_resp
    admin_mock.get.return_value = me_resp

    with (
        patch(
            "zndraw.cli_agent.connection.get_token_store",
            return_value=admin_token_store,
        ),
        patch(
            "zndraw.cli_agent.admin.get_token_store",
            return_value=admin_token_store,
        ),
        patch(
            "zndraw.cli_agent.admin.resolve_url", return_value="http://localhost:8000"
        ),
        patch("zndraw.cli_agent.connection.httpx.Client", return_value=mock),
        patch("zndraw.cli_agent.admin.httpx.Client", return_value=admin_mock),
    ):
        result = runner.invoke(app, ["admin", "users", "login", "target-id"])

    assert result.exit_code == 0, result.output
    stored = admin_token_store.get("http://localhost:8000")
    assert stored is not None
    assert stored.access_token == "impersonated.jwt.token"
    assert stored.email == "target@example.com"
