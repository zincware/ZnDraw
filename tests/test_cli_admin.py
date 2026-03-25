"""Tests for CLI admin commands (users list, users login)."""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from zndraw.cli_agent import app
from zndraw.cli_agent.connection import Connection
from zndraw.state_file import StateFile

runner = CliRunner()


@pytest.fixture
def state_file(tmp_path):
    return StateFile(directory=tmp_path)


def _make_connection(base_url="http://localhost:8000", token="admin.jwt.token"):  # noqa: S107
    """Create a mock Connection."""
    conn = MagicMock(spec=Connection)
    conn.base_url = base_url
    conn.token = token
    conn.client = MagicMock()
    return conn


# -- admin users list ----------------------------------------------------------


def test_admin_users_list():
    """admin users list should return user list."""
    list_resp = MagicMock()
    list_resp.status_code = 200
    list_resp.json.return_value = {
        "items": [
            {"id": "u1", "email": "user1@example.com", "is_superuser": False},
            {"id": "u2", "email": "admin@example.com", "is_superuser": True},
        ],
        "total": 2,
        "limit": 100,
        "offset": 0,
    }

    conn = _make_connection()
    conn.get.return_value = list_resp

    with patch(
        "zndraw.cli_agent.admin.get_connection",
        return_value=conn,
    ):
        result = runner.invoke(app, ["admin", "users", "list"])

    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout)
    assert "items" in data
    assert len(data["items"]) == 2


# -- admin users login ---------------------------------------------------------


def test_admin_users_login(state_file):
    """admin users login should mint token and save to state file."""
    mint_resp = MagicMock()
    mint_resp.status_code = 200
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

    conn = _make_connection()
    conn.post.return_value = mint_resp
    conn.client.get.return_value = me_resp

    with (
        patch("zndraw.cli_agent.admin.get_connection", return_value=conn),
        patch("zndraw.cli_agent.admin.StateFile", return_value=state_file),
    ):
        result = runner.invoke(app, ["admin", "users", "login", "target-id"])

    assert result.exit_code == 0, result.output
    stored = state_file.get_token("http://localhost:8000")
    assert stored is not None
    assert stored.access_token == "impersonated.jwt.token"
    assert stored.email == "target@example.com"
