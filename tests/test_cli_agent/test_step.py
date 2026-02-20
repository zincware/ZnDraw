"""Tests for CLI step commands."""

from __future__ import annotations

from typer.testing import CliRunner

from zndraw.schemas import StepResponse, StepUpdateResponse

from .conftest import invoke_cli


def test_step_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """step get should return a valid StepResponse."""
    data = invoke_cli(cli_runner, server_url, auth_token, ["step", "get", test_room])
    resp = StepResponse.model_validate(data)
    assert resp.step >= 0
    assert resp.total_frames >= 1


def test_step_set(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """step set should update the step and return StepUpdateResponse."""
    data = invoke_cli(
        cli_runner, server_url, auth_token, ["step", "set", test_room, "0"]
    )
    resp = StepUpdateResponse.model_validate(data)
    assert resp.success is True
    assert resp.step == 0


def test_step_set_then_get(
    cli_runner: CliRunner, server_url: str, auth_token: str, test_room: str
) -> None:
    """Setting step then getting it should reflect the change."""
    invoke_cli(cli_runner, server_url, auth_token, ["step", "set", test_room, "0"])
    data = invoke_cli(cli_runner, server_url, auth_token, ["step", "get", test_room])
    resp = StepResponse.model_validate(data)
    assert resp.step == 0
