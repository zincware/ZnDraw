"""Tests for CLI room naming behavior with --append and --room flags."""

import re

import pytest
from typer.testing import CliRunner

from zndraw import ZnDraw
from zndraw.cli import app
from zndraw.utils import path_to_room, sanitize_room_name

runner = CliRunner()


# =============================================================================
# Unit tests for room naming functions
# =============================================================================


@pytest.mark.parametrize(
    "input_name,expected",
    [
        ("tmp/s22.xyz", "tmp_s22_xyz"),
        ("data/file.pdb", "data_file_pdb"),
        ("./structure.xyz", "__structure_xyz"),
        ("file with spaces.xyz", "file_with_spaces_xyz"),
        ("file@special#chars.xyz", "file_special_chars_xyz"),
        ("simple", "simple"),
        ("my-file-name.xyz", "my-file-name_xyz"),
    ],
)
def test_sanitize_room_name(input_name, expected):
    """Test that sanitize_room_name correctly replaces non-alphanumeric characters."""
    assert sanitize_room_name(input_name) == expected


def test_path_to_room_unique_generates_different_names():
    """Test that path_to_room with unique=True generates different names each time."""
    room1 = path_to_room("tmp/s22.xyz", unique=True)
    room2 = path_to_room("tmp/s22.xyz", unique=True)

    assert room1 != room2
    assert room1.startswith("tmp_s22_xyz_")
    assert room2.startswith("tmp_s22_xyz_")


def test_path_to_room_unique_format():
    """Test that unique room names have correct format with 4-char hex suffix."""
    room = path_to_room("data/file.xyz", unique=True)

    # Should match pattern: sanitized_name_XXXX where XXXX is 4 hex chars
    pattern = r"^data_file_xyz_[a-f0-9]{4}$"
    assert re.match(pattern, room), f"Room '{room}' doesn't match expected pattern"


def test_path_to_room_not_unique_deterministic():
    """Test that path_to_room with unique=False returns same name each time."""
    room1 = path_to_room("tmp/s22.xyz", unique=False)
    room2 = path_to_room("tmp/s22.xyz", unique=False)

    assert room1 == room2
    assert room1 == "tmp_s22_xyz"


# =============================================================================
# CLI validation tests (no server needed)
# =============================================================================


def test_cli_append_and_room_mutually_exclusive(tmp_path):
    """Test that --append and --room cannot be used together."""
    test_file = tmp_path / "test.xyz"
    test_file.write_text("dummy")

    result = runner.invoke(app, [str(test_file), "--append", "--room", "myroom"])

    assert result.exit_code == 1
    assert "cannot be used together" in result.output


def test_cli_append_requires_file():
    """Test that --append requires file path(s) to be specified."""
    result = runner.invoke(app, ["--append"])

    assert result.exit_code == 1
    assert "require file path" in result.output


def test_cli_room_requires_file():
    """Test that --room requires file path(s) to be specified."""
    result = runner.invoke(app, ["--room", "myroom"])

    assert result.exit_code == 1
    assert "require file path" in result.output


# =============================================================================
# Integration tests (server needed)
# =============================================================================


def test_cli_default_creates_unique_room(server, s22_xyz, s22):
    """Test that default CLI behavior creates a unique room each time."""
    # Run CLI twice with same file
    result1 = runner.invoke(
        app, [s22_xyz, "--connect", server, "--no-browser"], catch_exceptions=False
    )
    assert result1.exit_code == 0

    result2 = runner.invoke(
        app, [s22_xyz, "--connect", server, "--no-browser"], catch_exceptions=False
    )
    assert result2.exit_code == 0

    # Extract room names from output
    room1 = _extract_room_from_output(result1.output)
    room2 = _extract_room_from_output(result2.output)

    assert room1 != room2, "Default behavior should create unique rooms"
    assert room1.startswith("tmp_") or "s22" in room1
    assert room2.startswith("tmp_") or "s22" in room2

    assert len(ZnDraw(room=room1, url=server, user="tester")) == len(s22)
    assert len(ZnDraw(room=room2, url=server, user="tester")) == len(s22)


def test_cli_append_creates_deterministic_room(server, s22_xyz, s22):
    """Test that --append creates the same room name each time."""
    # Run CLI twice with --append
    result1 = runner.invoke(
        app,
        [s22_xyz, "--append", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result1.exit_code == 0

    result2 = runner.invoke(
        app,
        [s22_xyz, "--append", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result2.exit_code == 0

    # Extract room names from output
    room1 = _extract_room_from_output(result1.output)
    room2 = _extract_room_from_output(result2.output)

    assert room1 == room2, "--append should create same room name"


def test_cli_append_appends_data_to_existing_room(server, s22_xyz, s22):
    """Test that --append actually appends data to an existing room."""
    # First upload
    result1 = runner.invoke(
        app,
        [s22_xyz, "--append", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result1.exit_code == 0
    room = _extract_room_from_output(result1.output)

    # Check initial frame count
    vis1 = ZnDraw(room=room, url=server, user="tester")
    initial_count = len(vis1)
    assert initial_count == len(s22)

    # Second upload (should append)
    result2 = runner.invoke(
        app,
        [s22_xyz, "--append", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result2.exit_code == 0

    # Check frame count after append
    vis2 = ZnDraw(room=room, url=server, user="tester")
    final_count = len(vis2)
    assert final_count == 2 * len(s22), "Data should be appended"


def test_cli_room_uses_explicit_name(server, s22_xyz, s22):
    """Test that --room uses the exact room name specified."""
    room_name = "my-custom-room"

    result = runner.invoke(
        app,
        [s22_xyz, "--room", room_name, "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0
    assert room_name in result.output

    # Verify data is in the specified room
    vis = ZnDraw(room=room_name, url=server, user="tester")
    assert len(vis) == len(s22)


def test_cli_room_sanitizes_name(server, s22_xyz, s22):
    """Test that --room sanitizes the room name (replaces special chars)."""
    # Use a name with special characters
    result = runner.invoke(
        app,
        [s22_xyz, "--room", "my/room@name", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0

    # Room name should be sanitized
    sanitized_name = "my_room_name"
    vis = ZnDraw(room=sanitized_name, url=server, user="tester")
    assert len(vis) == len(s22)


def test_cli_room_multiple_files_same_room(server, s22_xyz, s22, tmp_path):
    """Test that --room puts multiple files into the same room."""
    import ase.io

    # Create a second file
    second_file = tmp_path / "second.xyz"
    ase.io.write(second_file, s22[:5])  # Write first 5 structures

    room_name = "combined-room"
    result = runner.invoke(
        app,
        [
            s22_xyz,
            str(second_file),
            "--room",
            room_name,
            "--connect",
            server,
            "--no-browser",
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0

    # Both files should be in the same room
    vis = ZnDraw(room=room_name, url=server, user="tester")
    expected_count = len(s22) + 5  # Full s22 + 5 from second file
    assert len(vis) == expected_count


def test_cli_append_multiple_files_separate_rooms(server, s22_xyz, s22, tmp_path):
    """Test that --append with multiple files creates separate deterministic rooms."""
    import ase.io

    # Create a second file
    second_file = tmp_path / "second.xyz"
    ase.io.write(second_file, s22[:5])

    result = runner.invoke(
        app,
        [s22_xyz, str(second_file), "--append", "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0

    # Each file should have its own room
    room1 = path_to_room(s22_xyz, unique=False)
    room2 = path_to_room(str(second_file), unique=False)

    vis1 = ZnDraw(room=room1, url=server, user="tester")
    vis2 = ZnDraw(room=room2, url=server, user="tester")

    assert len(vis1) == len(s22)
    assert len(vis2) == 5


def test_cli_default_multiple_files_separate_unique_rooms(
    server, s22_xyz, s22, tmp_path
):
    """Test that default behavior with multiple files creates separate unique rooms."""
    import ase.io

    # Create a second file
    second_file = tmp_path / "second.xyz"
    ase.io.write(second_file, s22[:5])

    result = runner.invoke(
        app,
        [s22_xyz, str(second_file), "--connect", server, "--no-browser"],
        catch_exceptions=False,
    )
    assert result.exit_code == 0

    # Extract both room names from output
    rooms = _extract_all_rooms_from_output(result.output)
    assert len(rooms) == 2
    assert rooms[0] != rooms[1], "Each file should get a unique room"


# =============================================================================
# Helper functions
# =============================================================================


def _extract_room_from_output(output: str) -> str:
    """Extract the first room name from CLI output."""
    # Look for "Uploading file ... to room <room_name>"
    match = re.search(r"to room (\S+)", output)
    if match:
        return match.group(1)
    raise ValueError(f"Could not extract room name from output: {output}")


def _extract_all_rooms_from_output(output: str) -> list[str]:
    """Extract all room names from CLI output."""
    matches = re.findall(r"to room (\S+)", output)
    if matches:
        return matches
    raise ValueError(f"Could not extract room names from output: {output}")
