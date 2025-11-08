"""Tests for file browser functionality."""

import tempfile
from pathlib import Path

import pytest

from zndraw.app.file_browser import is_supported_file, validate_path
from zndraw.auth import create_jwt_token
from zndraw.server import create_app


def get_auth_headers(app):
    """Get authentication headers for test client requests."""
    with app.app_context():
        token = create_jwt_token("test-user", role="user")
        return {"Authorization": f"Bearer {token}"}


@pytest.mark.parametrize(
    "requested_path,expected_valid",
    [
        ("", True),  # Root directory
        ("subdir", True),  # Valid subdirectory
        ("subdir/file.txt", True),  # Valid file in subdirectory
        ("../../../etc/passwd", False),  # Path traversal attempt
        ("/etc/passwd", False),  # Absolute path
        ("subdir/../../../etc/passwd", False),  # Path traversal with ..
    ],
)
def test_validate_path(requested_path, expected_valid, tmp_path):
    """Test path validation against traversal attacks."""
    root = tmp_path
    # Create a subdirectory
    (root / "subdir").mkdir(exist_ok=True)

    result = validate_path(requested_path, str(root))

    if expected_valid:
        assert result is not None
        assert result.is_relative_to(root)
    else:
        assert result is None


@pytest.mark.parametrize(
    "filename,expected",
    [
        ("test.xyz", True),
        ("test.pdb", True),
        ("test.h5md", True),
        ("test.h5", True),
        ("test.traj", True),
        ("test.extxyz", True),
        ("test.cif", True),
        ("test.txt", False),
        ("test.py", False),
        ("test.jpg", False),
        ("TEST.XYZ", True),  # Case insensitive
    ],
)
def test_is_supported_file(filename, expected):
    """Test file type detection."""
    filepath = Path(filename)
    assert is_supported_file(filepath) == expected


def test_list_endpoint_disabled():
    """Test that list endpoint returns 403 when feature is disabled."""
    app = create_app()
    # Modify the config object stored in extensions
    app.extensions["config"].file_browser_enabled = False
    app.extensions["config"].file_browser_root = "."

    with app.test_client() as client:
        # Without auth token, should get 401
        response = client.get("/api/file-browser/list")
        assert response.status_code == 401

        # With auth token, should get 403 (feature disabled)
        response = client.get("/api/file-browser/list", headers=get_auth_headers(app))
        assert response.status_code == 403
        assert "not enabled" in response.json["error"]


def test_load_endpoint_disabled():
    """Test that load endpoint returns 403 when feature is disabled."""
    app = create_app()
    # Modify the config object stored in extensions
    app.extensions["config"].file_browser_enabled = False
    app.extensions["config"].file_browser_root = "."

    with app.test_client() as client:
        # Without auth token, should get 401
        response = client.post("/api/file-browser/load", json={"path": "test.xyz"})
        assert response.status_code == 401

        # With auth token, should get 403 (feature disabled)
        response = client.post(
            "/api/file-browser/load",
            json={"path": "test.xyz"},
            headers=get_auth_headers(app),
        )
        assert response.status_code == 403
        assert "not enabled" in response.json["error"]


def test_list_endpoint_enabled(tmp_path):
    """Test list endpoint when feature is enabled."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    # Create test structure
    (tmp_path / "test.xyz").write_text("test")
    (tmp_path / "test.txt").write_text("test")
    (tmp_path / "subdir").mkdir()

    with app.test_client() as client:
        response = client.get("/api/file-browser/list", headers=get_auth_headers(app))
        assert response.status_code == 200
        data = response.json
        assert "items" in data
        assert "current_path" in data
        assert len(data["items"]) == 3  # 2 files + 1 directory


def test_list_endpoint_path_traversal(tmp_path):
    """Test that path traversal is blocked in list endpoint."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.get(
            "/api/file-browser/list?path=../../../etc", headers=get_auth_headers(app)
        )
        assert response.status_code == 400
        assert "Invalid path" in response.json["error"]


def test_list_endpoint_nonexistent_path(tmp_path):
    """Test list endpoint with non-existent path."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.get(
            "/api/file-browser/list?path=nonexistent", headers=get_auth_headers(app)
        )
        assert response.status_code == 404
        assert "does not exist" in response.json["error"]


def test_supported_types_endpoint(tmp_path):
    """Test supported-types endpoint."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.get("/api/file-browser/supported-types")
        assert response.status_code == 200
        data = response.json
        assert "extensions" in data
        assert "descriptions" in data
        assert ".xyz" in data["extensions"]
        assert ".pdb" in data["extensions"]


def test_load_endpoint_missing_path(tmp_path):
    """Test load endpoint with missing path field."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.post(
            "/api/file-browser/load", json={}, headers=get_auth_headers(app)
        )
        assert response.status_code == 400
        assert "Missing required field" in response.json["error"]


def test_load_endpoint_path_traversal(tmp_path):
    """Test that path traversal is blocked in load endpoint."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.post(
            "/api/file-browser/load",
            json={"path": "../../../etc/passwd"},
            headers=get_auth_headers(app),
        )
        assert response.status_code == 400
        assert "Invalid path" in response.json["error"]


def test_load_endpoint_nonexistent_file(tmp_path):
    """Test load endpoint with non-existent file."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    with app.test_client() as client:
        response = client.post(
            "/api/file-browser/load",
            json={"path": "nonexistent.xyz"},
            headers=get_auth_headers(app),
        )
        assert response.status_code == 404
        assert "does not exist" in response.json["error"]


def test_load_endpoint_unsupported_file(tmp_path):
    """Test load endpoint with unknown file type.

    Note: The endpoint now accepts all files and lets the backend reader handle
    errors. Unknown formats are attempted via ASE.
    """
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    # Create unknown file format
    (tmp_path / "test.txt").write_text("test")

    with app.test_client() as client:
        response = client.post(
            "/api/file-browser/load",
            json={"path": "test.txt"},
            headers=get_auth_headers(app),
        )
        # File is accepted for processing (backend will handle read errors)
        assert response.status_code == 200
        assert "task_id" in response.json


def test_list_endpoint_hidden_files(tmp_path):
    """Test that hidden files are filtered out."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    # Create visible and hidden files
    (tmp_path / "visible.xyz").write_text("test")
    (tmp_path / ".hidden.xyz").write_text("test")

    with app.test_client() as client:
        response = client.get("/api/file-browser/list", headers=get_auth_headers(app))
        assert response.status_code == 200
        data = response.json
        names = [item["name"] for item in data["items"]]
        assert "visible.xyz" in names
        assert ".hidden.xyz" not in names


def test_list_endpoint_sorted(tmp_path):
    """Test that list endpoint returns items sorted (directories first)."""
    app = create_app()
    app.extensions["config"].file_browser_enabled = True
    app.extensions["config"].file_browser_root = str(tmp_path)

    # Create files and directories
    (tmp_path / "z_file.xyz").write_text("test")
    (tmp_path / "a_dir").mkdir()
    (tmp_path / "b_file.xyz").write_text("test")

    with app.test_client() as client:
        response = client.get("/api/file-browser/list", headers=get_auth_headers(app))
        assert response.status_code == 200
        data = response.json
        names = [item["name"] for item in data["items"]]
        # Directories should come first
        assert names[0] == "a_dir"
        assert names[1] == "b_file.xyz"
        assert names[2] == "z_file.xyz"
