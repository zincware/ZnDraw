"""Tests for version compatibility checking between client and server."""

from unittest.mock import Mock, patch

import pytest

from zndraw.version_utils import (
    check_version_compatibility,
    parse_version,
    validate_server_version,
)


@pytest.mark.parametrize(
    "version,expected",
    [
        ("0.6.0", (0, 6, 0)),
        ("0.6.0a4", (0, 6, 0)),
        ("1.2.3", (1, 2, 3)),
        ("2.0.0b1", (2, 0, 0)),
        ("0.10.15rc2", (0, 10, 15)),
        ("invalid", None),
        ("1.2", None),
        ("1.a.3", None),
    ],
)
def test_parse_version(version, expected):
    """Test version string parsing."""
    assert parse_version(version) == expected


@pytest.mark.parametrize(
    "client,server,expected_compatible,expected_severity",
    [
        # Exact match
        ("0.6.0", "0.6.0", True, "info"),
        # Patch mismatch - warning
        ("0.6.0", "0.6.1", True, "warning"),
        ("0.6.1", "0.6.0", True, "warning"),
        # Minor mismatch - error
        ("0.6.0", "0.7.0", False, "error"),
        ("0.7.0", "0.6.0", False, "error"),
        # Major mismatch - error
        ("1.0.0", "2.0.0", False, "error"),
        ("2.0.0", "1.0.0", False, "error"),
        # Pre-release versions (pre-release suffixes are ignored in comparison)
        ("0.6.0a4", "0.6.0a5", True, "info"),  # Both parse to 0.6.0
        ("0.6.0a4", "0.6.0", True, "info"),  # Both parse to 0.6.0
        ("0.6.0a4", "0.6.0a4", True, "info"),
        ("0.6.0a4", "0.7.0a1", False, "error"),
        ("0.6.1a1", "0.6.0", True, "warning"),  # Different patch versions
    ],
)
def test_check_version_compatibility(
    client, server, expected_compatible, expected_severity
):
    """Test version compatibility checking."""
    compatible, severity, message = check_version_compatibility(client, server)
    assert compatible == expected_compatible
    assert severity == expected_severity
    assert isinstance(message, str)
    assert len(message) > 0


def test_check_version_compatibility_invalid():
    """Test version compatibility with invalid version strings."""
    compatible, severity, message = check_version_compatibility("invalid", "0.6.0")
    assert not compatible
    assert severity == "error"
    assert "Invalid version format" in message


def test_validate_server_version_compatible():
    """Test validate_server_version with compatible versions."""
    mock_api = Mock()
    mock_api.get_version.return_value = "0.6.0"
    mock_api.url = "http://localhost:5000"

    # Should not raise
    validate_server_version(mock_api, "0.6.0")


def test_validate_server_version_patch_mismatch():
    """Test validate_server_version with patch version mismatch (warning)."""
    mock_api = Mock()
    mock_api.get_version.return_value = "0.6.1"
    mock_api.url = "http://localhost:5000"

    # Should warn but not raise
    with pytest.warns(UserWarning, match="Patch version mismatch"):
        validate_server_version(mock_api, "0.6.0")


def test_validate_server_version_minor_mismatch():
    """Test validate_server_version with minor version mismatch (error)."""
    mock_api = Mock()
    mock_api.get_version.return_value = "0.7.0"
    mock_api.url = "http://localhost:5000"

    # Should raise RuntimeError
    with pytest.raises(RuntimeError, match="Minor version mismatch"):
        validate_server_version(mock_api, "0.6.0")


def test_validate_server_version_major_mismatch():
    """Test validate_server_version with major version mismatch (error)."""
    mock_api = Mock()
    mock_api.get_version.return_value = "1.0.0"
    mock_api.url = "http://localhost:5000"

    # Should raise RuntimeError
    with pytest.raises(RuntimeError, match="Major version mismatch"):
        validate_server_version(mock_api, "0.6.0")


def test_validate_server_version_network_error():
    """Test validate_server_version handles network errors gracefully."""
    mock_api = Mock()
    mock_api.get.side_effect = Exception("Network error")
    mock_api.url = "http://localhost:5000"

    # Should not raise - just log warning
    validate_server_version(mock_api, "0.6.0")
