"""Tests for global settings API endpoint."""

import requests


def test_global_settings_simgen_disabled_by_default(server):
    """Test that simgen is disabled by default when no CLI flag is provided."""
    response = requests.get(f"{server}/api/config/global-settings")
    assert response.status_code == 200

    data = response.json()
    assert "simgen" in data
    assert data["simgen"]["enabled"] is False


def test_global_settings_endpoint_structure(server):
    """Test that the endpoint returns the expected structure."""
    response = requests.get(f"{server}/api/config/global-settings")
    assert response.status_code == 200

    data = response.json()

    # Verify the structure
    assert isinstance(data, dict)
    assert "simgen" in data
    assert isinstance(data["simgen"], dict)
    assert "enabled" in data["simgen"]
    assert isinstance(data["simgen"]["enabled"], bool)
