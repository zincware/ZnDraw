"""Tests for Presets accessor class.

Tests the Presets(MutableMapping[str, Preset]) accessor by mocking
APIManager to verify mapping interface and convenience methods.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from zndraw.accessors import Presets
from zndraw.schemas import Preset, PresetApplyResult


@pytest.fixture
def mock_api() -> MagicMock:
    """Create a mock APIManager with preset methods."""
    return MagicMock()


@pytest.fixture
def presets(mock_api: MagicMock) -> Presets:
    """Create a Presets accessor with mock API."""
    return Presets(mock_api)


# =============================================================================
# Mapping interface
# =============================================================================


def test_getitem_returns_preset(presets: Presets, mock_api: MagicMock) -> None:
    """__getitem__ returns a Preset from API data."""
    mock_api.get_preset.return_value = {
        "name": "matt",
        "description": "Matt finish",
        "rules": [{"pattern": "fog", "config": {"active": True}}],
    }

    result = presets["matt"]

    assert isinstance(result, Preset)
    assert result.name == "matt"
    assert len(result.rules) == 1
    mock_api.get_preset.assert_called_once_with("matt")


def test_getitem_raises_keyerror_when_not_found(
    presets: Presets, mock_api: MagicMock
) -> None:
    """__getitem__ raises KeyError when API returns None."""
    mock_api.get_preset.return_value = None

    with pytest.raises(KeyError):
        presets["nonexistent"]


def test_setitem_calls_update(presets: Presets, mock_api: MagicMock) -> None:
    """__setitem__ calls APIManager.update_preset."""
    preset = Preset(name="custom", rules=[])
    presets["custom"] = preset

    mock_api.update_preset.assert_called_once_with("custom", preset.model_dump())


def test_delitem_calls_delete(presets: Presets, mock_api: MagicMock) -> None:
    """__delitem__ calls APIManager.delete_preset."""
    del presets["custom"]

    mock_api.delete_preset.assert_called_once_with("custom")


def test_iter_returns_names(presets: Presets, mock_api: MagicMock) -> None:
    """__iter__ yields preset names from list."""
    mock_api.list_presets.return_value = [
        {"name": "matt"},
        {"name": "flat"},
    ]

    names = list(presets)
    assert names == ["matt", "flat"]


def test_len_returns_count(presets: Presets, mock_api: MagicMock) -> None:
    """__len__ returns number of presets."""
    mock_api.list_presets.return_value = [{"name": "a"}, {"name": "b"}, {"name": "c"}]

    assert len(presets) == 3


# =============================================================================
# Convenience methods
# =============================================================================


def test_apply_returns_result(presets: Presets, mock_api: MagicMock) -> None:
    """apply() returns a PresetApplyResult."""
    mock_api.apply_preset.return_value = {"geometries_updated": ["fog", "particles"]}

    result = presets.apply("matt")

    assert isinstance(result, PresetApplyResult)
    assert result.geometries_updated == ["fog", "particles"]
    mock_api.apply_preset.assert_called_once_with("matt")


def test_load_reads_json_and_creates(presets: Presets, mock_api: MagicMock) -> None:
    """load() reads JSON from file and calls create_preset."""
    preset_data = {
        "name": "loaded",
        "description": "From file",
        "rules": [{"pattern": "fog", "config": {"active": True}}],
    }
    mock_api.create_preset.return_value = preset_data

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(preset_data, f)
        f.flush()
        result = presets.load(Path(f.name))

    assert isinstance(result, Preset)
    assert result.name == "loaded"
    mock_api.create_preset.assert_called_once()


def test_export_writes_json_file(presets: Presets, mock_api: MagicMock) -> None:
    """export() writes preset JSON to file."""
    mock_api.get_preset.return_value = {
        "name": "matt",
        "description": "Matt finish",
        "rules": [],
        "created_at": "2026-01-01T00:00:00",
        "updated_at": "2026-01-01T00:00:00",
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    presets.export("matt", output_path)

    data = json.loads(output_path.read_text())
    assert data["name"] == "matt"
    # Timestamps should be stripped
    assert "created_at" not in data
    assert "updated_at" not in data
