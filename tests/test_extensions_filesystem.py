"""Unit tests for LoadFile extension."""

import pytest

from zndraw.extensions.filesystem import LoadFile


def test_load_file_schema():
    schema = LoadFile.model_json_schema()
    assert "provider_name" in schema["properties"]
    assert "path" in schema["properties"]
    assert "start" in schema["properties"]
    assert "stop" in schema["properties"]
    assert "step" in schema["properties"]


def test_load_file_requires_provider():
    ext = LoadFile(provider_name="missing", path="/data/file.xyz")
    with pytest.raises(ValueError, match="not found"):
        ext.run(None, providers={})


def test_load_file_category():
    from zndraw_joblib.client import Category

    assert LoadFile.category == Category.MODIFIER
