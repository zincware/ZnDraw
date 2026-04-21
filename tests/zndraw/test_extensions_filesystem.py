"""Unit tests for LoadFile extension."""

from pathlib import Path

import fsspec
import pytest
from fsspec.implementations.dirfs import DirFileSystem

from zndraw.extensions.filesystem import LoadFile, _require_local_path


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


@pytest.mark.parametrize(
    "attack",
    [
        "/../etc/passwd",
        "/a/../../etc/passwd",
        "/./../x",
        "/subdir/../../outside",
    ],
)
def test_require_local_path_rejects_traversal(tmp_path: Path, attack: str) -> None:
    fs = DirFileSystem(path=str(tmp_path), fs=fsspec.filesystem("file"))
    with pytest.raises(PermissionError, match="escapes provider root"):
        _require_local_path(fs, attack)


def test_require_local_path_accepts_subpath(tmp_path: Path) -> None:
    fs = DirFileSystem(path=str(tmp_path), fs=fsspec.filesystem("file"))
    resolved = _require_local_path(fs, "/sub/file.xyz")
    assert Path(resolved) == (tmp_path / "sub" / "file.xyz").resolve()
