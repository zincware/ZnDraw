"""Unit tests for FilesystemRead provider."""

from pathlib import Path
from unittest.mock import MagicMock

import fsspec
import pytest

from zndraw.providers.filesystem import FileItem, FilesystemRead, _from_info


@pytest.fixture
def nested_fs(tmp_path: Path):
    """Real fsspec local filesystem with nested .xyz files.

    Structure::

        tmp_path/
            top.xyz
            sub/
                s22.xyz
                other.txt
                deep/
                    nested.xyz
    """
    (tmp_path / "top.xyz").touch()
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "s22.xyz").touch()
    (sub / "other.txt").touch()
    deep = sub / "deep"
    deep.mkdir()
    (deep / "nested.xyz").touch()

    return fsspec.filesystem("file"), str(tmp_path)


def test_filesystem_read_ls():
    """ls(detail=True) returns info dicts directly."""
    fs = MagicMock()
    fs.ls.return_value = [
        {"name": "/data/a.xyz", "size": 100, "type": "file"},
        {"name": "/data/b.xyz", "size": 200, "type": "file"},
    ]

    provider = FilesystemRead(path="/data")
    result = provider.read(fs)

    fs.ls.assert_called_once_with("/data", detail=True)
    assert len(result) == 2
    assert result[0]["name"] == "a.xyz"
    assert result[0]["type"] == "file"


def test_filesystem_read_glob():
    fs = MagicMock()
    fs.glob.return_value = ["/data/a.xyz", "/data/sub/b.xyz"]
    fs.info.side_effect = lambda p: {"name": p, "size": 50, "type": "file"}

    provider = FilesystemRead(path="/data", glob="**/*.xyz")
    result = provider.read(fs)

    fs.glob.assert_called_once_with("/data/**/*.xyz")
    assert len(result) == 2


def test_filesystem_read_directory_type():
    fs = MagicMock()
    fs.ls.return_value = [{"name": "/data/subdir", "size": 0, "type": "directory"}]

    provider = FilesystemRead(path="/data")
    result = provider.read(fs)

    assert result[0]["type"] == "directory"


def test_filesystem_read_category():
    assert FilesystemRead.category == "filesystem"


def test_filesystem_read_schema():
    schema = FilesystemRead.model_json_schema()
    assert "path" in schema["properties"]
    assert "glob" in schema["properties"]


def test_from_info():
    """_from_info normalizes an fsspec info dict."""
    info = {"name": "/some/path/file.xyz", "size": 42, "type": "file"}
    result = _from_info(info)
    assert result == FileItem(
        name="file.xyz",
        path="/some/path/file.xyz",
        size=42,
        type="file",
    )


def test_from_info_missing_size():
    """Missing size defaults to 0."""
    info = {"name": "/a/b.txt", "type": "file"}
    assert _from_info(info).size == 0


def test_from_info_none_size():
    """None size defaults to 0."""
    info = {"name": "/a/b.txt", "size": None, "type": "file"}
    assert _from_info(info).size == 0


# --- RED tests: real fsspec glob behavior ---


def test_glob_recursive_finds_all_xyz(nested_fs):
    """**/*.xyz from base should find all .xyz files recursively."""
    fs, base = nested_fs
    provider = FilesystemRead(path=base, glob="**/*.xyz")
    result = provider.read(fs)

    names = {r["name"] for r in result}
    assert names == {"top.xyz", "s22.xyz", "nested.xyz"}


def test_glob_recursive_finds_specific_file(nested_fs):
    """**/s22.xyz should find s22.xyz in any subdirectory."""
    fs, base = nested_fs
    provider = FilesystemRead(path=base, glob="**/s22.xyz")
    result = provider.read(fs)

    assert len(result) == 1
    assert result[0]["name"] == "s22.xyz"


def test_glob_scoped_to_subdirectory(nested_fs):
    """Glob from a subdirectory should not find files above it."""
    fs, base = nested_fs
    provider = FilesystemRead(path=f"{base}/sub", glob="**/*.xyz")
    result = provider.read(fs)

    names = {r["name"] for r in result}
    assert "s22.xyz" in names
    assert "nested.xyz" in names
    assert "top.xyz" not in names


def test_glob_non_recursive_only_matches_immediate(nested_fs):
    """*.xyz (no **) should only match files in the immediate directory."""
    fs, base = nested_fs
    provider = FilesystemRead(path=base, glob="*.xyz")
    result = provider.read(fs)

    names = {r["name"] for r in result}
    assert names == {"top.xyz"}


def test_glob_result_paths_are_absolute(nested_fs):
    """Glob results should have full absolute paths usable for opening."""
    fs, base = nested_fs
    provider = FilesystemRead(path=base, glob="**/*.xyz")
    result = provider.read(fs)

    for item in result:
        assert item["path"].startswith(base), (
            f"Path {item['path']!r} should start with base {base!r}"
        )
        assert fs.exists(item["path"]), f"Path {item['path']!r} should be a valid path"


def test_glob_memory_filesystem():
    """Glob works with fsspec memory filesystem."""
    fs = fsspec.filesystem("memory")
    fs.pipe("/data/top.xyz", b"content")
    fs.pipe("/data/sub/s22.xyz", b"content")
    fs.pipe("/data/sub/other.txt", b"content")

    provider = FilesystemRead(path="/data", glob="**/*.xyz")
    result = provider.read(fs)

    names = {r["name"] for r in result}
    assert "top.xyz" in names
    assert "s22.xyz" in names
    assert "other.txt" not in names


def test_glob_memory_filesystem_specific_file():
    """Find a specific file by name across subdirectories."""
    fs = fsspec.filesystem("memory")
    fs.pipe("/data/sub/s22.xyz", b"content")
    fs.pipe("/data/deep/nested/s22.xyz", b"content")

    provider = FilesystemRead(path="/data", glob="**/s22.xyz")
    result = provider.read(fs)

    assert len(result) == 2
    assert all(r["name"] == "s22.xyz" for r in result)


# --- RED tests: root path glob bug (relative filesystem) ---


@pytest.fixture
def relative_fs(tmp_path: Path):
    """DirFileSystem rooted at tmp_path — paths are relative.

    This mimics how register_filesystem() works in practice:
    the fsspec handler is relative to a base directory, so ls("/")
    returns items with relative paths (no leading slash).
    """
    (tmp_path / "top.xyz").touch()
    sub = tmp_path / "sub"
    sub.mkdir()
    (sub / "s22.xyz").touch()
    (sub / "other.txt").touch()
    deep = sub / "deep"
    deep.mkdir()
    (deep / "nested.xyz").touch()

    from fsspec.implementations.dirfs import DirFileSystem

    return DirFileSystem(path=str(tmp_path), fs=fsspec.filesystem("file"))


def test_glob_from_root_slash_with_relative_fs(relative_fs):
    """**/*.xyz from path='/' should find files in a relative filesystem.

    This is the actual user scenario: the filesystem browser starts at '/',
    the user searches for '**/*.xyz', but the pattern becomes '/**/*.xyz'
    which doesn't match relative paths like 'sub/s22.xyz'.
    """
    provider = FilesystemRead(path="/", glob="**/*.xyz")
    result = provider.read(relative_fs)

    names = {r["name"] for r in result}
    assert "top.xyz" in names
    assert "s22.xyz" in names
    assert "nested.xyz" in names


def test_glob_specific_file_from_root_slash(relative_fs):
    """**/s22.xyz from path='/' should find s22.xyz in subdirectories."""
    provider = FilesystemRead(path="/", glob="**/s22.xyz")
    result = provider.read(relative_fs)

    assert len(result) == 1
    assert result[0]["name"] == "s22.xyz"


def test_glob_from_subdirectory_with_relative_fs(relative_fs):
    """Glob from a relative subdirectory should work."""
    provider = FilesystemRead(path="sub", glob="**/*.xyz")
    result = provider.read(relative_fs)

    names = {r["name"] for r in result}
    assert "s22.xyz" in names
    assert "nested.xyz" in names
    assert "top.xyz" not in names


def test_ls_root_slash_returns_relative_paths(relative_fs):
    """ls('/') on a DirFileSystem returns relative paths (no leading /)."""
    provider = FilesystemRead(path="/")
    result = provider.read(relative_fs)

    # Verify that we actually get results (the fs is browseable from /)
    names = {r["name"] for r in result}
    assert "top.xyz" in names
    assert "sub" in names
