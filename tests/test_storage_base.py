"""Tests for StorageBackend abstract base class."""

import inspect
from abc import ABC

import pytest

from zndraw.storage import StorageBackend


def test_backend_cannot_instantiate_abc() -> None:
    """StorageBackend cannot be instantiated directly."""
    with pytest.raises(TypeError, match="abstract"):
        StorageBackend()  # type: ignore[abstract]


def test_backend_is_abstract_base_class() -> None:
    """StorageBackend is an ABC."""
    assert issubclass(StorageBackend, ABC)


def test_backend_has_required_methods() -> None:
    """StorageBackend defines required async methods."""
    assert hasattr(StorageBackend, "get")
    assert hasattr(StorageBackend, "get_range")
    assert hasattr(StorageBackend, "get_many")
    assert hasattr(StorageBackend, "extend")
    assert hasattr(StorageBackend, "set_item")
    assert hasattr(StorageBackend, "merge_item")
    assert hasattr(StorageBackend, "get_length")
    assert hasattr(StorageBackend, "delete_range")
    assert hasattr(StorageBackend, "clear")
    assert hasattr(StorageBackend, "close")


def test_backend_methods_are_abstract() -> None:
    """All required methods are abstract."""
    abstract_methods = StorageBackend.__abstractmethods__
    expected_methods = {
        "get",
        "get_range",
        "get_many",
        "extend",
        "set_item",
        "merge_item",
        "get_length",
        "delete_range",
        "clear",
        "reserve",
        "remove_items",
        "close",
    }
    assert abstract_methods == expected_methods


def test_backend_methods_are_coroutines() -> None:
    """All methods are async (coroutine functions)."""
    methods = [
        "get",
        "get_range",
        "get_many",
        "extend",
        "set_item",
        "merge_item",
        "get_length",
        "delete_range",
        "clear",
        "reserve",
        "remove_items",
        "close",
    ]
    for method_name in methods:
        method = getattr(StorageBackend, method_name)
        assert inspect.iscoroutinefunction(method), f"{method_name} should be async"


def test_backend_incomplete_implementation_raises() -> None:
    """Implementing class without all methods raises TypeError."""

    class IncompleteStorage(StorageBackend):
        async def get(self, room_id: str, index: int) -> dict | None:
            return None

    with pytest.raises(TypeError, match="abstract"):
        IncompleteStorage()  # type: ignore[abstract]


def test_backend_complete_implementation_works() -> None:
    """Implementing class with all methods can be instantiated."""
    from typing import Any

    class CompleteStorage(StorageBackend):
        async def get(self, room_id: str, index: int) -> dict[str, Any] | None:
            return None

        async def get_range(
            self, room_id: str, start: int, stop: int | None
        ) -> list[dict[str, Any]]:
            return []

        async def get_many(
            self, room_id: str, indices: list[int]
        ) -> list[dict[str, Any]]:
            return []

        async def extend(self, room_id: str, frames: list[dict[str, Any]]) -> int:
            return 0

        async def set_item(
            self, room_id: str, index: int, frame: dict[str, Any]
        ) -> None:
            pass

        async def merge_item(
            self, room_id: str, index: int, partial: dict[bytes, bytes]
        ) -> None:
            pass

        async def get_length(self, room_id: str) -> int:
            return 0

        async def delete_range(self, room_id: str, start: int, stop: int) -> None:
            pass

        async def clear(self, room_id: str) -> None:
            pass

        async def reserve(self, room_id: str, count: int) -> None:
            pass

        async def remove_items(self, room_id: str, indices: list[int]) -> None:
            pass

        async def close(self) -> None:
            pass

    # Should not raise
    storage = CompleteStorage()
    assert isinstance(storage, StorageBackend)
