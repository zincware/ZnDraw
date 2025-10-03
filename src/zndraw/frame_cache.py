import typing as t

from cachetools import LRUCache


class FrameCache:
    def __init__(self, maxsize=100):
        self._cache = LRUCache(maxsize=maxsize)

    def keys(self):
        return self._cache.keys()

    def __len__(self):
        return len(self._cache)

    def get(self, key: int) -> t.Any:
        return self._cache.get(key)

    def set(self, key: int, value: t.Any):
        self._cache[key] = value

    def __contains__(self, key: int) -> bool:
        return key in self._cache

    def pop(self, key: int, default: t.Any = None) -> t.Any:
        return self._cache.pop(key, default)

    def invalidate_from(self, from_idx: int):
        new_cache = LRUCache(maxsize=self._cache.maxsize)
        for k, v in self._cache.items():
            if k < from_idx:
                new_cache[k] = v
        self._cache = new_cache

    def clear(self):
        self._cache.clear()

    @property
    def maxsize(self) -> int:
        return int(self._cache.maxsize)

    def items(self):
        return self._cache.items()
