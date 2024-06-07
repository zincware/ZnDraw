import dataclasses
import os

import znsocket
from redis import Redis

KEY: str = "config"


def _get_client_from_url(url: str) -> Redis | znsocket.Client:
    if url.startswith("redis"):
        return Redis.from_url(url, decode_responses=True)
    elif url.startswith("znsocket"):
        return znsocket.Client.from_url(url)
    else:
        raise ValueError(f"Unknown storage type: {url}")


@dataclasses.dataclass
class Config:
    port: int
    address: str = "http://127.0.0.1"

    def __post_init__(self):
        storage = os.environ["ZNDRAW_STORAGE"]
        self._r = _get_client_from_url(storage)

    @classmethod
    def from_env(cls) -> "Config":
        """Load the Config using the `ZNDRAW_STORAGE` environment variable."""
        storage = os.environ["ZNDRAW_STORAGE"]
        r = _get_client_from_url(storage)

        shared = znsocket.Dict(r, KEY)
        return cls(**shared)

    def save(self) -> None:
        shared = znsocket.Dict(self._r, KEY)
        shared.update(dataclasses.asdict(self))

    def get_url(self) -> str:
        return f"{self.address}:{self.port}"
