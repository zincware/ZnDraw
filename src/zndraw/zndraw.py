import dataclasses
from zndraw.client import Client
from zndraw.extensions import Extension, ExtensionType
import typing as t
import requests




@dataclasses.dataclass
class ZnDraw:
    url: str
    room: str
    user: str
    auto_pickup_jobs: bool = True

    def __post_init__(self):
        self.client = Client(url=self.url, room=self.room, user=self.user, auto_pickup_jobs=self.auto_pickup_jobs)
        self.client.connect()

    @property
    def step(self) -> int:
        return self.client.step

    @step.setter
    def step(self, value: int):
        self.client.step = value

    def __len__(self) -> int:
        return len(self.client)

    @property
    def settings(self):
        return self.client.settings

    def register_extension(self, *args, **kwargs):
        self.client.register_extension(*args, **kwargs)
