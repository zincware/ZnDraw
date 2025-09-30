import dataclasses
from zndraw.client import Client


@dataclasses.dataclass
class ZnDraw:
    url: str
    room: str
    user: str

    def __post_init__(self):
        self.client = Client(url=self.url, room=self.room, user=self.user)
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

