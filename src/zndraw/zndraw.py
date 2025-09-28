import dataclasses
from zndraw.client import Client

@dataclasses.dataclass
class ZnDraw:
    url: str
    room: str

    def __post_init__(self):
        self.client = Client(self.url, self.room)
        self.client.connect()

    @property
    def step(self) -> int|None:
        return self.client.step

    @step.setter
    def step(self, value: int):
        self.client.step = value

    def __len__(self) -> int:
        return len(self.client)

