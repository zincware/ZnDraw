import dataclasses


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None


config = Config()

graph = None
