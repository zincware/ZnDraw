import dataclasses


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None
    sphere_size: float = None
    bond_size: float = None
    max_fps: int = None


# TODO set defaults here and load in typer?


config = Config()

graph = None
