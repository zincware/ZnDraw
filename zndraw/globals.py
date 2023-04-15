import dataclasses


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None
    sphere_size: float = 1.0
    bond_size: float = 1.0


config = Config()

graph = None
