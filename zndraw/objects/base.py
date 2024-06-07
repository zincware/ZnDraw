"""Base class for all ZnDraw objects."""

import dataclasses


@dataclasses.dataclass
class Material:
    color: str = "#ffffff"
    opacity: float = 1
    wireframe: bool = False


@dataclasses.dataclass
class Object3D:
    position: list[float] | tuple = dataclasses.field(default=(0, 0, 0))
    scale: list[float] | tuple = dataclasses.field(default=(1, 1, 1))
    material: Material = dataclasses.field(default_factory=Material)


@dataclasses.dataclass
class Box(Object3D):
    width: float = 1
    height: float = 1
    depth: float = 1
