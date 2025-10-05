from zndraw.geometries.arrow import Arrow
from zndraw.geometries.sphere import Sphere

geometries = {
    "Sphere": Sphere,
    "Arrow": Arrow,
}

__all__ = ["Sphere", "Arrow", "geometries"]
