from zndraw.geometries.sphere import Sphere
from zndraw.geometries.arrow import Arrow

geometries = {
    "Sphere": Sphere,
    "Arrow": Arrow,
}

__all__ = ["Sphere", "Arrow", "geometries"]
