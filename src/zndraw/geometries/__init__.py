from zndraw.geometries.arrow import Arrow
from zndraw.geometries.sphere import Sphere
from zndraw.geometries.base import InteractionSettings

geometries = {
    "Sphere": Sphere,
    "Arrow": Arrow,
}

__all__ = ["Sphere", "Arrow", "geometries", "InteractionSettings"]
