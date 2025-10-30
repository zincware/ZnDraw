from zndraw.geometries.arrow import Arrow
from zndraw.geometries.base import (
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    SizeProp,
)
from zndraw.geometries.bonds import Bond
from zndraw.geometries.box import Box
from zndraw.geometries.camera import Camera, CameraType
from zndraw.geometries.cell import Cell
from zndraw.geometries.curve import Curve, CurveMarker
from zndraw.geometries.floor import Floor
from zndraw.geometries.plane import Plane
from zndraw.geometries.sphere import Sphere

# Rebuild Pydantic models after Transform is fully defined
# This resolves forward references to Transform in type hints
try:
    from zndraw.transformations import Transform

    Sphere.model_rebuild()
    Arrow.model_rebuild()
    Bond.model_rebuild()
    Box.model_rebuild()
    Cell.model_rebuild()
    Curve.model_rebuild()
    Floor.model_rebuild()
    Plane.model_rebuild()
except ImportError:
    # Transform not available yet during initial imports
    pass

geometries = {
    "Sphere": Sphere,
    "Arrow": Arrow,
    "Bond": Bond,
    "Curve": Curve,
    "Cell": Cell,
    "Floor": Floor,
    "Camera": Camera,
    "Box": Box,
    "Plane": Plane,
}

__all__ = [
    "Sphere",
    "Arrow",
    "Bond",
    "Curve",
    "geometries",
    "InteractionSettings",
    "CurveMarker",
    "Cell",
    "Floor",
    "Camera",
    "CameraType",
    "Box",
    "Plane",
    "PositionProp",
    "ColorProp",
    "SizeProp",
    "RotationProp",
]
