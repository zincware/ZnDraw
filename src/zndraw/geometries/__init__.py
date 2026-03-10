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
from zndraw.geometries.circle_curve import CircleCurve
from zndraw.geometries.curve import Curve, CurveMarker
from zndraw.geometries.floor import Floor
from zndraw.geometries.fog import Fog
from zndraw.geometries.lights import (
    AmbientLight,
    DirectionalLight,
    HemisphereLight,
    LightPosition,
)
from zndraw.geometries.pathtracing import EnvironmentPreset, PathTracing
from zndraw.geometries.plane import Plane
from zndraw.geometries.property_inspector import PropertyInspector
from zndraw.geometries.shape import Shape
from zndraw.geometries.sphere import Sphere
from zndraw.transformations import InArrayTransform, Transform

Sphere.model_rebuild()
Arrow.model_rebuild()
Bond.model_rebuild()
Box.model_rebuild()
Cell.model_rebuild()
CircleCurve.model_rebuild()
Curve.model_rebuild()
Floor.model_rebuild()
Plane.model_rebuild()
Shape.model_rebuild()


geometries = {
    "Sphere": Sphere,
    "Arrow": Arrow,
    "Bond": Bond,
    "Curve": Curve,
    "CircleCurve": CircleCurve,
    "Cell": Cell,
    "Floor": Floor,
    "Camera": Camera,
    "Box": Box,
    "Plane": Plane,
    "Shape": Shape,
    # New scene object types
    "DirectionalLight": DirectionalLight,
    "AmbientLight": AmbientLight,
    "HemisphereLight": HemisphereLight,
    "Fog": Fog,
    "PathTracing": PathTracing,
    "PropertyInspector": PropertyInspector,
}

__all__ = [
    "AmbientLight",
    "Arrow",
    "Bond",
    "Box",
    "Camera",
    "CameraType",
    "Cell",
    "CircleCurve",
    "ColorProp",
    "Curve",
    "CurveMarker",
    "DirectionalLight",
    "EnvironmentPreset",
    "Floor",
    "Fog",
    "HemisphereLight",
    "InArrayTransform",
    "InteractionSettings",
    "LightPosition",
    "PathTracing",
    "Plane",
    "PositionProp",
    "PropertyInspector",
    "RotationProp",
    "Shape",
    "SizeProp",
    "Sphere",
    "Transform",
    "geometries",
]
