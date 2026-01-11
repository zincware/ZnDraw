"""Base geometry class for all ZnDraw geometries."""

import typing as t

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

if t.TYPE_CHECKING:
    from zndraw.transformations import Transform

from zndraw.materials import MaterialProp

# Type aliases for geometry properties
Vec3 = tuple[float, float, float]
Vec2 = tuple[float, float]

# Position: always per-instance list of 3D coordinates
PositionProp = t.Union[str, list[Vec3], "Transform"]

# Color: always per-instance list of hex strings
ColorProp = t.Union[str, list[str], "Transform"]

# Size/radius: always per-instance list of floats
SizeProp = t.Union[str, list[float], "Transform"]

# 2D size for plane (width, height): always per-instance list
Size2DProp = t.Union[str, list[Vec2]]

# 3D size for box (width, height, depth): always per-instance list
Size3DProp = t.Union[str, list[Vec3]]

# Rotation: always per-instance list of Euler angles (radians)
RotationProp = t.Union[str, list[Vec3]]

# Scale: always per-instance list of 3D scale factors
ScaleProp = t.Union[str, list[Vec3]]

# Connectivity for bonds: per-bond indices and properties
ConnectivityProp = t.Union[str, list[tuple[float, float, float | None]]]


def coerce_numpy_arrays(data: dict[str, t.Any]) -> dict[str, t.Any]:
    """Convert numpy arrays in data dict to lists for Pydantic validation."""
    result = {}
    for key, value in data.items():
        if isinstance(value, np.ndarray):
            result[key] = value.tolist()
        else:
            result[key] = value
    return result


class InteractionSettings(BaseModel):
    """Settings for geometry interaction (selection, hover)."""

    model_config = ConfigDict(frozen=True)

    enabled: bool = Field(default=True)
    color: str = Field(
        default="#FF6600",
        json_schema_extra={
            "format": "color",
            "x-custom-type": "dynamic-enum",
            "x-features": ["color-picker", "dynamic-atom-props", "free-solo"],
        },
    )
    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        json_schema_extra={"format": "range", "step": 0.01},
    )


class BaseGeometry(BaseModel):
    """Base class for all geometries with common properties."""

    model_config = ConfigDict(frozen=True)

    @model_validator(mode="before")
    @classmethod
    def coerce_numpy(cls, data: t.Any) -> t.Any:
        """Convert numpy arrays to lists before Pydantic validation."""
        if isinstance(data, dict):
            return coerce_numpy_arrays(data)
        return data

    active: bool = Field(
        default=True,
        description="Whether this geometry should be rendered.",
    )

    protected: bool = Field(
        default=False,
        description="Whether this geometry is protected from deletion.",
    )

    position: PositionProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Position coordinates [(x,y,z), ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array", "transform"],
        },
    )

    color: ColorProp = Field(
        default=["#FFA200"],
        description="Color values ['#RRGGBB', ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": [
                "color-picker",
                "dynamic-atom-props",
                "free-solo",
                "editable-array",
            ],
        },
    )

    material: MaterialProp = Field(
        default="MeshPhysicalMaterial_matt",
        description="Material type or object.",
        json_schema_extra={"x-custom-type": "three-material"},
    )
