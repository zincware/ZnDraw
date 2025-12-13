"""Base geometry class for all ZnDraw geometries."""

import math
from typing import TYPE_CHECKING, Any, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator

if TYPE_CHECKING:
    from zndraw.transformations import Transform

from zndraw.materials import MaterialProp

# Type aliases for geometry properties
# Position is ALWAYS per-instance (list of tuples or dynamic key or transform)
PositionProp = Union[str, list[tuple[float, float, float]], "Transform"]

# Color is always hex strings - single hex or list of hex strings, or dynamic reference, or transform
ColorProp = Union[str, list[str], "Transform"]

# Size/radius can be shared (single value) or per-instance (list) or dynamic, or transform
SizeProp = Union[str, float, list[float], "Transform"]

# 2D size for plane (width, height)
Size2DProp = Union[str, tuple[float, float], list[tuple[float, float]]]

# 3D size for box (width, height, depth)
Size3DProp = Union[str, tuple[float, float, float], list[tuple[float, float, float]]]

# Rotation can be shared (single tuple) or per-instance (list of tuples) or dynamic
RotationProp = Union[str, tuple[float, float, float], list[tuple[float, float, float]]]

# Scale can be uniform (float), per-instance uniform (list[float]),
# per-instance anisotropic (list[tuple]), or dynamic (str)
ScaleProp = Union[
    str,  # Dynamic reference
    float,  # Uniform for all instances (isotropic)
    list[float],  # Per-instance uniform (isotropic)
    tuple[float, float, float],  # Shared anisotropic
    list[tuple[float, float, float]],  # Per-instance anisotropic
]

# Connectivity for bonds: per-bond indices and properties (atom_a, atom_b, bond_order, etc)
# Bond order can be None (defaults to 1), int, or float (e.g., 1.5 for aromatic)
ConnectivityProp = Union[str, list[tuple[float, float, float | None]]]


# Validation constants
SCALE_MIN = 1e-6  # Minimum scale to prevent zero/negative scales
SCALE_MAX = 1e6  # Maximum scale to prevent WebGL issues
ROTATION_MAX = 1000.0  # Maximum rotation in radians (reasonable upper bound)


def validate_numeric_value(value: float, name: str = "value") -> float:
    """Validate a numeric value for NaN, Infinity, and reasonable bounds.

    Parameters
    ----------
    value
        The numeric value to validate
    name
        The name of the value (for error messages)

    Returns
    -------
    float
        The validated value

    Raises
    ------
    ValueError
        If the value is NaN, Infinity, or outside reasonable bounds
    """
    if math.isnan(value):
        raise ValueError(f"{name} cannot be NaN")
    if math.isinf(value):
        raise ValueError(f"{name} cannot be Infinity")
    return value


def validate_scale_value(value: float) -> float:
    """Validate a scale value.

    Parameters
    ----------
    value
        The scale value to validate

    Returns
    -------
    float
        The validated scale value

    Raises
    ------
    ValueError
        If the scale is NaN, Infinity, or outside bounds [SCALE_MIN, SCALE_MAX]
    """
    value = validate_numeric_value(value, "scale")
    if value < SCALE_MIN:
        raise ValueError(f"scale must be >= {SCALE_MIN}, got {value}")
    if value > SCALE_MAX:
        raise ValueError(f"scale must be <= {SCALE_MAX}, got {value}")
    return value


def validate_rotation_value(value: float) -> float:
    """Validate a rotation value.

    Parameters
    ----------
    value
        The rotation value (in radians) to validate

    Returns
    -------
    float
        The validated rotation value

    Raises
    ------
    ValueError
        If the rotation is NaN, Infinity, or too large
    """
    value = validate_numeric_value(value, "rotation")
    if abs(value) > ROTATION_MAX:
        raise ValueError(
            f"rotation must be in range [-{ROTATION_MAX}, {ROTATION_MAX}], got {value}"
        )
    return value


class InteractionSettings(BaseModel):
    model_config = ConfigDict(frozen=True)

    enabled: bool = Field(default=True)
    color: str = Field(default="#FF6600")
    opacity: float = Field(default=1.0, ge=0.0, le=1.0)


def apply_schema_feature(
    schema: dict[str, Any],
    property_name: str,
    features: list[str],
    force_string_type: bool = True,
    definition_path: str | None = None,
) -> None:
    """Apply x-custom-type and x-features to a property in the schema.

    This helper consolidates common schema manipulation patterns used across
    geometry types for UI customization (dropdowns, color pickers, etc.).

    Parameters
    ----------
    schema
        The JSON schema dictionary to modify
    property_name
        Name of the property to customize
    features
        List of feature flags (e.g. ["color-picker", "dynamic-atom-props"])
    force_string_type
        If True, forces property type to "string" and removes "anyOf"
    definition_path
        Optional path to nested definition (e.g. "InteractionSettings")
    """
    if definition_path:
        target = schema["$defs"][definition_path]["properties"][property_name]
    else:
        target = schema["properties"][property_name]

    target["x-custom-type"] = "dynamic-enum"
    target["x-features"] = features

    if force_string_type:
        target["type"] = "string"
        target.pop("anyOf", None)


class BaseGeometry(BaseModel):
    """Base class for all geometries with common properties."""

    model_config = ConfigDict(frozen=True)

    active: bool = Field(
        default=True,
        description="Whether this geometry should be rendered. Inactive geometries are hidden.",
    )

    position: PositionProp = Field(
        default=[(0, 0, 0)],
        description="Position coordinates. String for dynamic data key (e.g. 'arrays.positions'), list of tuples for static per-instance positions [(x,y,z), ...].",
    )

    color: ColorProp = Field(
        default=["#FFA200"],
        description="Color values. String for dynamic key (e.g. 'arrays.colors') or list of hex colors ['#FF0000', '#00FF00', ...]. Single hex colors are automatically normalized to lists.",
    )

    material: MaterialProp = Field(
        default="MeshPhysicalMaterial_matt",
        description="Material type or object. String for preset (e.g. 'MeshPhysicalMaterial_matt') or material object for full customization.",
    )

    @field_validator("position", mode="before")
    @classmethod
    def normalize_position(cls, v):
        """Normalize position to list of tuples."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, dict):  # Transform object
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v  # Already a list or Transform instance

    @field_validator("color", mode="before")
    @classmethod
    def normalize_color(cls, v):
        """Normalize color to list of hex strings."""
        if v is None:
            return []
        # Dynamic reference -> pass through (strings not starting with #)
        if isinstance(v, str) and not v.startswith("#"):
            return v
        # Transform object -> pass through
        if isinstance(v, dict):
            return v
        # Single hex color -> wrap in list
        if isinstance(v, str) and v.startswith("#"):
            return [v]
        # Already a list or Transform instance -> return as is
        return v

    @classmethod
    def model_json_schema(cls, **kwargs: Any) -> dict[str, Any]:
        """Customize JSON schema for UI rendering.

        Applies x-custom-type: "three-material" to the material field
        to enable the custom material editor in JSON Forms.
        """
        schema = super().model_json_schema(**kwargs)

        # Apply three-material custom type to material field
        if "material" in schema["properties"]:
            schema["properties"]["material"]["x-custom-type"] = "three-material"

        return schema
