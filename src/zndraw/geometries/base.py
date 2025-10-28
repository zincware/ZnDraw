"""Base geometry class for all ZnDraw geometries."""

from enum import Enum
from typing import Any, Literal, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator

# Type aliases for geometry properties
# Position is ALWAYS per-instance (list of tuples or dynamic key)
PositionProp = Union[str, list[tuple[float, float, float]]]

# Color is always hex strings - single hex or list of hex strings, or dynamic reference
ColorProp = Union[str, list[str]]

# Size/radius can be shared (single value) or per-instance (list) or dynamic
SizeProp = Union[str, float, list[float]]

# 2D size for plane (width, height)
Size2DProp = Union[str, tuple[float, float], list[tuple[float, float]]]

# 3D size for box (width, height, depth)
Size3DProp = Union[str, tuple[float, float, float], list[tuple[float, float, float]]]

# Rotation can be shared (single tuple) or per-instance (list of tuples) or dynamic
RotationProp = Union[str, tuple[float, float, float], list[tuple[float, float, float]]]

# Connectivity for bonds: per-bond indices and properties (atom_a, atom_b, bond_order, etc)
# Bond order can be None (defaults to 1), int, or float (e.g., 1.5 for aromatic)
ConnectivityProp = Union[str, list[tuple[float, float, float | None]]]


class Material(str, Enum):
    # --- PHYSICAL MATERIALS ---
    MeshPhysicalMaterial_matt = "MeshPhysicalMaterial (matt)"
    MeshPhysicalMaterial_semi_gloss = "MeshPhysicalMaterial (semi-gloss)"
    MeshPhysicalMaterial_shiny = "MeshPhysicalMaterial (shiny)"
    MeshPhysicalMaterial_transparent = "MeshPhysicalMaterial (transparent)"
    MeshPhysicalMaterial_glass = "MeshPhysicalMaterial (glass)"

    # --- STANDARD MATERIALS ---
    MeshStandardMaterial_matt = "MeshStandardMaterial (matt)"
    MeshStandardMaterial_metallic = "MeshStandardMaterial (metallic)"

    # --- SIMPLE / STYLISED MATERIALS ---
    MeshBasicMaterial = "MeshBasicMaterial"
    MeshToonMaterial = "MeshToonMaterial"
    MeshLambertMaterial_matt = "MeshLambertMaterial (matt)"
    MeshPhongMaterial_classic = "MeshPhongMaterial (classic)"


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

    material: Material = Field(
        default=Material.MeshPhysicalMaterial_matt,
        description="Material type (static config, not fetched from server)",
    )

    @field_validator("position", mode="before")
    @classmethod
    def normalize_position(cls, v):
        """Normalize position to list of tuples."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v  # Already a list

    @field_validator("color", mode="before")
    @classmethod
    def normalize_color(cls, v):
        """Normalize color to list of hex strings."""
        if v is None:
            return []
        # Dynamic reference -> pass through (strings not starting with #)
        if isinstance(v, str) and not v.startswith("#"):
            return v
        # Single hex color -> wrap in list
        if isinstance(v, str) and v.startswith("#"):
            return [v]
        # Already a list -> return as is
        return v
