"""Base geometry class for all ZnDraw geometries."""

from typing import Literal, Union, Any
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

# Type aliases for geometry properties
# Position is ALWAYS per-instance (list of tuples or dynamic key)
PositionProp = Union[str, list[tuple[float, float, float]]]

# Color can be shared (single hex) or per-instance (list of hex) or dynamic
ColorProp = Union[str, list[str]]

# Size/radius can be shared (single value) or per-instance (list) or dynamic
SizeProp = Union[str, float, list[float]]

# Rotation can be shared (single tuple) or per-instance (list of tuples) or dynamic
RotationProp = Union[str, tuple[float, float, float], list[tuple[float, float, float]]]

# Generic data property (for backward compatibility)
DataProp = Union[str, float, tuple[float, float, float], list[tuple[float, float, float]]]

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
        default="arrays.positions",
        description="Position coordinates. String for dynamic data key (e.g. 'arrays.positions'), list of tuples for static per-instance positions [(x,y,z), ...].",
    )

    color: ColorProp = Field(
        default="arrays.colors",
        description="Color values. String for dynamic key (e.g. 'arrays.colors') or shared hex color (e.g. '#FF0000'), list of hex colors for per-instance ['#FF0000', '#00FF00', ...].",
    )

    material: Material = Field(
        default=Material.MeshPhysicalMaterial_matt,
        description="Material type (static config, not fetched from server)",
    )
