"""Base geometry class for all ZnDraw geometries."""

from typing import Literal, Union
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field

# Type alias for data properties that can be dynamic or static
DataProp = Union[
    str, float, tuple[float, float, float], list[tuple[float, float, float]]
]

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

class BaseGeometry(BaseModel):
    """Base class for all geometries with common properties."""

    model_config = ConfigDict(frozen=True)

    active: bool = Field(
        default=True,
        description="Whether this geometry should be rendered. Inactive geometries are hidden.",
    )

    position: DataProp = Field(
        default="arrays.positions",
        description="Position [x,y,z]. String for dynamic data key, tuple/list for static values.",
    )

    color: DataProp = Field(
        default="arrays.colors",
        description="Color [r,g,b]. String for dynamic data key, tuple/list for static values.",
    )

    material: Material = Field(
        default=Material.MeshPhysicalMaterial_matt,
        description="Material type (static config, not fetched from server)",
    )
