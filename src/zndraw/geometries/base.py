"""Base geometry class for all ZnDraw geometries."""

from typing import Literal, Union

from pydantic import BaseModel, ConfigDict, Field

# Type alias for data properties that can be dynamic or static
DataProp = Union[
    str, float, tuple[float, float, float], list[tuple[float, float, float]]
]

class InteractionSettings(BaseModel):
    color: str | None = Field(None)
    opacity: float = Field(1.0, ge=0.0, le=1.0)

class BaseGeometry(BaseModel):
    """Base class for all geometries with common properties."""

    model_config = ConfigDict(frozen=True)

    position: DataProp = Field(
        default="arrays.positions",
        description="Position [x,y,z]. String for dynamic data key, tuple/list for static values.",
    )

    color: DataProp = Field(
        default="arrays.colors",
        description="Color [r,g,b]. String for dynamic data key, tuple/list for static values.",
    )

    material: Literal[
        "MeshPhysicalMaterial",
        "MeshStandardMaterial",
        "MeshBasicMaterial",
        "MeshToonMaterial",
    ] = Field(
        default="MeshStandardMaterial",
        description="Material type (static config, not fetched from server)",
    )
