from pydantic import Field, BaseModel
import typing as t

from .base import BaseGeometry

class CurveMarker(BaseModel):
    """Settings for markers on the curve control points."""

    size: float = Field(
        default=0.1,
        description="Size of the markers",
        gt=0,
    )
    color: t.Optional[str] = Field(
        default=None,
        description="Color of the markers. If None, uses the curve color",
    )


class Curve(BaseGeometry):
    """A curve geometry."""

    material: t.Literal[
        "LineBasicMaterial",
        "LineDashedMaterial",
    ] = Field(
        default="LineBasicMaterial",
        description="Material type (static config, not fetched from server)",
    )
    variant: t.Literal["CatmullRomCurve3"] = Field(
        default="CatmullRomCurve3",
        description="Curve variant type (static config, not fetched from server)",
    )
    divisions: int = Field(
        default=50,
        description="Number of divisions along the curve",
        ge=1,
    )
    thickness: float = Field(
        default=1.0,
        description="Thickness of the line (not implemented in Three.js LineBasicMaterial)",
        gt=0,
    )
    marker: CurveMarker | None = Field(
        default_factory=CurveMarker,
        description="Settings for markers at the control points",
    )