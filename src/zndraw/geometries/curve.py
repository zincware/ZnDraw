from pydantic import Field, BaseModel
import typing as t

from .base import BaseGeometry, DataProp

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
    opacity: float = Field(
        default=1.0,
        description="Opacity of the markers",
        ge=0,
        le=1,
    )


class Curve(BaseGeometry):
    """A curve geometry."""

    position: DataProp = Field(
        default=[],
        description="Position [x,y,z]. String for dynamic data key, tuple/list for static values.",
    )

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
    virtual_marker: CurveMarker | None = Field(
        default_factory=lambda: CurveMarker(size=0.08, color="gray", opacity=0.5),
        description="Virtual marker between two existing markers (for adding new points)",
    )
    color: str = Field(
        default="#001722",
        description="Curve color",
    )