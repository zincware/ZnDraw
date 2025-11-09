import typing as t

import numpy as np
import splines
from pydantic import BaseModel, Field

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    apply_schema_feature,
)


class CurveMarker(BaseModel):
    """Settings for markers on the curve control points."""

    enabled: bool = Field(
        default=True,
        description="Whether to show markers at the control points",
    )
    size: float = Field(
        default=0.1,
        description="Size of the markers",
        gt=0,
    )
    color: str = Field(
        default="default",
        description="Color of the markers. If None, uses the curve color",
    )
    opacity: float = Field(
        default=1.0,
        description="Opacity of the markers",
        ge=0,
        le=1,
    )
    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings for markers.",
    )
    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings for markers.",
    )


class Curve(BaseGeometry):
    """A curve geometry."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
        )

        # Color picker for CurveMarker properties
        if "CurveMarker" in schema.get("$defs", {}):
            schema["$defs"]["CurveMarker"]["properties"]["color"]["x-custom-type"] = (
                "dynamic-enum"
            )
            schema["$defs"]["CurveMarker"]["properties"]["color"]["x-features"] = [
                "color-picker",
                "dynamic-atom-props",
                "free-solo",
                "editable-array",
            ]

        return schema

    position: PositionProp = Field(
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
        default=2.0,
        description="Thickness of the line (not implemented in Three.js LineBasicMaterial)",
        gt=0,
    )
    marker: CurveMarker = Field(
        default_factory=CurveMarker,
        description="Settings for markers at the control points",
    )
    virtual_marker: CurveMarker = Field(
        default_factory=lambda: CurveMarker(size=0.08, color="default", opacity=0.5),
        description="Virtual marker between two existing markers (for adding new points)",
    )
    color: ColorProp = Field(
        default="default",
        description="Curve color",
    )

    def get_interpolated_points(self) -> np.ndarray:
        """Get interpolated points along the curve using CatmullRom spline.

        Returns
        -------
        np.ndarray
            Array of shape (n_points, 3) with interpolated points along the curve.
            If less than 2 control points, returns the control points as-is.
        """
        points = np.array(self.position)
        if points.shape[0] <= 1:
            return points
        t = np.linspace(0, len(points) - 1, len(points) * self.divisions)
        return splines.CatmullRom(points).evaluate(t)
