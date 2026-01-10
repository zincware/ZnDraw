"""Curve geometry for ZnDraw."""

import typing as t

import numpy as np
import splines
from pydantic import BaseModel, Field

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
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
        ge=0.01,
        le=1,
        json_schema_extra={"format": "range", "step": 0.01},
    )
    color: str = Field(
        default="default",
        description="Color of the markers. If 'default', uses the curve color.",
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
    opacity: float = Field(
        default=1.0,
        description="Opacity of the markers",
        ge=0,
        le=1,
        json_schema_extra={"format": "range", "step": 0.01},
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

    position: PositionProp = Field(
        default_factory=list,
        description="Control points: list of [x,y,z]. String for dynamic data key, list of tuples/lists for static values.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    color: ColorProp = Field(
        default="default",
        description="Curve color",
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
        le=200,
        json_schema_extra={"format": "range", "step": 1},
    )

    thickness: float = Field(
        default=2.0,
        description="Thickness of the line (not implemented in Three.js LineBasicMaterial)",
        ge=0.5,
        le=10,
        json_schema_extra={"format": "range", "step": 0.5},
    )

    marker: CurveMarker = Field(
        default_factory=CurveMarker,
        description="Settings for markers at the control points",
    )

    virtual_marker: CurveMarker = Field(
        default_factory=lambda: CurveMarker(size=0.08, color="default", opacity=0.5),
        description="Virtual marker between two existing markers (for adding new points)",
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
