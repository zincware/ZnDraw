"""Cell geometry for ZnDraw."""

import typing as t

from pydantic import Field

from .base import BaseGeometry


class Cell(BaseGeometry):
    """A geometry to draw the simulation cell."""

    position: str = Field(
        default="cell",
        description="Cell vectors. Should be a 3x3 matrix.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    color: str = Field(
        default="default",
        description="Color of the cell",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["color-picker", "editable-array"],
        },
    )

    thickness: float = Field(
        default=2.0,
        ge=0.0,
        le=10,
        description="Line thickness for the cell vectors.",
        json_schema_extra={"format": "range", "step": 0.5},
    )

    material: t.Literal[
        "LineBasicMaterial",
        "LineDashedMaterial",
    ] = Field(
        default="LineBasicMaterial",
        description="Material type (static config, not fetched from server)",
    )
