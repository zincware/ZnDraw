import typing as t

from pydantic import Field

from .base import BaseGeometry, apply_schema_feature


class Cell(BaseGeometry):
    """A geometry to draw the simulation cell."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(schema, "color", ["color-picker", "editable-array"])
        # TODO: add dropdown with "default" option, where default means follow color scheme

        return schema

    color: str = Field(default="default", description="Color of the cell")

    position: str = Field(
        default="cell",
        description="Cell vectors. Should be a 3x3 matrix.",
    )

    thickness: float = Field(
        default=2.0,
        ge=0.0,
        description="Line thickness for the cell vectors.",
    )

    material: t.Literal[
        "LineBasicMaterial",
        "LineDashedMaterial",
    ] = Field(
        default="LineBasicMaterial",
        description="Material type (static config, not fetched from server)",
    )
