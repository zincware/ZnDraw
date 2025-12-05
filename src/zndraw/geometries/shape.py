"""Shape geometry for ZnDraw - 2D polygons rendered in 3D space."""

import typing as t

from pydantic import Field, field_validator

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    apply_schema_feature,
)


class Shape(BaseGeometry):
    """A 2D polygon shape rendered in 3D space.

    Vertices define the shape outline in 2D (XY plane).
    Position, rotation, scale place instances in 3D space.
    Supports instancing: one shape definition, multiple placements.

    Note: Vertex editing is not supported because vertices are 2D
    but transform controls operate in 3D.

    Example
    -------
    >>> from zndraw import ZnDraw, Shape
    >>> vis = ZnDraw()
    >>> vis.geometries["triangle"] = Shape(
    ...     vertices=[(0, 0), (1, 0), (0.5, 0.866)],
    ...     position=[(0, 0, 0)],
    ...     color=["#FF0000"],
    ... )
    """

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        # Note: vertices intentionally NOT editable-array - 2D tuples not supported by UI
        apply_schema_feature(
            schema, "rotation", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
        )

        return schema

    # Shape outline vertices (2D) - THE SHAPE TEMPLATE
    vertices: list[tuple[float, float]] = Field(
        default=[(0, 0), (1, 0), (0.5, 0.866)],  # Default equilateral triangle
        description="Shape vertices in 2D [(x,y), ...]. Connected by lines, auto-closed.",
        json_schema_extra={"x-custom-type": "vertices-2d"},
    )

    # Instance positions (where to place shape copies)
    position: PositionProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Instance center positions in 3D [(x,y,z), ...]",
    )

    # Instance rotations
    rotation: RotationProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Instance rotations as Euler angles [x, y, z] in radians",
    )

    # Scale (uniform)
    scale: float = Field(
        default=1.0,
        ge=0.0,
        description="Uniform scale factor for all instances",
    )

    # Appearance
    color: ColorProp = Field(
        default=["#808080"],
        description="Fill color(s). Single color or per-instance list.",
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Shape opacity",
    )

    # Interaction
    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )
    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )

    @field_validator("vertices", mode="before")
    @classmethod
    def normalize_vertices(cls, v):
        """Ensure vertices is a list of tuples."""
        if v is None:
            return []
        return v

    @field_validator("rotation", mode="before")
    @classmethod
    def normalize_rotation(cls, v):
        """Normalize rotation to list of tuples."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v
