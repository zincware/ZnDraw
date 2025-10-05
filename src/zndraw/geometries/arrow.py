"""Arrow geometry for ZnDraw."""

from pydantic import Field

from .base import BaseGeometry, DataProp


class Arrow(BaseGeometry):
    """Arrow geometry with direction vector.

    Arrows are defined by a start position and a direction vector.
    The direction vector determines both the arrow's orientation and its length.
    """

    # Use 'start' as the primary field with 'position' alias for compatibility
    start: DataProp = Field(
        default="arrays.positions",
        description="Arrow start position [x,y,z]. String for dynamic data key, tuple for static value.",
        alias="position",
    )

    direction: DataProp = Field(
        default="arrays.forces",
        description="Direction vector [x,y,z]. Defines arrow orientation and base length. String for dynamic data key, tuple for static value.",
    )

    radius: DataProp = Field(
        default=0.05,
        description="Arrow shaft radius. String for dynamic data key, float for static value.",
    )

    scale: DataProp = Field(
        default=1.0,
        description="Length scale multiplier applied to direction vector length. String for dynamic data key, float for static value.",
    )
