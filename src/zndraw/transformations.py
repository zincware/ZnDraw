"""Data transformation models for filtering geometry properties.

Transformations allow dynamic filtering of geometry data based on frame data.
For example, filtering positions to only show atoms that are constrained.

Also includes CurveAttachment for referencing Curve geometries with progress.
"""

from typing import Literal

from pydantic import BaseModel, Field


class CurveAttachment(BaseModel):
    """Reference to a Curve geometry with progress along the curve.

    Used for camera position/target that follows a curve path.
    Similar to InArrayTransform, this allows referencing external geometry data.

    Parameters
    ----------
    type : str
        Must be "curve_attachment" for this type.
    geometry_key : str
        Key of the Curve geometry in vis.geometries.
    progress : float
        Progress along the curve (0.0 to 1.0).

    Examples
    --------
    Attach camera position to a flight path:

    >>> CurveAttachment(
    ...     geometry_key="flight_path",
    ...     progress=0.5
    ... )
    """

    type: Literal["curve_attachment"] = "curve_attachment"
    geometry_key: str = Field(description="Key of Curve geometry in vis.geometries")
    progress: float = Field(
        default=0.0, ge=0.0, le=1.0, description="Progress along curve (0.0 to 1.0)"
    )


def is_curve_attachment(value) -> bool:
    """Check if value is a CurveAttachment dict.

    Parameters
    ----------
    value : Any
        Value to check.

    Returns
    -------
    bool
        True if value is a dict with type="curve_attachment".
    """
    return isinstance(value, dict) and value.get("type") == "curve_attachment"


class InArrayTransform(BaseModel):
    """Filter array data to only include elements at specified indices.

    This transform fetches indices from one source and uses them to filter
    data from another source. Commonly used for constraint visualization.

    Parameters
    ----------
    type : str
        Must be "in_array" for this transform type.
    source : str
        Frame data key containing the indices to extract (e.g., "constraints").
    path : str
        Dot-separated path to extract indices from source data
        (e.g., "FixAtoms.indices" to get constraints["FixAtoms"]["indices"]).
    filter : str
        Frame data key to filter using the extracted indices
        (e.g., "arrays.positions" to filter atom positions).

    Examples
    --------
    Filter positions to only show fixed atoms:

    >>> InArrayTransform(
    ...     source="constraints",
    ...     path="FixAtoms.indices",
    ...     filter="arrays.positions"
    ... )

    Filter radii to only show fixed atoms:

    >>> InArrayTransform(
    ...     source="constraints",
    ...     path="FixAtoms.indices",
    ...     filter="arrays.radii"
    ... )
    """

    type: Literal["in_array"] = "in_array"
    source: str = Field(description="Frame data key containing indices")
    path: str = Field(description="Dot-separated path to extract indices")
    filter: str = Field(description="Frame data key to filter")


# Union type for all transforms (currently just InArrayTransform)
Transform = InArrayTransform
