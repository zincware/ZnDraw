"""Data transformation models for filtering geometry properties.

Transformations allow dynamic filtering of geometry data based on frame data.
For example, filtering positions to only show atoms that are constrained.
"""

from typing import Literal

from pydantic import BaseModel, Field


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
