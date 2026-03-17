"""Isosurface geometry for volumetric data visualization."""

from typing import Self

from pydantic import BaseModel, ConfigDict, Field, model_validator


class Isosurface(BaseModel):
    """An isosurface extracted from a 3D volumetric grid.

    References a key in ``atoms.info`` that contains a dict with:
    - ``grid``: 3D float array (Nx, Ny, Nz) of scalar values
    - ``origin``: 3-vector, world-space origin of the grid
    - ``cell``: (3, 3) matrix, axis vectors spanning the grid

    The frontend fetches mesh data from a dedicated endpoint that runs
    marching cubes server-side.
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the geometry owner. None means unowned.",
        json_schema_extra={"x-custom-type": "ownership-toggle"},
    )

    cube_key: str = Field(
        default="",
        description=(
            "Frame info key for the volumetric data dict. "
            "Must contain grid, origin, and cell entries."
        ),
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props"],
        },
    )

    isovalue: float = Field(
        default=0.02,
        description="Scalar threshold for surface extraction.",
        json_schema_extra={
            "x-custom-type": "editable-range",
            "step": 0.001,
            "x-min-field": "isovalue_min",
            "x-max-field": "isovalue_max",
        },
    )

    isovalue_min: float = Field(
        default=-0.25,
        description="Minimum isovalue for the slider.",
        json_schema_extra={"x-hidden": True},
    )

    isovalue_max: float = Field(
        default=0.25,
        description="Maximum isovalue for the slider.",
        json_schema_extra={"x-hidden": True},
    )

    sigma: float = Field(
        default=0.0,
        ge=0.0,
        le=5.0,
        description="Gaussian smoothing sigma. 0 = disabled.",
        json_schema_extra={"format": "range", "step": 0.1},
    )

    color: str = Field(
        default="#2244CC",
        description="Surface color.",
        json_schema_extra={"format": "color"},
    )

    resolution: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Mesh resolution. 0 = coarse (fast), 1 = fine (slow).",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    opacity: float = Field(
        default=0.6,
        ge=0.0,
        le=1.0,
        description="Surface opacity.",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    active: bool = Field(default=True, description="Show or hide this isosurface.")

    @model_validator(mode="after")
    def _check_isovalue_range(self) -> Self:
        if self.isovalue_min > self.isovalue_max:
            raise ValueError(
                f"isovalue_min ({self.isovalue_min}) must be <= "
                f"isovalue_max ({self.isovalue_max})"
            )
        if not (self.isovalue_min <= self.isovalue <= self.isovalue_max):
            raise ValueError(
                f"isovalue ({self.isovalue}) must be within "
                f"[{self.isovalue_min}, {self.isovalue_max}]"
            )
        return self
