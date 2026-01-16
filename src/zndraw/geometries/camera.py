"""Camera geometry with flexible position/target specification.

Cameras can use either direct coordinates or CurveAttachment for position/target.
For session cameras (viewport state), direct coordinates are typical.
For cinematic cameras (in vis.geometries), CurveAttachment enables path animation.
"""

import typing as t
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, field_validator

from zndraw.transformations import CurveAttachment


class CameraType(str, Enum):
    """Camera projection types."""

    PERSPECTIVE = "PerspectiveCamera"
    ORTHOGRAPHIC = "OrthographicCamera"


# Union type for position/target - either direct coordinates or curve reference
PositionType = tuple[float, float, float] | CurveAttachment


class Camera(BaseModel):
    """A camera geometry with flexible position/target specification.

    Used for both geometry cameras (stored in vis.geometries) and the interactive
    camera (stored in user settings). The only difference is that interactive
    camera has helper_visible=False by default.

    Position and target can be either:
    - Direct coordinates: (10, 5, 10) - for static cameras or session viewport state
    - CurveAttachment: reference to a Curve geometry for path animation

    Parameters
    ----------
    position : tuple[float, float, float] | CurveAttachment
        Camera position, either as direct coordinates or a CurveAttachment.
    target : tuple[float, float, float] | CurveAttachment
        Camera look-at target, either as direct coordinates or a CurveAttachment.
    up : tuple[float, float, float]
        Camera up vector [x,y,z]. Defaults to [0, 1, 0] (Y-axis up).
    camera_type : CameraType
        Type of camera projection (PERSPECTIVE or ORTHOGRAPHIC).
    fov : float
        Field of view in degrees (perspective only). Range: 1-179 degrees.
    near : float
        Near clipping plane distance. Must be positive.
    far : float
        Far clipping plane distance. Must be greater than near.
    zoom : float
        Camera zoom factor. Default 1.0 (no zoom).
    helper_visible : bool
        Whether to show the camera helper (cone visualization) in the scene.
    helper_color : str
        Color of the camera helper (hex or named color).

    Examples
    --------
    Create a static camera with direct coordinates:

    >>> Camera(
    ...     position=(-10, 10, 30),
    ...     target=(0, 0, 0),
    ...     fov=60
    ... )

    Create an animated camera using CurveAttachment:

    >>> from zndraw.transformations import CurveAttachment
    >>> vis.geometries["flight_path"] = Curve(position=[
    ...     [0, 0, 10], [5, 5, 10], [10, 0, 10]
    ... ])
    >>> Camera(
    ...     position=CurveAttachment(geometry_key="flight_path", progress=0.5),
    ...     target=(0, 0, 0),  # Can mix direct and CurveAttachment
    ...     fov=60
    ... )

    Animate camera by updating progress:

    >>> from zndraw.transformations import CurveAttachment
    >>> cam = vis.geometries["camera"]
    >>> # Update position to follow curve at 80%
    >>> cam_copy = cam.model_copy()
    >>> cam_copy.position = CurveAttachment(geometry_key="flight_path", progress=0.8)
    >>> vis.geometries["camera"] = cam_copy
    """

    model_config = ConfigDict(frozen=True)

    active: bool = Field(
        default=True,
        description="Whether this geometry should be rendered.",
    )

    protected: bool = Field(
        default=False,
        description="Whether this geometry is protected from deletion.",
    )

    # Position and target with union type
    position: PositionType = Field(
        default=(-10.0, 10.0, 30.0),
        description="Camera position - direct coordinates or CurveAttachment",
        json_schema_extra={
            "x-custom-type": "position-attachment",
            "x-features": ["dynamic-geometries"],
            "x-geometry-filter": "Curve",
        },
    )

    target: PositionType = Field(
        default=(0.0, 0.0, 0.0),
        description="Camera look-at target - direct coordinates or CurveAttachment",
        json_schema_extra={
            "x-custom-type": "position-attachment",
            "x-features": ["dynamic-geometries"],
            "x-geometry-filter": "Curve",
        },
    )

    # Camera orientation
    up: tuple[float, float, float] = Field(
        default=(0.0, 1.0, 0.0),
        description="Camera up vector [x,y,z].",
    )

    # Camera type and projection parameters
    camera_type: CameraType = Field(
        default=CameraType.PERSPECTIVE, description="Type of camera projection"
    )

    fov: float = Field(
        default=50.0,
        description="Field of view in degrees (perspective only)",
        ge=1,
        le=179,
        json_schema_extra={"format": "range", "step": 1},
    )

    near: float = Field(
        default=0.1,
        description="Near clipping plane",
        ge=0.01,
        le=100,
        json_schema_extra={"format": "range", "step": 0.01},
    )

    far: float = Field(
        default=1000.0,
        description="Far clipping plane",
        ge=1,
        le=10000,
        json_schema_extra={"format": "range", "step": 1},
    )

    zoom: float = Field(
        default=1.0,
        description="Camera zoom factor",
        ge=0.1,
        le=10,
        json_schema_extra={"format": "range", "step": 0.1},
    )

    # Visual helper settings
    helper_visible: bool = Field(
        default=False, description="Show camera helper (cone visualization)"
    )

    helper_color: str = Field(
        default="#00ff00",
        description="Color of the camera helper",
        json_schema_extra={
            "format": "color",
            "x-custom-type": "color-picker",
            "x-features": ["color-picker"],
        },
    )

    # Rendering settings
    show_crosshair: bool = Field(
        default=False, description="Show a crosshair at the screen center"
    )

    preserve_drawing_buffer: bool = Field(
        default=False,
        description="Enable screenshot capture (WARNING: reduces rendering performance)",
    )

    # Override material and color (not applicable for cameras)
    material: t.Any = Field(default=None, description="Not applicable for cameras")

    color: str = Field(
        default="#00ff00",
        description="Camera helper color (visualization only)",
        json_schema_extra={
            "format": "color",
            "x-custom-type": "color-picker",
            "x-features": ["color-picker"],
        },
    )

    @field_validator("far")
    @classmethod
    def validate_far_greater_than_near(cls, v, info):
        """Ensure far plane is greater than near plane."""
        if "near" in info.data and v <= info.data["near"]:
            raise ValueError(
                f"far ({v}) must be greater than near ({info.data['near']})"
            )
        return v

    @field_validator("up")
    @classmethod
    def validate_up_vector(cls, v):
        """Ensure up vector is not zero."""
        if abs(v[0]) < 1e-10 and abs(v[1]) < 1e-10 and abs(v[2]) < 1e-10:
            raise ValueError("up vector cannot be zero")
        return v
