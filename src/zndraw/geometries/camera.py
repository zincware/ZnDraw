"""Camera geometry for controlling viewport with curve attachment support.

Cameras always reference curves for position and target.
For static cameras, create single-point curves.
Moving curve markers with TransformControls updates camera immediately.
"""

import typing as t
from enum import Enum

from pydantic import BaseModel, ConfigDict, Field, field_validator

from .base import BaseGeometry


class CameraType(str, Enum):
    """Camera projection types."""

    PERSPECTIVE = "PerspectiveCamera"
    ORTHOGRAPHIC = "OrthographicCamera"


class Camera(BaseModel):
    """A camera geometry for controlling viewport using curve references.

    Cameras ALWAYS reference curves for position and target.
    For static cameras, create single-point curves.
    Moving curve markers with TransformControls updates camera immediately.

    Parameters
    ----------
    position_curve_key : str
        Geometry key of the curve defining camera position path.
    position_progress : float
        Position along the position curve (0.0 to 1.0).
    target_curve_key : str
        Geometry key of the curve defining camera target (look-at) path.
    target_progress : float
        Position along the target curve (0.0 to 1.0).
    up : tuple[float, float, float]
        Camera up vector [x,y,z]. Defines which direction is "up" for the camera.
        Defaults to [0, 1, 0] (Y-axis up).
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
    Create a static camera (using single-point curves):

    >>> vis.geometries["cam_pos"] = Curve(position=[[0, 0, 10]])
    >>> vis.geometries["cam_target"] = Curve(position=[[0, 0, 0]])
    >>> vis.geometries["camera"] = Camera(
    ...     position_curve_key="cam_pos",
    ...     target_curve_key="cam_target",
    ... )

    Create an animated flythrough camera:

    >>> vis.geometries["flight_path"] = Curve(position=[
    ...     [0, 0, 10], [5, 5, 10], [10, 0, 10]
    ... ])
    >>> vis.geometries["focus_path"] = Curve(position=[
    ...     [0, 0, 0], [2, 2, 0], [4, 0, 0]
    ... ])
    >>> vis.geometries["camera"] = Camera(
    ...     position_curve_key="flight_path",
    ...     position_progress=0.5,  # Halfway along path
    ...     target_curve_key="focus_path",
    ...     target_progress=0.5,
    ...     fov=60
    ... )

    Animate camera by updating progress:

    >>> camera = vis.geometries["camera"]
    >>> dump = camera.model_dump()
    >>> dump["position_progress"] = 0.8
    >>> vis.geometries["camera"] = Camera.model_validate(dump)
    """

    model_config = ConfigDict(frozen=True)

    active: bool = Field(
        default=True,
        description="Whether this geometry should be rendered. Inactive geometries are hidden.",
    )

    # Curve references (required)
    position_curve_key: str = Field(
        description="Geometry key of curve for camera position"
    )

    position_progress: float = Field(
        default=0.0,
        description="Position along position curve (0.0 to 1.0)",
        ge=0.0,
        le=1.0,
    )

    target_curve_key: str = Field(
        description="Geometry key of curve for camera target (look-at point)"
    )

    target_progress: float = Field(
        default=0.0,
        description="Position along target curve (0.0 to 1.0)",
        ge=0.0,
        le=1.0,
    )

    # Override material and color (not applicable for cameras)
    material: t.Any = Field(default=None, description="Not applicable for cameras")

    color: str = Field(
        default="#00ff00", description="Camera helper color (visualization only)"
    )

    # Camera orientation
    up: tuple[float, float, float] = Field(
        default=(0.0, 1.0, 0.0),
        description="Camera up vector [x,y,z]. Defines which direction is 'up'. Should be normalized.",
    )

    # Camera type and projection parameters
    camera_type: CameraType = Field(
        default=CameraType.PERSPECTIVE, description="Type of camera projection"
    )

    fov: float = Field(
        default=75.0,
        description="Field of view in degrees (perspective only)",
        gt=0,
        lt=180,
    )

    near: float = Field(
        default=0.1,
        description="Near clipping plane",
        gt=0,
    )

    far: float = Field(
        default=1000.0,
        description="Far clipping plane",
        gt=0,
    )

    zoom: float = Field(
        default=1.0,
        description="Camera zoom factor",
        gt=0,
    )

    # Visual helper settings (no size field - managed by Three.js)
    helper_visible: bool = Field(
        default=True, description="Show camera helper (cone visualization)"
    )

    helper_color: str = Field(
        default="#00ff00", description="Color of the camera helper (hex or named color)"
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

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        """Generate JSON schema with custom UI hints."""
        schema = super().model_json_schema(**kwargs)

        # Position curve key - dropdown of available curves
        if "position_curve_key" in schema["properties"]:
            schema["properties"]["position_curve_key"]["x-custom-type"] = "dynamic-enum"
            schema["properties"]["position_curve_key"]["x-features"] = [
                "dynamic-geometries"
            ]
            schema["properties"]["position_curve_key"]["x-geometry-filter"] = "Curve"

        # Target curve key - dropdown of available curves
        if "target_curve_key" in schema["properties"]:
            schema["properties"]["target_curve_key"]["x-custom-type"] = "dynamic-enum"
            schema["properties"]["target_curve_key"]["x-features"] = [
                "dynamic-geometries"
            ]
            schema["properties"]["target_curve_key"]["x-geometry-filter"] = "Curve"

        # Helper color customization
        if "helper_color" in schema["properties"]:
            schema["properties"]["helper_color"]["x-custom-type"] = "color-picker"
            schema["properties"]["helper_color"]["x-features"] = ["color-picker"]

        # Position progress - slider control
        if "position_progress" in schema["properties"]:
            schema["properties"]["position_progress"]["format"] = "range"
            schema["properties"]["position_progress"]["step"] = 0.01

        # Target progress - slider control
        if "target_progress" in schema["properties"]:
            schema["properties"]["target_progress"]["format"] = "range"
            schema["properties"]["target_progress"]["step"] = 0.01

        return schema
