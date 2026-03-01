"""CircleCurve geometry for ZnDraw."""

import math
import typing as t

import numpy as np
from pydantic import Field
from scipy.spatial.transform import Rotation

from .base import BaseGeometry, Vec3


class CircleCurve(BaseGeometry):
    """A circular/elliptical curve for camera orbit paths.

    The curve lies in the XY plane by default. Use ``rotation`` to orient it
    in 3D and ``scale`` (X/Y) to create an ellipse.  Angles are given as
    a percentage of a full circle (0-100).
    """

    position: Vec3 = Field(
        default=(0.0, 0.0, 0.0),
        description="Center point of the circle.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    radius: float = Field(
        default=5.0,
        ge=0.1,
        le=100.0,
        description="Base radius of the circle.",
        json_schema_extra={"format": "range", "step": 0.1},
    )

    start_angle: float = Field(
        default=0.0,
        ge=0.0,
        le=100.0,
        description="Start angle as percentage of full circle (0-100).",
        json_schema_extra={"format": "range", "step": 1.0},
    )

    end_angle: float = Field(
        default=100.0,
        ge=0.0,
        le=100.0,
        description="End angle as percentage of full circle (0-100).",
        json_schema_extra={"format": "range", "step": 1.0},
    )

    rotation: Vec3 = Field(
        default=(0.0, 0.0, 0.0),
        description="Euler angles [x, y, z] in radians to orient the circle plane.",
        json_schema_extra={"x-features": ["editable-array"]},
    )

    scale: Vec3 = Field(
        default=(1.0, 1.0, 1.0),
        description="Scale factors. X/Y create ellipse, Z has no visual effect.",
        json_schema_extra={"x-features": ["editable-array"]},
    )

    color: str = Field(
        default="#FFA200",
        description="Curve color.",
        json_schema_extra={"format": "color"},
    )

    material: t.Literal["LineBasicMaterial", "LineDashedMaterial"] = Field(
        default="LineBasicMaterial",
        description="Line material type.",
    )

    divisions: int = Field(
        default=50,
        ge=1,
        le=200,
        description="Number of segments to approximate the curve.",
        json_schema_extra={"format": "range", "step": 1},
    )

    thickness: float = Field(
        default=2.0,
        ge=0.5,
        le=10.0,
        description="Line thickness.",
        json_schema_extra={"format": "range", "step": 0.5},
    )

    def get_interpolated_points(self) -> np.ndarray:
        """Get interpolated points along the circle/ellipse.

        Returns
        -------
        np.ndarray
            Array of shape (divisions + 1, 3) with points along the curve.
            Points are in world space (center + rotation + scale applied).
        """
        start_rad = (self.start_angle / 100.0) * 2 * math.pi
        end_rad = (self.end_angle / 100.0) * 2 * math.pi

        t = np.linspace(start_rad, end_rad, self.divisions + 1)

        # Generate points in XY plane
        points = np.column_stack(
            [
                self.radius * np.cos(t),
                self.radius * np.sin(t),
                np.zeros_like(t),
            ]
        )

        # Apply scale (X/Y for ellipse, Z is a no-op for a flat curve)
        sx, sy, sz = self.scale
        points[:, 0] *= sx
        points[:, 1] *= sy
        points[:, 2] *= sz

        # Apply rotation
        rx, ry, rz = self.rotation
        if rx != 0.0 or ry != 0.0 or rz != 0.0:
            rot = Rotation.from_euler("xyz", [rx, ry, rz])
            points = rot.apply(points)

        # Apply translation
        points += np.array(self.position)

        return points
