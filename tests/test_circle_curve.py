"""Tests for CircleCurve geometry model."""

import math

import numpy as np
import pytest
from pydantic import ValidationError

from zndraw.geometries.circle_curve import CircleCurve


def test_circle_curve_defaults():
    """CircleCurve has sensible defaults."""
    c = CircleCurve()
    assert c.position == (0.0, 0.0, 0.0)
    assert c.radius == 5.0
    assert c.start_angle == 0.0
    assert c.end_angle == 100.0
    assert c.rotation == (0.0, 0.0, 0.0)
    assert c.scale == (1.0, 1.0, 1.0)
    assert c.color == "#FFA200"
    assert c.material == "LineBasicMaterial"
    assert c.divisions == 50
    assert c.thickness == 2.0


def test_circle_curve_custom_values():
    """CircleCurve accepts custom values."""
    c = CircleCurve(
        position=(1.0, 2.0, 3.0),
        radius=10.0,
        start_angle=25.0,
        end_angle=75.0,
        rotation=(0.1, 0.2, 0.3),
        scale=(2.0, 1.0, 1.0),
        color="#FF0000",
        material="LineDashedMaterial",
        divisions=100,
        thickness=5.0,
    )
    assert c.position == (1.0, 2.0, 3.0)
    assert c.radius == 10.0
    assert c.start_angle == 25.0
    assert c.end_angle == 75.0


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("radius", 0.0),
        ("radius", 101.0),
        ("start_angle", -1.0),
        ("start_angle", 101.0),
        ("end_angle", -1.0),
        ("end_angle", 101.0),
        ("divisions", 0),
        ("divisions", 201),
        ("thickness", 0.4),
        ("thickness", 10.1),
    ],
)
def test_circle_curve_validation_rejects_out_of_range(field: str, value: float):
    """Fields with ge/le constraints reject out-of-range values."""
    with pytest.raises(ValidationError):
        CircleCurve(**{field: value})


def test_get_interpolated_points_full_circle():
    """Full circle (0-100%) returns points forming a closed loop."""
    c = CircleCurve(radius=1.0, start_angle=0.0, end_angle=100.0, divisions=36)
    points = c.get_interpolated_points()
    assert isinstance(points, np.ndarray)
    assert points.shape == (37, 3)  # divisions + 1 points
    # All points should be at radius=1 from center in XY plane
    distances = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    np.testing.assert_allclose(distances, 1.0, atol=1e-10)
    # Z should be 0 (XY plane before rotation)
    np.testing.assert_allclose(points[:, 2], 0.0, atol=1e-10)
    # First and last point should be the same (full circle)
    np.testing.assert_allclose(points[0], points[-1], atol=1e-10)


def test_get_interpolated_points_half_circle():
    """Half circle (0-50%) returns a semicircle."""
    c = CircleCurve(radius=2.0, start_angle=0.0, end_angle=50.0, divisions=20)
    points = c.get_interpolated_points()
    assert points.shape == (21, 3)
    distances = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    np.testing.assert_allclose(distances, 2.0, atol=1e-10)


def test_get_interpolated_points_with_position():
    """Center offset translates all points."""
    c = CircleCurve(
        position=(10.0, 20.0, 30.0),
        radius=1.0,
        start_angle=0.0,
        end_angle=100.0,
        divisions=4,
    )
    points = c.get_interpolated_points()
    # All points should be offset by center
    center = np.array([10.0, 20.0, 30.0])
    distances_from_center = np.sqrt(((points - center) ** 2).sum(axis=1))
    np.testing.assert_allclose(distances_from_center, 1.0, atol=1e-10)


def test_get_interpolated_points_with_scale():
    """Scale creates an ellipse."""
    c = CircleCurve(
        radius=1.0,
        scale=(2.0, 1.0, 1.0),
        start_angle=0.0,
        end_angle=100.0,
        divisions=36,
    )
    points = c.get_interpolated_points()
    # X range should be [-2, 2], Y range should be [-1, 1]
    assert points[:, 0].max() == pytest.approx(2.0, abs=0.1)
    assert points[:, 1].max() == pytest.approx(1.0, abs=0.1)


def test_get_interpolated_points_with_rotation():
    """Rotation orients the circle plane."""
    # Rotate 90 degrees around X axis: XY plane -> XZ plane
    c = CircleCurve(
        radius=1.0,
        rotation=(math.pi / 2, 0.0, 0.0),
        start_angle=0.0,
        end_angle=100.0,
        divisions=36,
    )
    points = c.get_interpolated_points()
    # After rotating XY plane 90 deg around X: Y->Z, Z->-Y
    # So Y should be ~0 and Z should have extent
    np.testing.assert_allclose(points[:, 1], 0.0, atol=1e-10)
    assert abs(points[:, 2].max() - 1.0) < 0.1


def test_circle_curve_numpy_coercion():
    """Numpy arrays are coerced to tuples."""
    c = CircleCurve(position=np.array([1.0, 2.0, 3.0]))
    assert c.position == (1.0, 2.0, 3.0)


def test_circle_curve_in_geometry_registry():
    """CircleCurve is registered in the geometries dict."""
    from zndraw.geometries import geometries

    assert "CircleCurve" in geometries


def test_circle_curve_json_schema():
    """JSON schema includes expected format metadata."""
    schema = CircleCurve.model_json_schema()
    props = schema["properties"]
    assert props["radius"]["format"] == "range"
    assert props["start_angle"]["format"] == "range"
    assert props["end_angle"]["format"] == "range"
    assert props["color"]["format"] == "color"
    assert props["divisions"]["format"] == "range"
    assert props["thickness"]["format"] == "range"
