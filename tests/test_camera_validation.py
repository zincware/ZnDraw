"""Tests for Camera model validation."""

import pytest
from pydantic import ValidationError

from zndraw.geometries.camera import Camera, CameraType


def test_defaults_are_valid():
    """Default Camera() creates a valid instance."""
    cam = Camera()
    assert cam.fov == 50.0
    assert cam.near == 0.1
    assert cam.far == 1000.0
    assert cam.zoom == 1.0
    assert cam.camera_type == CameraType.PERSPECTIVE


def test_far_greater_than_near():
    """far <= near raises ValidationError."""
    with pytest.raises(ValidationError, match="far.*must be greater than near"):
        Camera(near=10.0, far=5.0)


def test_far_equal_to_near():
    """far == near raises ValidationError."""
    with pytest.raises(ValidationError, match="far.*must be greater than near"):
        Camera(near=10.0, far=10.0)


@pytest.mark.parametrize("fov", [0, -1, 180, 200])
def test_fov_out_of_range(fov: float):
    """FOV outside 1-179 raises ValidationError."""
    with pytest.raises(ValidationError):
        Camera(fov=fov)


@pytest.mark.parametrize("fov", [1, 50, 179])
def test_fov_valid_range(fov: float):
    """FOV within 1-179 is accepted."""
    cam = Camera(fov=fov)
    assert cam.fov == fov


@pytest.mark.parametrize("zoom", [0.0, -1.0])
def test_zoom_invalid(zoom: float):
    """Zero or negative zoom raises ValidationError."""
    with pytest.raises(ValidationError):
        Camera(zoom=zoom)


def test_zoom_valid():
    """Positive zoom within range is accepted."""
    cam = Camera(zoom=2.5)
    assert cam.zoom == 2.5


def test_frozen_model():
    """Camera is frozen — assignment raises."""
    cam = Camera()
    with pytest.raises(ValidationError):
        cam.fov = 90.0  # type: ignore[misc]


def test_serialization_roundtrip():
    """model_dump → model_validate preserves all fields."""
    cam = Camera(
        position=(1.0, 2.0, 3.0),
        target=(4.0, 5.0, 6.0),
        fov=75.0,
        near=0.5,
        far=500.0,
        zoom=2.0,
        camera_type=CameraType.ORTHOGRAPHIC,
    )
    data = cam.model_dump()
    restored = Camera.model_validate(data)
    assert restored == cam


def test_up_vector_zero_raises():
    """Zero up vector raises ValidationError."""
    with pytest.raises(ValidationError, match="up vector cannot be zero"):
        Camera(up=(0.0, 0.0, 0.0))


def test_near_below_minimum():
    """near < 0.01 raises ValidationError."""
    with pytest.raises(ValidationError):
        Camera(near=0.001)


def test_orthographic_camera():
    """Orthographic camera type is valid."""
    cam = Camera(camera_type=CameraType.ORTHOGRAPHIC)
    assert cam.camera_type == CameraType.ORTHOGRAPHIC
