"""Unit tests for Camera geometry model.

Tests Camera Pydantic model validation, defaults, and serialization.
"""

import pytest
from pydantic import ValidationError

from zndraw.geometries import Camera
from zndraw.geometries.camera import CameraType
from zndraw.transformations import CurveAttachment


def test_camera_defaults():
    """Camera() should have correct default values."""
    cam = Camera()

    assert cam.position == (-10.0, 10.0, 30.0)
    assert cam.target == (0.0, 0.0, 0.0)
    assert cam.up == (0.0, 1.0, 0.0)
    assert cam.camera_type == CameraType.PERSPECTIVE
    assert cam.fov == 50.0
    assert cam.near == 0.1
    assert cam.far == 1000.0
    assert cam.zoom == 1.0
    assert cam.helper_visible is False
    assert cam.helper_color == "#00ff00"
    assert cam.show_crosshair is False
    assert cam.preserve_drawing_buffer is False
    assert cam.active is True


def test_camera_with_custom_values():
    """Camera with custom values should store them correctly."""
    cam = Camera(
        position=(10.0, 20.0, 30.0),
        target=(1.0, 2.0, 3.0),
        up=(0.0, 0.0, 1.0),
        camera_type=CameraType.ORTHOGRAPHIC,
        fov=60.0,
        near=0.5,
        far=500.0,
        zoom=2.0,
        helper_visible=False,
        helper_color="#ff0000",
    )

    assert cam.position == (10.0, 20.0, 30.0)
    assert cam.target == (1.0, 2.0, 3.0)
    assert cam.up == (0.0, 0.0, 1.0)
    assert cam.camera_type == CameraType.ORTHOGRAPHIC
    assert cam.fov == 60.0
    assert cam.near == 0.5
    assert cam.far == 500.0
    assert cam.zoom == 2.0
    assert cam.helper_visible is False
    assert cam.helper_color == "#ff0000"


def test_camera_validation_far_greater_than_near():
    """Camera should reject far <= near."""
    with pytest.raises(ValidationError) as exc_info:
        Camera(near=10.0, far=5.0)

    assert "far" in str(exc_info.value).lower()
    assert "near" in str(exc_info.value).lower()


def test_camera_validation_far_equals_near():
    """Camera should reject far == near."""
    with pytest.raises(ValidationError) as exc_info:
        Camera(near=10.0, far=10.0)

    assert "far" in str(exc_info.value).lower()


def test_camera_validation_zero_up_vector():
    """Camera should reject zero up vector."""
    with pytest.raises(ValidationError) as exc_info:
        Camera(up=(0.0, 0.0, 0.0))

    assert "up" in str(exc_info.value).lower() or "zero" in str(exc_info.value).lower()


def test_camera_validation_near_positive():
    """Camera near must be positive."""
    with pytest.raises(ValidationError):
        Camera(near=-1.0)


def test_camera_validation_far_positive():
    """Camera far must be positive."""
    with pytest.raises(ValidationError):
        Camera(far=-1.0)


def test_camera_validation_fov_range():
    """Camera fov must be between 1 and 179 (inclusive)."""
    with pytest.raises(ValidationError):
        Camera(fov=0.0)

    with pytest.raises(ValidationError):
        Camera(fov=180.0)

    # Valid fov values
    cam1 = Camera(fov=1.0)
    assert cam1.fov == 1.0

    cam2 = Camera(fov=179.0)
    assert cam2.fov == 179.0


def test_camera_validation_zoom_positive():
    """Camera zoom must be positive."""
    with pytest.raises(ValidationError):
        Camera(zoom=0.0)

    with pytest.raises(ValidationError):
        Camera(zoom=-1.0)


def test_camera_with_curve_attachment_position():
    """Camera can use CurveAttachment for position."""
    attachment = CurveAttachment(geometry_key="my_curve", progress=0.5)
    cam = Camera(position=attachment, target=(0.0, 0.0, 0.0))

    assert isinstance(cam.position, CurveAttachment)
    assert cam.position.geometry_key == "my_curve"
    assert cam.position.progress == 0.5


def test_camera_with_curve_attachment_target():
    """Camera can use CurveAttachment for target."""
    attachment = CurveAttachment(geometry_key="target_curve", progress=0.8)
    cam = Camera(position=(10.0, 10.0, 10.0), target=attachment)

    assert isinstance(cam.target, CurveAttachment)
    assert cam.target.geometry_key == "target_curve"
    assert cam.target.progress == 0.8


def test_camera_with_both_curve_attachments():
    """Camera can use CurveAttachment for both position and target."""
    pos_attachment = CurveAttachment(geometry_key="pos_curve", progress=0.3)
    target_attachment = CurveAttachment(geometry_key="target_curve", progress=0.7)
    cam = Camera(position=pos_attachment, target=target_attachment)

    assert isinstance(cam.position, CurveAttachment)
    assert isinstance(cam.target, CurveAttachment)


def test_camera_type_enum():
    """CameraType enum values should match Three.js names."""
    assert CameraType.PERSPECTIVE.value == "PerspectiveCamera"
    assert CameraType.ORTHOGRAPHIC.value == "OrthographicCamera"


def test_camera_serialization():
    """Camera should serialize to JSON correctly."""
    cam = Camera(
        position=(1.0, 2.0, 3.0),
        target=(4.0, 5.0, 6.0),
        fov=60.0,
    )

    data = cam.model_dump()

    assert data["position"] == (1.0, 2.0, 3.0)
    assert data["target"] == (4.0, 5.0, 6.0)
    assert data["fov"] == 60.0
    assert data["camera_type"] == "PerspectiveCamera"


def test_camera_serialization_with_curve_attachment():
    """Camera with CurveAttachment should serialize correctly."""
    attachment = CurveAttachment(geometry_key="curve1", progress=0.5)
    cam = Camera(position=attachment)

    data = cam.model_dump()

    assert isinstance(data["position"], dict)
    assert data["position"]["geometry_key"] == "curve1"
    assert data["position"]["progress"] == 0.5


def test_camera_deserialization():
    """Camera should deserialize from dict correctly."""
    data = {
        "position": (10.0, 20.0, 30.0),
        "target": (0.0, 0.0, 0.0),
        "fov": 90.0,
        "near": 0.5,
        "far": 2000.0,
    }

    cam = Camera.model_validate(data)

    assert cam.position == (10.0, 20.0, 30.0)
    assert cam.fov == 90.0
    assert cam.near == 0.5
    assert cam.far == 2000.0


def test_camera_frozen():
    """Camera model should be frozen (immutable)."""
    cam = Camera()

    with pytest.raises(ValidationError):
        cam.fov = 60.0


def test_camera_model_copy():
    """Camera can be copied with modifications using model_copy."""
    cam = Camera(fov=60.0)
    cam2 = cam.model_copy(update={"fov": 90.0})

    assert cam.fov == 60.0
    assert cam2.fov == 90.0


def test_camera_json_schema():
    """Camera should generate valid JSON schema."""
    schema = Camera.model_json_schema()

    assert "properties" in schema
    assert "position" in schema["properties"]
    assert "target" in schema["properties"]
    assert "fov" in schema["properties"]
    assert "camera_type" in schema["properties"]


def test_camera_helper_color_schema():
    """Camera helper_color should have color-picker hint in schema."""
    schema = Camera.model_json_schema()

    helper_color_schema = schema["properties"]["helper_color"]
    assert helper_color_schema.get("x-custom-type") == "color-picker"


def test_camera_deserialization_with_curve_attachment_dict():
    """Camera should deserialize CurveAttachment from dict (JSON API format).

    This is the critical test for union type handling - when data comes from
    JSON (frontend or API), CurveAttachment is a dict, not an object.
    Pydantic must correctly convert it to CurveAttachment.
    """
    data = {
        "position": {
            "type": "curve_attachment",
            "geometry_key": "flight_path",
            "progress": 0.5,
        },
        "target": (0.0, 0.0, 0.0),
        "fov": 60.0,
    }

    cam = Camera.model_validate(data)

    # Verify position was converted to CurveAttachment (not left as dict)
    assert isinstance(cam.position, CurveAttachment)
    assert cam.position.geometry_key == "flight_path"
    assert cam.position.progress == 0.5

    # Verify target remained as tuple
    assert isinstance(cam.target, tuple)
    assert cam.target == (0.0, 0.0, 0.0)


def test_camera_deserialization_both_curve_attachments_as_dicts():
    """Camera should handle both position and target as CurveAttachment dicts."""
    data = {
        "position": {
            "type": "curve_attachment",
            "geometry_key": "pos_curve",
            "progress": 0.3,
        },
        "target": {
            "type": "curve_attachment",
            "geometry_key": "target_curve",
            "progress": 0.7,
        },
    }

    cam = Camera.model_validate(data)

    assert isinstance(cam.position, CurveAttachment)
    assert cam.position.geometry_key == "pos_curve"
    assert isinstance(cam.target, CurveAttachment)
    assert cam.target.geometry_key == "target_curve"
