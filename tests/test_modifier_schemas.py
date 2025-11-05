"""Test modifier JSON schemas for dynamic enum curve fields."""

import pytest

from zndraw.extensions import modifiers


@pytest.mark.parametrize(
    "modifier_class",
    [
        modifiers.Connect,
        modifiers.Rotate,
        modifiers.Translate,
        modifiers.AddLineParticles,
    ],
)
def test_modifier_curve_field_schema(modifier_class):
    """Test that modifiers with curve fields have correct JSON schema annotations."""
    schema = modifier_class.model_json_schema()

    # Verify curve field exists
    assert "curve" in schema["properties"], f"{modifier_class.__name__} missing 'curve' field"

    curve_schema = schema["properties"]["curve"]

    # Verify x-custom-type is set to dynamic-enum
    assert (
        curve_schema.get("x-custom-type") == "dynamic-enum"
    ), f"{modifier_class.__name__} curve field missing x-custom-type: dynamic-enum"

    # Verify x-features contains dynamic-geometries
    assert "x-features" in curve_schema, f"{modifier_class.__name__} curve field missing x-features"
    assert (
        "dynamic-geometries" in curve_schema["x-features"]
    ), f"{modifier_class.__name__} curve field x-features missing 'dynamic-geometries'"

    # Verify x-geometry-filter is set to Curve
    assert (
        curve_schema.get("x-geometry-filter") == "Curve"
    ), f"{modifier_class.__name__} curve field missing x-geometry-filter: Curve"


def test_connect_schema():
    """Test Connect modifier schema structure."""
    schema = modifiers.Connect.model_json_schema()

    assert "properties" in schema
    assert "curve" in schema["properties"]
    assert schema["properties"]["curve"]["type"] == "string"
    assert schema["properties"]["curve"]["description"] == "Name of the curve geometry to use"


def test_rotate_schema():
    """Test Rotate modifier schema structure."""
    schema = modifiers.Rotate.model_json_schema()

    assert "properties" in schema
    assert "curve" in schema["properties"]
    assert "angle" in schema["properties"]
    assert "direction" in schema["properties"]
    assert "steps" in schema["properties"]


def test_translate_schema():
    """Test Translate modifier schema structure."""
    schema = modifiers.Translate.model_json_schema()

    assert "properties" in schema
    assert "curve" in schema["properties"]
    assert "steps" in schema["properties"]


def test_add_line_particles_schema():
    """Test AddLineParticles modifier schema structure."""
    schema = modifiers.AddLineParticles.model_json_schema()

    assert "properties" in schema
    assert "curve" in schema["properties"]
    assert "symbol" in schema["properties"]
