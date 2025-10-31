"""Tests for the Three.js materials system.

This test suite verifies that material Pydantic models work correctly, including
preset strings, full material objects, integration with geometries, and schema
customization for the UI.
"""

import pytest

from zndraw.geometries import Sphere
from zndraw.materials import (
    PRESETS,
    MeshBasicMaterial,
    MeshLambertMaterial,
    MeshPhongMaterial,
    MeshPhysicalMaterial,
    MeshStandardMaterial,
    MeshToonMaterial,
    get_all_preset_names,
    get_preset,
)

# --- Basic Material Creation Tests ---


def test_mesh_basic_material_creation():
    """Test creating a MeshBasicMaterial with default values."""
    material = MeshBasicMaterial()

    assert material.to_three_type() == "meshBasicMaterial"
    assert material.wireframe is False
    assert material.flatShading is False
    assert material.transparent is False


def test_mesh_standard_material_creation():
    """Test creating a MeshStandardMaterial with custom values."""
    material = MeshStandardMaterial(
        roughness=0.5, metalness=0.8, emissive="#FF0000", emissiveIntensity=0.5
    )

    assert material.to_three_type() == "meshStandardMaterial"
    assert material.roughness == 0.5
    assert material.metalness == 0.8
    assert material.emissive == "#FF0000"
    assert material.emissiveIntensity == 0.5


def test_mesh_physical_material_creation():
    """Test creating a MeshPhysicalMaterial with advanced properties."""
    material = MeshPhysicalMaterial(
        roughness=0.1,
        metalness=0.0,
        transmission=0.9,
        thickness=0.5,
        ior=1.5,
        clearcoat=1.0,
        clearcoatRoughness=0.1,
    )

    assert material.to_three_type() == "meshPhysicalMaterial"
    assert material.transmission == 0.9
    assert material.thickness == 0.5
    assert material.ior == 1.5
    assert material.clearcoat == 1.0
    assert material.clearcoatRoughness == 0.1


def test_mesh_phong_material_creation():
    """Test creating a MeshPhongMaterial with specular properties."""
    material = MeshPhongMaterial(
        emissive="#0000FF", emissiveIntensity=0.3, specular="#FFFFFF", shininess=50.0
    )

    assert material.to_three_type() == "meshPhongMaterial"
    assert material.emissive == "#0000FF"
    assert material.emissiveIntensity == 0.3
    assert material.specular == "#FFFFFF"
    assert material.shininess == 50.0


def test_mesh_lambert_material_creation():
    """Test creating a MeshLambertMaterial."""
    material = MeshLambertMaterial(emissive="#00FF00", emissiveIntensity=0.8)

    assert material.to_three_type() == "meshLambertMaterial"
    assert material.emissive == "#00FF00"
    assert material.emissiveIntensity == 0.8


def test_mesh_toon_material_creation():
    """Test creating a MeshToonMaterial."""
    material = MeshToonMaterial(wireframe=True)

    assert material.to_three_type() == "meshToonMaterial"
    assert material.wireframe is True


# --- Property Validation Tests ---


def test_roughness_validation():
    """Test that roughness is validated to be between 0 and 1."""
    # Valid values
    MeshStandardMaterial(roughness=0.0)
    MeshStandardMaterial(roughness=1.0)
    MeshStandardMaterial(roughness=0.5)

    # Invalid values should raise ValidationError
    with pytest.raises(Exception):  # Pydantic ValidationError
        MeshStandardMaterial(roughness=-0.1)

    with pytest.raises(Exception):
        MeshStandardMaterial(roughness=1.1)


def test_metalness_validation():
    """Test that metalness is validated to be between 0 and 1."""
    # Valid values
    MeshStandardMaterial(metalness=0.0)
    MeshStandardMaterial(metalness=1.0)
    MeshStandardMaterial(metalness=0.5)

    # Invalid values
    with pytest.raises(Exception):
        MeshStandardMaterial(metalness=-0.1)

    with pytest.raises(Exception):
        MeshStandardMaterial(metalness=1.5)


def test_transmission_validation():
    """Test that transmission is validated to be between 0 and 1."""
    # Valid values
    MeshPhysicalMaterial(transmission=0.0)
    MeshPhysicalMaterial(transmission=1.0)
    MeshPhysicalMaterial(transmission=0.5)

    # Invalid values
    with pytest.raises(Exception):
        MeshPhysicalMaterial(transmission=-0.1)

    with pytest.raises(Exception):
        MeshPhysicalMaterial(transmission=1.1)


def test_ior_validation():
    """Test that IOR is validated to be within valid range."""
    # Valid values
    MeshPhysicalMaterial(ior=1.0)
    MeshPhysicalMaterial(ior=1.5)
    MeshPhysicalMaterial(ior=2.333)

    # Invalid values
    with pytest.raises(Exception):
        MeshPhysicalMaterial(ior=0.5)

    with pytest.raises(Exception):
        MeshPhysicalMaterial(ior=3.0)


# --- Serialization Tests ---


def test_material_serialization():
    """Test that materials serialize to dict correctly."""
    material = MeshStandardMaterial(
        roughness=0.7,
        metalness=0.3,
        emissive="#FF0000",
        wireframe=True,
        flatShading=True,
    )

    data = material.model_dump()
    assert data["roughness"] == 0.7
    assert data["metalness"] == 0.3
    assert data["emissive"] == "#FF0000"
    assert data["wireframe"] is True
    assert data["flatShading"] is True


def test_material_json_serialization():
    """Test that materials serialize to JSON correctly."""
    material = MeshPhysicalMaterial(transmission=0.9, thickness=0.5)

    json_str = material.model_dump_json()
    assert "transmission" in json_str
    assert "0.9" in json_str
    assert "thickness" in json_str
    assert "0.5" in json_str


def test_material_deserialization():
    """Test that materials can be reconstructed from dict."""
    data = {
        "roughness": 0.5,
        "metalness": 0.8,
        "emissive": "#00FF00",
        "emissiveIntensity": 0.6,
        "wireframe": False,
        "flatShading": False,
        "transparent": False,
        "polygonOffset": False,
        "polygonOffsetFactor": 0,
    }

    material = MeshStandardMaterial(**data)
    assert material.roughness == 0.5
    assert material.metalness == 0.8
    assert material.emissive == "#00FF00"
    assert material.emissiveIntensity == 0.6


# --- Integration with Geometries Tests ---


def test_sphere_with_material_preset():
    """Test that Sphere accepts material preset strings."""
    sphere = Sphere(
        position=[(0, 0, 0)],
        color=["#FF0000"],
        material="MeshPhysicalMaterial_matt",
    )

    assert sphere.material == "MeshPhysicalMaterial_matt"
    assert isinstance(sphere.material, str)


def test_sphere_with_material_object():
    """Test that Sphere accepts full material objects."""
    material = MeshStandardMaterial(roughness=0.3, metalness=0.7, wireframe=True)
    sphere = Sphere(position=[(0, 0, 0)], color=["#FF0000"], material=material)

    assert isinstance(sphere.material, MeshStandardMaterial)
    assert sphere.material.roughness == 0.3
    assert sphere.material.metalness == 0.7
    assert sphere.material.wireframe is True


def test_sphere_with_physical_material():
    """Test that Sphere accepts MeshPhysicalMaterial objects."""
    material = MeshPhysicalMaterial(
        roughness=0.0, transmission=1.0, thickness=1.0, ior=1.5
    )
    sphere = Sphere(position=[(0, 0, 0)], color=["#FF0000"], material=material)

    assert isinstance(sphere.material, MeshPhysicalMaterial)
    assert sphere.material.transmission == 1.0
    assert sphere.material.ior == 1.5


def test_sphere_serialization_with_material_object():
    """Test that Sphere with material object serializes correctly."""
    material = MeshStandardMaterial(roughness=0.5, metalness=0.5)
    sphere = Sphere(position=[(0, 0, 0)], color=["#FF0000"], material=material)

    data = sphere.model_dump()
    assert "material" in data
    assert data["material"]["roughness"] == 0.5
    assert data["material"]["metalness"] == 0.5


def test_sphere_deserialization_with_material_object():
    """Test that Sphere can be reconstructed from dict with material object."""
    data = {
        "position": [(0, 0, 0)],
        "color": ["#FF0000"],
        "material": {
            "material_type": "MeshStandardMaterial",
            "roughness": 0.7,
            "metalness": 0.3,
            "emissive": "#000000",
            "emissiveIntensity": 1.0,
            "wireframe": False,
            "flatShading": False,
            "transparent": False,
            "polygonOffset": False,
            "polygonOffsetFactor": 0,
        },
    }

    sphere = Sphere(**data)
    assert isinstance(sphere.material, MeshStandardMaterial)
    assert sphere.material.roughness == 0.7
    assert sphere.material.metalness == 0.3


def test_sphere_roundtrip_with_material():
    """Test full serialization roundtrip of Sphere with material."""
    original_sphere = Sphere(
        position=[(1, 2, 3)],
        color=["#FF0000"],
        material=MeshPhysicalMaterial(
            roughness=0.2, transmission=0.8, thickness=0.5, clearcoat=1.0
        ),
    )

    # Serialize to dict
    data = original_sphere.model_dump()

    # Deserialize back to Sphere
    restored_sphere = Sphere(**data)

    # Verify material is preserved
    assert isinstance(restored_sphere.material, MeshPhysicalMaterial)
    assert restored_sphere.material.roughness == 0.2
    assert restored_sphere.material.transmission == 0.8
    assert restored_sphere.material.thickness == 0.5
    assert restored_sphere.material.clearcoat == 1.0


# --- Schema Customization Tests ---


def test_sphere_schema_has_material_custom_type():
    """Test that Sphere schema has x-custom-type for material field."""
    schema = Sphere.model_json_schema()

    assert "material" in schema["properties"]
    assert "x-custom-type" in schema["properties"]["material"]
    assert schema["properties"]["material"]["x-custom-type"] == "three-material"


def test_material_schema_includes_all_types():
    """Test that material field schema includes all material types."""
    schema = Sphere.model_json_schema()
    material_schema = schema["properties"]["material"]

    # Should have anyOf with preset strings and material object types
    assert "anyOf" in material_schema


# --- Preset Tests ---


def test_get_all_preset_names():
    """Test that get_all_preset_names returns all preset names."""
    presets = get_all_preset_names()

    assert isinstance(presets, list)
    assert len(presets) == 11
    assert "MeshPhysicalMaterial_matt" in presets
    assert "MeshPhysicalMaterial_glass" in presets
    assert "MeshStandardMaterial_metallic" in presets
    assert "MeshBasicMaterial" in presets


def test_get_preset():
    """Test that get_preset returns correct material instances."""
    matt = get_preset("MeshPhysicalMaterial_matt")
    assert isinstance(matt, MeshPhysicalMaterial)
    assert matt.roughness == 0.9
    assert matt.reflectivity == 0.1

    glass = get_preset("MeshPhysicalMaterial_glass")
    assert isinstance(glass, MeshPhysicalMaterial)
    assert glass.transmission == 1.0
    assert glass.transparent is True

    metallic = get_preset("MeshStandardMaterial_metallic")
    assert isinstance(metallic, MeshStandardMaterial)
    assert metallic.metalness == 1.0

    basic = get_preset("MeshBasicMaterial")
    assert isinstance(basic, MeshBasicMaterial)
    assert basic.toneMapped is False


def test_get_preset_nonexistent():
    """Test that get_preset returns None for nonexistent presets."""
    result = get_preset("NonExistentPreset")
    assert result is None


def test_presets_dict_structure():
    """Test that PRESETS dict has correct structure."""
    assert isinstance(PRESETS, dict)
    assert len(PRESETS) == 11

    # All keys are strings
    for key in PRESETS.keys():
        assert isinstance(key, str)

    # All values are MaterialType instances
    for value in PRESETS.values():
        assert isinstance(
            value,
            (
                MeshBasicMaterial,
                MeshStandardMaterial,
                MeshPhysicalMaterial,
                MeshToonMaterial,
                MeshLambertMaterial,
                MeshPhongMaterial,
            ),
        )


# --- Parametrized Tests ---


@pytest.mark.parametrize(
    "material_class,expected_type",
    [
        (MeshBasicMaterial, "meshBasicMaterial"),
        (MeshStandardMaterial, "meshStandardMaterial"),
        (MeshPhysicalMaterial, "meshPhysicalMaterial"),
        (MeshToonMaterial, "meshToonMaterial"),
        (MeshLambertMaterial, "meshLambertMaterial"),
        (MeshPhongMaterial, "meshPhongMaterial"),
    ],
)
def test_material_types(material_class, expected_type):
    """Test that all material types return correct Three.js type."""
    material = material_class()
    assert material.to_three_type() == expected_type


@pytest.mark.parametrize(
    "property_name,valid_values,invalid_values",
    [
        ("roughness", [0.0, 0.5, 1.0], [-0.1, 1.1]),
        ("metalness", [0.0, 0.5, 1.0], [-0.1, 1.5]),
        ("transmission", [0.0, 0.5, 1.0], [-0.1, 1.2]),
        ("clearcoat", [0.0, 0.5, 1.0], [-0.1, 1.3]),
    ],
)
def test_property_range_validation(property_name, valid_values, invalid_values):
    """Test that material properties validate ranges correctly."""
    # Test valid values
    for value in valid_values:
        if property_name in ["transmission", "clearcoat"]:
            MeshPhysicalMaterial(**{property_name: value})
        else:
            MeshStandardMaterial(**{property_name: value})

    # Test invalid values
    for value in invalid_values:
        with pytest.raises(Exception):
            if property_name in ["transmission", "clearcoat"]:
                MeshPhysicalMaterial(**{property_name: value})
            else:
                MeshStandardMaterial(**{property_name: value})


@pytest.mark.parametrize("wireframe", [True, False])
@pytest.mark.parametrize("flatShading", [True, False])
def test_common_properties(wireframe, flatShading):
    """Test that common properties work across all material types."""
    material = MeshStandardMaterial(wireframe=wireframe, flatShading=flatShading)
    assert material.wireframe == wireframe
    assert material.flatShading == flatShading


# --- Edge Cases ---


def test_mesh_physical_inherits_standard_properties():
    """Test that MeshPhysicalMaterial inherits MeshStandardMaterial properties."""
    material = MeshPhysicalMaterial(
        roughness=0.5,  # From MeshStandardMaterial
        metalness=0.3,  # From MeshStandardMaterial
        transmission=0.9,  # MeshPhysicalMaterial-specific
    )

    assert material.roughness == 0.5
    assert material.metalness == 0.3
    assert material.transmission == 0.9


def test_default_values():
    """Test that all materials have sensible default values."""
    basic = MeshBasicMaterial()
    assert basic.wireframe is False
    assert basic.transparent is False

    standard = MeshStandardMaterial()
    assert standard.roughness == 1.0
    assert standard.metalness == 0.0
    assert standard.emissive == "#000000"

    physical = MeshPhysicalMaterial()
    assert physical.transmission == 0.0
    assert physical.ior == 1.5
    assert physical.clearcoat == 0.0


def test_polygon_offset_properties():
    """Test polygon offset properties on materials."""
    material = MeshStandardMaterial(polygonOffset=True, polygonOffsetFactor=1.0)

    assert material.polygonOffset is True
    assert material.polygonOffsetFactor == 1.0
