"""Tests for the transformation system.

This test suite verifies that transformations work correctly for filtering
geometry data based on frame data. Focuses on InArrayTransform which filters
arrays by indices extracted from frame data (e.g., constraint indices).
"""

import ase
import ase.constraints
import pytest

from asebytes import encode
from zndraw.geometries import Sphere
from zndraw.transformations import InArrayTransform


def _decode_msgpack_dict(obj):
    """Recursively decode msgpack dict[bytes, bytes] to regular Python dict."""
    import msgpack
    import msgpack_numpy as m

    if isinstance(obj, dict):
        result = {}
        for k, v in obj.items():
            # Decode key
            key = k.decode() if isinstance(k, bytes) else k
            # Decode value if it's bytes (msgpack-encoded)
            if isinstance(v, bytes):
                try:
                    value = msgpack.unpackb(v, object_hook=m.decode, strict_map_key=False)
                    value = _decode_msgpack_dict(value)
                except:
                    value = v
            else:
                value = _decode_msgpack_dict(v)
            result[key] = value
        return result
    elif isinstance(obj, list):
        return [_decode_msgpack_dict(item) for item in obj]
    else:
        return obj


def test_in_array_transform_creation():
    """Test creating an InArrayTransform with required fields."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )

    assert transform.type == "in_array"
    assert transform.source == "constraints"
    assert transform.path == "0.kwargs.indices"
    assert transform.filter == "arrays.positions"


def test_in_array_transform_serialization():
    """Test that InArrayTransform serializes to JSON correctly."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )

    # Test Pydantic model_dump
    data = transform.model_dump()
    assert data["type"] == "in_array"
    assert data["source"] == "constraints"
    assert data["path"] == "0.kwargs.indices"
    assert data["filter"] == "arrays.positions"

    # Test Pydantic model_dump_json
    json_str = transform.model_dump_json()
    assert "in_array" in json_str
    assert "constraints" in json_str
    assert "arrays.positions" in json_str


def test_in_array_transform_deserialization():
    """Test that InArrayTransform can be reconstructed from dict."""
    data = {
        "type": "in_array",
        "source": "constraints",
        "path": "0.kwargs.indices",
        "filter": "arrays.positions",
    }

    transform = InArrayTransform(**data)
    assert transform.type == "in_array"
    assert transform.source == "constraints"
    assert transform.path == "0.kwargs.indices"
    assert transform.filter == "arrays.positions"


def test_sphere_with_transform_position():
    """Test that Sphere geometry accepts InArrayTransform for position."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )

    sphere = Sphere(position=transform, radius=1.0, color=["#FF0000"])

    # The position field should contain the transform
    assert isinstance(sphere.position, InArrayTransform)
    assert sphere.position.source == "constraints"
    assert sphere.position.filter == "arrays.positions"


def test_sphere_with_transform_radius():
    """Test that Sphere geometry accepts InArrayTransform for radius."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.radii"
    )

    sphere = Sphere(position=[(0, 0, 0)], radius=transform, color=["#FF0000"])

    # The radius field should contain the transform
    assert isinstance(sphere.radius, InArrayTransform)
    assert sphere.radius.source == "constraints"
    assert sphere.radius.filter == "arrays.radii"


def test_sphere_with_transform_color():
    """Test that Sphere geometry accepts InArrayTransform for color."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.colors"
    )

    sphere = Sphere(position=[(0, 0, 0)], radius=1.0, color=transform)

    # The color field should contain the transform
    assert isinstance(sphere.color, InArrayTransform)
    assert sphere.color.source == "constraints"
    assert sphere.color.filter == "arrays.colors"


def test_sphere_with_multiple_transforms():
    """Test that Sphere geometry can have transforms on multiple properties."""
    pos_transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )
    radius_transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.radii"
    )

    sphere = Sphere(position=pos_transform, radius=radius_transform, color=["#FF0000"])

    assert isinstance(sphere.position, InArrayTransform)
    assert isinstance(sphere.radius, InArrayTransform)
    assert sphere.position.filter == "arrays.positions"
    assert sphere.radius.filter == "arrays.radii"


def test_sphere_serialization_with_transform():
    """Test that Sphere with transforms serializes correctly."""
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )

    sphere = Sphere(position=transform, radius=1.0, color=["#FF0000"])

    # Serialize to dict
    data = sphere.model_dump()
    assert data["position"]["type"] == "in_array"
    assert data["position"]["source"] == "constraints"
    assert data["position"]["filter"] == "arrays.positions"

    # Serialize to JSON
    json_str = sphere.model_dump_json()
    assert "in_array" in json_str
    assert "constraints" in json_str


def test_sphere_deserialization_with_transform():
    """Test that Sphere can be reconstructed from dict with transforms."""
    data = {
        "position": {
            "type": "in_array",
            "source": "constraints",
            "path": "0.kwargs.indices",
            "filter": "arrays.positions",
        },
        "radius": 1.0,
        "color": ["#FF0000"],
    }

    sphere = Sphere(**data)
    assert isinstance(sphere.position, InArrayTransform)
    assert sphere.position.source == "constraints"
    assert sphere.position.filter == "arrays.positions"


def test_transform_with_fixatoms_constraint():
    """Integration test: Create atoms with FixAtoms constraint and verify transform path."""
    import msgpack
    import msgpack_numpy as m

    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraint = ase.constraints.FixAtoms(indices=[0, 2])
    atoms.set_constraint(constraint)

    # Serialize atoms to bytes and unpack for testing
    frame_bytes = encode(atoms)
    # Convert to regular Python dict for testing (simulates frontend decoding)
    frame_data = _decode_msgpack_dict(frame_bytes)

    # Verify constraint structure matches transform path
    assert "constraints" in frame_data
    assert len(frame_data["constraints"]) == 1
    assert frame_data["constraints"][0]["name"] == "FixAtoms"
    assert frame_data["constraints"][0]["kwargs"]["indices"] == [0, 2]

    # Create transform that would extract these indices
    transform = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )

    # Verify the path matches the data structure
    # In frontend: frame_data[transform.source][0]["kwargs"]["indices"]
    constraint_data = frame_data[transform.source][0]
    indices = constraint_data["kwargs"]["indices"]
    assert indices == [0, 2]


def test_transform_with_multiple_constraints():
    """Test transform path with multiple constraints."""
    import msgpack
    import msgpack_numpy as m

    atoms = ase.Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    constraints = [
        ase.constraints.FixAtoms(indices=[0]),
        ase.constraints.FixBondLength(1, 2),
    ]
    atoms.set_constraint(constraints)

    # Serialize and unpack for testing
    frame_bytes = encode(atoms)
    # Convert to regular Python dict for testing (simulates frontend decoding)
    frame_data = _decode_msgpack_dict(frame_bytes)

    # First constraint (FixAtoms)
    transform_0 = InArrayTransform(
        source="constraints", path="0.kwargs.indices", filter="arrays.positions"
    )
    assert frame_data["constraints"][0]["kwargs"]["indices"] == [0]

    # Second constraint (FixBondLengths)
    # Note: FixBondLength becomes FixBondLengths after serialization
    assert frame_data["constraints"][1]["name"] == "FixBondLengths"


def test_transform_type_literal():
    """Test that transform type is enforced as literal 'in_array'."""
    # Should work with correct type
    transform = InArrayTransform(
        type="in_array",
        source="constraints",
        path="0.kwargs.indices",
        filter="arrays.positions",
    )
    assert transform.type == "in_array"

    # Should fail with wrong type (Pydantic validation)
    with pytest.raises(Exception):  # ValidationError
        InArrayTransform(
            type="wrong_type",
            source="constraints",
            path="0.kwargs.indices",
            filter="arrays.positions",
        )


def test_transform_with_nested_path():
    """Test transform with deeply nested path."""
    # Test various path formats
    paths = [
        "0.kwargs.indices",  # Array index + object keys
        "FixAtoms.kwargs.indices",  # Object keys only
        "0.indices",  # Shorter path
    ]

    for path in paths:
        transform = InArrayTransform(
            source="constraints", path=path, filter="arrays.positions"
        )
        assert transform.path == path


def test_transform_roundtrip_in_sphere():
    """Test full serialization roundtrip of Sphere with transforms."""
    original_sphere = Sphere(
        position=InArrayTransform(
            source="constraints", path="0.kwargs.indices", filter="arrays.positions"
        ),
        radius=InArrayTransform(
            source="constraints", path="0.kwargs.indices", filter="arrays.radii"
        ),
        color=["#FF0000"],
    )

    # Serialize to dict
    data = original_sphere.model_dump()

    # Deserialize back to Sphere
    restored_sphere = Sphere(**data)

    # Verify transforms are preserved
    assert isinstance(restored_sphere.position, InArrayTransform)
    assert isinstance(restored_sphere.radius, InArrayTransform)
    assert restored_sphere.position.source == "constraints"
    assert restored_sphere.position.filter == "arrays.positions"
    assert restored_sphere.radius.filter == "arrays.radii"
