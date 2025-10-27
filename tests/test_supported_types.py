"""Tests for supported file types and upload endpoints.

This module tests that the supported-types and upload endpoints are always available
regardless of whether the file browser feature is enabled, since they're needed for
drag/drop upload which is a core feature.
"""

import io

import requests


def test_supported_types_without_file_browser(server):
    """Test that supported-types endpoint works when file browser is disabled."""
    # The server fixture doesn't enable file browser by default
    # This test verifies that we can still get supported types
    response = requests.get(f"{server}/api/file-browser/supported-types")

    assert response.status_code == 200

    data = response.json()

    # Verify response structure
    assert "extensions" in data
    assert "descriptions" in data
    assert "backends" in data

    # Verify extensions is a list
    assert isinstance(data["extensions"], list)
    assert len(data["extensions"]) > 0

    # Verify common extensions are present
    assert ".xyz" in data["extensions"]
    assert ".pdb" in data["extensions"]
    assert ".cif" in data["extensions"]
    assert ".h5" in data["extensions"]

    # Verify descriptions is a dict
    assert isinstance(data["descriptions"], dict)
    assert ".xyz" in data["descriptions"]

    # Verify backends is a dict
    assert isinstance(data["backends"], dict)
    assert "xyz" in data["backends"]
    assert "ASE" in data["backends"]["xyz"]


def test_supported_types_includes_all_common_formats(server):
    """Test that supported-types includes all commonly used file formats."""
    response = requests.get(f"{server}/api/file-browser/supported-types")

    assert response.status_code == 200

    data = response.json()
    extensions = data["extensions"]

    # Common molecular structure formats
    common_formats = [
        ".xyz",  # XYZ format
        ".extxyz",  # Extended XYZ
        ".pdb",  # Protein Data Bank
        ".cif",  # Crystallographic Information File
        ".h5",  # HDF5/H5MD
        ".h5md",  # H5MD
        ".gro",  # GROMACS
        ".mol",  # MDL Molfile
        ".sdf",  # Structure Data File
        ".db",  # ASE database
        ".json",  # JSON database
    ]

    for fmt in common_formats:
        assert fmt in extensions, f"Format {fmt} should be supported"


def test_supported_types_descriptions_format(server):
    """Test that descriptions follow expected format."""
    response = requests.get(f"{server}/api/file-browser/supported-types")

    assert response.status_code == 200

    data = response.json()
    descriptions = data["descriptions"]

    # Each extension should have a description
    for ext in data["extensions"]:
        assert ext in descriptions
        description = descriptions[ext]
        assert isinstance(description, str)
        assert len(description) > 0
        # Description should mention the backend
        assert "Supported by" in description


def test_upload_without_file_browser(server):
    """Test that upload endpoint works when file browser is disabled."""
    # Create a simple XYZ file content
    xyz_content = """2
Test molecule
H 0.0 0.0 0.0
H 1.0 0.0 0.0
"""

    # Prepare file upload
    files = {
        "file": ("test.xyz", io.BytesIO(xyz_content.encode("utf-8")), "text/plain")
    }

    # Try to upload the file
    response = requests.post(f"{server}/api/file-browser/upload", files=files)

    # Should successfully accept the upload (returns 200 or 202)
    assert response.status_code in [200, 202]

    data = response.json()

    # Should return success status
    assert "status" in data
    assert data["status"] in ["success", "processing", "queued"]

    # Should return room information
    assert "room" in data


def test_upload_sets_drag_drop_description(server, s22):
    """Test that drag/drop uploads include source information in room description."""
    import time

    from zndraw import ZnDraw

    # Create a simple XYZ file with a molecule
    xyz_content = """2
Test molecule
H 0.0 0.0 0.0
H 1.0 0.0 0.0
"""

    # Prepare file upload with a specific filename
    filename = "my_molecule.xyz"
    files = {"file": (filename, io.BytesIO(xyz_content.encode("utf-8")), "text/plain")}

    # Upload the file
    response = requests.post(f"{server}/api/file-browser/upload", files=files)

    assert response.status_code in [200, 202]
    data = response.json()
    room_id = data["room"]

    # Wait a bit for the async task to complete
    time.sleep(2)

    # Get room details to check description
    room_response = requests.get(f"{server}/api/rooms/{room_id}")

    if room_response.status_code == 200:
        room_data = room_response.json()
        description = room_data.get("description", "")

        # Description should mention drag & drop
        assert (
            "drag & drop" in description.lower()
            or "drag and drop" in description.lower()
        )
        # Description should include the filename
        assert filename in description
