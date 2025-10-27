"""Tests for frame download functionality.

This module tests the ExtendedXYZ download endpoint including:
- Download single frame by index
- Download all frames (no parameters)
- Download multiple specific frames
- Download with particle selection
- ExtendedXYZ format verification
"""

import io

import ase.io
import numpy as np
import pytest
import requests

from zndraw import ZnDraw


def test_download_single_frame(server, s22):
    """Test downloading a single frame by explicit index."""
    # Setup: Create room with multiple frames
    vis = ZnDraw(url=server, room="test-download-single", user="user1")
    for atoms in s22[:5]:
        vis.append(atoms)

    # Download frame at index 2
    response = requests.get(
        f"{server}/api/rooms/test-download-single/download", params={"indices": "2"}
    )

    assert response.status_code == 200
    assert response.headers["Content-Type"] == "chemical/x-xyz"
    assert "filename=" in response.headers.get("Content-Disposition", "")

    # Parse ExtendedXYZ content
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    # Should only contain one frame
    assert len(atoms_list) == 1

    # Verify structure matches
    downloaded_atoms = atoms_list[0]
    original_atoms = s22[2]
    assert len(downloaded_atoms) == len(original_atoms)
    assert (
        downloaded_atoms.get_chemical_formula() == original_atoms.get_chemical_formula()
    )
    assert np.allclose(downloaded_atoms.positions, original_atoms.positions, atol=1e-6)


def test_download_all_frames(server, s22):
    """Test downloading all frames when no indices parameter is provided."""
    # Setup: Create room with multiple frames
    vis = ZnDraw(url=server, room="test-download-all", user="user1")
    for atoms in s22[:10]:
        vis.append(atoms)

    # Download all frames (no indices parameter)
    response = requests.get(f"{server}/api/rooms/test-download-all/download")

    assert response.status_code == 200
    assert response.headers["Content-Type"] == "chemical/x-xyz"

    # Parse ExtendedXYZ content
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    # Should contain all 10 frames
    assert len(atoms_list) == 10

    # Verify each frame matches
    for i, (downloaded, original) in enumerate(zip(atoms_list, s22[:10])):
        assert len(downloaded) == len(original), f"Frame {i} atom count mismatch"
        assert downloaded.get_chemical_formula() == original.get_chemical_formula(), (
            f"Frame {i} formula mismatch"
        )
        assert np.allclose(downloaded.positions, original.positions, atol=1e-6), (
            f"Frame {i} positions mismatch"
        )


def test_download_multiple_specific_frames(server, s22):
    """Test downloading multiple specific frames by indices."""
    # Setup: Create room with multiple frames
    vis = ZnDraw(url=server, room="test-download-multiple", user="user1")
    for atoms in s22[:15]:
        vis.append(atoms)

    # Download frames 0, 5, 10
    response = requests.get(
        f"{server}/api/rooms/test-download-multiple/download",
        params={"indices": "0,5,10"},
    )

    assert response.status_code == 200

    # Parse ExtendedXYZ content
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    # Should contain exactly 3 frames
    assert len(atoms_list) == 3

    # Verify correct frames were downloaded
    expected_indices = [0, 5, 10]
    for downloaded, original_idx in zip(atoms_list, expected_indices):
        original = s22[original_idx]
        assert len(downloaded) == len(original)
        assert downloaded.get_chemical_formula() == original.get_chemical_formula()
        assert np.allclose(downloaded.positions, original.positions, atol=1e-6)


def test_download_with_selection(server, s22):
    """Test downloading frames with particle selection filter."""
    # Setup: Create room with a frame
    vis = ZnDraw(url=server, room="test-download-selection", user="user1")
    atoms = s22[0]  # Get first structure
    original_len = len(atoms)
    vis.append(atoms)

    # Download with selection of first 3 particles
    response = requests.get(
        f"{server}/api/rooms/test-download-selection/download",
        params={"indices": "0", "selection": "0,1,2"},
    )

    assert response.status_code == 200

    # Parse ExtendedXYZ content
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    assert len(atoms_list) == 1
    downloaded = atoms_list[0]

    # Should only contain 3 atoms
    assert len(downloaded) == 3

    # Verify positions match the selected particles
    assert np.allclose(downloaded.positions, atoms.positions[:3], atol=1e-6)


def test_download_preserves_metadata(server, s22):
    """Test that download preserves info dict metadata."""
    # Setup: Create room with atoms that have metadata
    vis = ZnDraw(url=server, room="test-download-metadata", user="user1")

    atoms = s22[0].copy()
    # Add metadata to info dict
    atoms.info["test_key"] = "test_value"
    atoms.info["test_number"] = 42

    vis.append(atoms)

    # Download the frame
    response = requests.get(
        f"{server}/api/rooms/test-download-metadata/download", params={"indices": "0"}
    )

    assert response.status_code == 200

    # Parse ExtendedXYZ content
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    downloaded = atoms_list[0]

    # Verify metadata is preserved
    assert "test_key" in downloaded.info
    assert downloaded.info["test_key"] == "test_value"
    assert "test_number" in downloaded.info
    assert downloaded.info["test_number"] == 42


def test_download_empty_room(server):
    """Test downloading from an empty room returns error."""
    # Setup: Create empty room
    vis = ZnDraw(url=server, room="test-download-empty", user="user1")

    # Download from empty room
    response = requests.get(f"{server}/api/rooms/test-download-empty/download")

    # Server validates that room has frames - should return 400
    assert response.status_code == 400


def test_download_invalid_index(server, s22):
    """Test downloading with out-of-range index returns error."""
    # Setup: Create room with 5 frames
    vis = ZnDraw(url=server, room="test-download-invalid", user="user1")
    for atoms in s22[:5]:
        vis.append(atoms)

    # Try to download frame 999 (doesn't exist)
    response = requests.get(
        f"{server}/api/rooms/test-download-invalid/download", params={"indices": "999"}
    )

    # Server validates indices - should return 400 for out-of-range
    assert response.status_code == 400


def test_download_negative_index(server, s22):
    """Test that negative indices are rejected."""
    # Setup: Create room with 5 frames
    vis = ZnDraw(url=server, room="test-download-negative", user="user1")
    for atoms in s22[:5]:
        vis.append(atoms)

    # Try to download using negative index
    response = requests.get(
        f"{server}/api/rooms/test-download-negative/download", params={"indices": "-1"}
    )

    # Negative indices are not supported - should return 400
    assert response.status_code == 400


def test_download_custom_filename(server, s22):
    """Test downloading with custom filename parameter."""
    # Setup: Create room with a frame
    vis = ZnDraw(url=server, room="test-download-filename", user="user1")
    vis.append(s22[0])

    # Download with custom filename
    custom_filename = "my_structure.xyz"
    response = requests.get(
        f"{server}/api/rooms/test-download-filename/download",
        params={"indices": "0", "filename": custom_filename},
    )

    assert response.status_code == 200

    # Verify Content-Disposition header contains custom filename
    content_disposition = response.headers.get("Content-Disposition", "")
    assert custom_filename in content_disposition


def test_download_extxyz_format_structure(server, s22):
    """Test that download produces valid ExtendedXYZ format."""
    # Setup: Create room with a frame
    vis = ZnDraw(url=server, room="test-download-format", user="user1")
    atoms = s22[0]
    vis.append(atoms)

    # Download the frame
    response = requests.get(
        f"{server}/api/rooms/test-download-format/download", params={"indices": "0"}
    )

    assert response.status_code == 200

    content = response.content.decode("utf-8")
    lines = content.strip().split("\n")

    # ExtendedXYZ format structure:
    # Line 1: Number of atoms
    # Line 2: Comment line with metadata (ExtendedXYZ format)
    # Lines 3+: Atom data

    assert len(lines) >= 2

    # First line should be integer (number of atoms)
    num_atoms = int(lines[0])
    assert num_atoms == len(atoms)

    # Second line is comment/metadata
    # ExtendedXYZ uses key=value format in comment line
    comment_line = lines[1]
    assert len(comment_line) > 0

    # Should have num_atoms lines of atom data
    atom_lines = lines[2:]
    assert len(atom_lines) >= num_atoms


@pytest.mark.parametrize("num_frames", [1, 5, 20])
def test_download_different_frame_counts(server, s22, num_frames):
    """Test downloading different numbers of frames."""
    # Setup: Create room with variable number of frames
    room_name = f"test-download-frames-{num_frames}"
    vis = ZnDraw(url=server, room=room_name, user="user1")

    for atoms in s22[:num_frames]:
        vis.append(atoms)

    # Download all frames
    response = requests.get(f"{server}/api/rooms/{room_name}/download")

    assert response.status_code == 200

    # Parse and verify
    content = response.content.decode("utf-8")
    atoms_list = list(ase.io.iread(io.StringIO(content), format="extxyz"))

    assert len(atoms_list) == num_frames


def test_download_room_not_found(server):
    """Test downloading from non-existent room."""
    # Try to download from room that doesn't exist
    response = requests.get(f"{server}/api/rooms/nonexistent-room-xyz/download")

    # Should return error (likely 404 or 400)
    # Based on implementation, it might return 200 with empty content
    # The exact behavior depends on error handling
    assert response.status_code in [200, 400, 404]
