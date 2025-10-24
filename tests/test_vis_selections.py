import time

import pytest

from zndraw.zndraw import ZnDraw


# ==================== Per-Geometry Selections Tests ====================


def test_selections_init_empty(s22, server):
    """Test that selections are initially empty."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis.selections) == 0


def test_selections_set_get(s22, server):
    """Test setting and getting selections for a specific geometry."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Set selection for "particles"
    vis.selections["particles"] = [0, 1, 2]
    assert vis.selections["particles"] == (0, 1, 2)

    # Set selection for another geometry
    vis.selections["forces"] = [2, 3, 4]
    assert vis.selections["forces"] == (2, 3, 4)

    # Verify both selections exist
    assert "particles" in vis.selections
    assert "forces" in vis.selections


def test_selections_delete(s22, server):
    """Test deleting selections for a geometry."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.selections["particles"] = [0, 1, 2]
    assert vis.selections["particles"] == (0, 1, 2)

    # Delete the selection
    del vis.selections["particles"]
    assert vis.selections["particles"] == ()


def test_selections_iteration(s22, server):
    """Test iterating over geometries with selections."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.selections["particles"] = [0, 1]
    vis.selections["forces"] = [2, 3]
    vis.selections["velocities"] = [4, 5]

    geometries = list(vis.selections)
    assert "particles" in geometries
    assert "forces" in geometries
    assert "velocities" in geometries


def test_selections_len(s22, server):
    """Test getting the number of geometries with selections."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    assert len(vis.selections) == 0

    vis.selections["particles"] = [0, 1]
    assert len(vis.selections) == 1

    vis.selections["forces"] = [2, 3]
    assert len(vis.selections) == 2


def test_selections_persistence(s22, server):
    """Test that selections persist across clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)

    # Set selection in vis1
    vis1.selections["particles"] = [0, 1, 2]
    time.sleep(0.2)  # Wait for invalidation

    # Check selection in vis2
    assert vis2.selections["particles"] == (0, 1, 2)

    # Update in vis2
    vis2.selections["particles"] = [3, 4]
    time.sleep(0.2)

    # Check update in vis1
    assert vis1.selections["particles"] == (3, 4)


def test_selections_different_rooms(s22, server):
    """Test that selections are isolated between different rooms."""
    vis1 = ZnDraw(url=server, room="room1", user="testuser1")
    vis2 = ZnDraw(url=server, room="room2", user="testuser2")
    vis1.extend(s22)
    vis2.extend(s22)

    vis1.selections["particles"] = [0, 1]
    vis2.selections["particles"] = [2, 3]

    assert vis1.selections["particles"] == (0, 1)
    assert vis2.selections["particles"] == (2, 3)


# ==================== Selection Groups Tests ====================


def test_selection_groups_init_empty(s22, server):
    """Test that selection groups are initially empty."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis.selection_groups) == 0


def test_selection_groups_create(s22, server):
    """Test creating a selection group."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create a group
    group_data = {
        "particles": [0, 1, 2],
        "forces": [3, 4, 5]
    }
    vis.selection_groups["group1"] = group_data

    # Verify the group
    retrieved = vis.selection_groups["group1"]
    assert retrieved["particles"] == [0, 1, 2]
    assert retrieved["forces"] == [3, 4, 5]


def test_selection_groups_update(s22, server):
    """Test updating a selection group."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create a group
    vis.selection_groups["group1"] = {"particles": [0, 1]}

    # Update the group
    vis.selection_groups["group1"] = {"particles": [2, 3, 4]}

    # Verify the update
    retrieved = vis.selection_groups["group1"]
    assert retrieved["particles"] == [2, 3, 4]


def test_selection_groups_delete(s22, server):
    """Test deleting a selection group."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create a group
    vis.selection_groups["group1"] = {"particles": [0, 1]}
    assert "group1" in vis.selection_groups

    # Delete the group
    del vis.selection_groups["group1"]

    # Verify it's gone
    with pytest.raises(Exception):  # Should raise an error when accessing non-existent group
        _ = vis.selection_groups["group1"]


def test_selection_groups_iteration(s22, server):
    """Test iterating over selection groups."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.selection_groups["group1"] = {"particles": [0, 1]}
    vis.selection_groups["group2"] = {"forces": [2, 3]}
    vis.selection_groups["group3"] = {"velocities": [4, 5]}

    groups = list(vis.selection_groups)
    assert "group1" in groups
    assert "group2" in groups
    assert "group3" in groups


def test_selection_groups_len(s22, server):
    """Test getting the number of selection groups."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    assert len(vis.selection_groups) == 0

    vis.selection_groups["group1"] = {"particles": [0, 1]}
    assert len(vis.selection_groups) == 1

    vis.selection_groups["group2"] = {"forces": [2, 3]}
    assert len(vis.selection_groups) == 2


def test_selection_groups_persistence(s22, server):
    """Test that selection groups persist across clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)

    # Create group in vis1
    vis1.selection_groups["group1"] = {"particles": [0, 1, 2]}
    time.sleep(0.1)

    # Check in vis2
    group_data = vis2.selection_groups["group1"]
    assert group_data["particles"] == [0, 1, 2]


# ==================== Active Selection Group Tests ====================


def test_active_selection_group_init_none(s22, server):
    """Test that active selection group is initially None."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert vis.active_selection_group is None


def test_load_selection_group(s22, server):
    """Test loading a selection group."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create a group
    group_data = {
        "particles": [0, 1, 2],
        "forces": [3, 4, 5]
    }
    vis.selection_groups["mygroup"] = group_data

    # Load the group
    vis.load_selection_group("mygroup")
    time.sleep(0.2)  # Wait for server update and invalidation

    # Verify selections were updated
    assert vis.selections["particles"] == (0, 1, 2)
    assert vis.selections["forces"] == (3, 4, 5)

    # Verify active group was set
    assert vis.active_selection_group == "mygroup"


def test_load_selection_group_clears_active_on_manual_change(s22, server):
    """Test that manually changing a selection clears the active group."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create and load a group
    vis.selection_groups["mygroup"] = {"particles": [0, 1, 2]}
    vis.load_selection_group("mygroup")
    time.sleep(0.2)

    assert vis.active_selection_group == "mygroup"

    # Manually change a selection
    vis.selections["particles"] = [5, 6, 7]
    time.sleep(0.2)

    # Active group should be cleared
    assert vis.active_selection_group is None


def test_load_selection_group_persistence(s22, server):
    """Test that loaded selection groups sync across clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)

    # Create and load group in vis1
    vis1.selection_groups["group1"] = {"particles": [0, 1, 2]}
    vis1.load_selection_group("group1")
    time.sleep(0.3)

    # Check selections in vis2
    assert vis2.selections["particles"] == (0, 1, 2)
    assert vis2.active_selection_group == "group1"


def test_selection_particles(s22, server):
    """Test that old 'selection' property still works and maps to 'particles'."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Set using old property
    vis.selection = [0, 1, 2]
    time.sleep(0.1)

    # Should update both old and new
    assert vis.selection == (0, 1, 2)
    assert vis.selections["particles"] == (0, 1, 2)

    # Set using new property
    vis.selections["particles"] = [3, 4, 5]
    time.sleep(0.1)

    # Should update old property too
    assert vis.selection == (3, 4, 5)


def test_selection_particles_persistence(s22, server):
    """Test that old and new selection APIs work together across clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)

    # Set using old API in vis1
    vis1.selection = [0, 1, 2]
    time.sleep(0.2)

    # Check using new API in vis2
    assert vis2.selections["particles"] == (0, 1, 2)

    # Set using new API in vis2
    vis2.selections["particles"] = [5, 6]
    time.sleep(0.2)

    # Check using old API in vis1
    assert vis1.selection == (5, 6)
