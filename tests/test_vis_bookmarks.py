import time

import pytest

from zndraw.zndraw import ZnDraw

# ==================== Basic Bookmark Operations ====================


def test_vis_bookmarks_init_empty(server):
    """Test that bookmarks are initially empty."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis.bookmarks) == 0
    assert dict(vis.bookmarks) == {}


def test_vis_bookmarks_set(server, s22):
    """Test setting a bookmark at a frame index."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "First Frame"
    assert vis.bookmarks[1] == "First Frame"


def test_vis_bookmarks_set_multiple(server, s22):
    """Test setting multiple bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "First Frame"
    vis.bookmarks[5] = "Middle Frame"
    vis.bookmarks[9] = "Last Frame"

    assert dict(vis.bookmarks) == {1: "First Frame", 5: "Middle Frame", 9: "Last Frame"}


def test_vis_bookmarks_overwrite(server, s22):
    """Test overwriting an existing bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "Original"
    vis.bookmarks[1] = "Updated"

    assert vis.bookmarks[1] == "Updated"


def test_vis_bookmarks_delete(server, s22):
    """Test deleting a bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "First Frame"
    assert 1 in vis.bookmarks

    del vis.bookmarks[1]
    assert 1 not in vis.bookmarks


def test_vis_bookmarks_len(server, s22):
    """Test bookmark count."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    assert len(vis.bookmarks) == 0

    vis.bookmarks[0] = "First"
    assert len(vis.bookmarks) == 1

    vis.bookmarks[1] = "Second"
    assert len(vis.bookmarks) == 2

    del vis.bookmarks[0]
    assert len(vis.bookmarks) == 1


# ==================== Index Adjustment Tests ====================


def test_vis_bookmarks_adjust_on_frame_delete(server, s22):
    """Test that bookmark indices shift when a frame before them is deleted."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "First Frame"
    vis.bookmarks[5] = "Middle Frame"
    vis.bookmarks[9] = "Last Frame"

    del vis[2]

    assert dict(vis.bookmarks) == {1: "First Frame", 4: "Middle Frame", 8: "Last Frame"}


def test_vis_bookmarks_adjust_on_early_frame_delete(server, s22):
    """Test that deleting a frame before all bookmarks shifts all indices."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[3] = "Middle Frame"
    vis.bookmarks[7] = "Last Frame"

    del vis[1]

    assert dict(vis.bookmarks) == {2: "Middle Frame", 6: "Last Frame"}


def test_vis_bookmarks_removed_on_bookmarked_frame_delete(server, s22):
    """Test that deleting a bookmarked frame removes that bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[1] = "First Frame"
    vis.bookmarks[5] = "Middle Frame"

    del vis[1]

    assert 1 not in vis.bookmarks
    assert dict(vis.bookmarks) == {4: "Middle Frame"}


def test_vis_bookmarks_adjust_on_frame_insert(server, s22):
    """Test that bookmark indices shift when a frame is inserted before them."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[4] = "Middle Frame"
    vis.bookmarks[8] = "Last Frame"

    vis.insert(0, s22[0])

    assert dict(vis.bookmarks) == {5: "Middle Frame", 9: "Last Frame"}


def test_vis_bookmarks_removed_on_frame_replace(server, s22):
    """Test that replacing a bookmarked frame removes its bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[4] = "Middle Frame"
    vis.bookmarks[8] = "Last Frame"

    vis[4] = s22[3]

    assert 4 not in vis.bookmarks
    assert dict(vis.bookmarks) == {8: "Last Frame"}


# ==================== Error Handling Tests ====================


@pytest.mark.parametrize(
    "index,label,expected_error",
    [
        (999, "Out of bounds", IndexError),
        (0, 123, TypeError),
        (0, "", ValueError),
    ],
    ids=["out_of_bounds", "non_string_label", "empty_label"],
)
def test_vis_bookmarks_set_errors(server, s22, index, label, expected_error):
    """Test that setting bookmarks validates inputs correctly."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(expected_error):
        vis.bookmarks[index] = label


def test_vis_bookmarks_get_nonexistent(server, s22):
    """Test that getting a nonexistent bookmark raises KeyError."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(KeyError):
        _ = vis.bookmarks[0]


# ==================== Multi-Client Sync Tests ====================


def test_vis_bookmarks_sync_to_second_client(server, s22):
    """Test that bookmarks set by one client are visible to another."""
    vis1 = ZnDraw(url=server, room="testroom", user="user1")
    vis2 = ZnDraw(url=server, room="testroom", user="user2")

    vis1.extend(s22)
    vis1.bookmarks[1] = "First Frame"
    vis1.bookmarks[5] = "Middle Frame"

    time.sleep(0.1)

    assert dict(vis2.bookmarks) == {1: "First Frame", 5: "Middle Frame"}


def test_vis_bookmarks_sync_bidirectional(server, s22):
    """Test that bookmarks sync bidirectionally between clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="user1")
    vis2 = ZnDraw(url=server, room="testroom", user="user2")

    vis1.extend(s22)
    vis1.bookmarks[1] = "From User1"

    time.sleep(0.1)

    vis2.bookmarks[2] = "From User2"

    time.sleep(0.1)

    expected = {1: "From User1", 2: "From User2"}
    assert dict(vis1.bookmarks) == expected
    assert dict(vis2.bookmarks) == expected


def test_vis_bookmarks_overwrite_syncs(server, s22):
    """Test that overwriting a bookmark syncs to other clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="user1")
    vis2 = ZnDraw(url=server, room="testroom", user="user2")

    vis1.extend(s22)
    vis1.bookmarks[1] = "Original"

    time.sleep(0.1)

    vis2.bookmarks[1] = "Updated"

    time.sleep(0.1)

    assert vis1.bookmarks[1] == "Updated"
    assert vis2.bookmarks[1] == "Updated"


# ==================== Update Method Tests ====================


def test_vis_bookmarks_update_with_dict(server, s22):
    """Test update() with a dictionary."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 1: "Second", 2: "Third"})

    assert dict(vis.bookmarks) == {0: "First", 1: "Second", 2: "Third"}


def test_vis_bookmarks_update_with_tuples(server, s22):
    """Test update() with an iterable of tuples."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update([(0, "First"), (1, "Second")])

    assert dict(vis.bookmarks) == {0: "First", 1: "Second"}


def test_vis_bookmarks_update_overwrites(server, s22):
    """Test that update() overwrites existing bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[0] = "Original"
    vis.bookmarks.update({0: "Updated"})

    assert vis.bookmarks[0] == "Updated"


def test_vis_bookmarks_update_empty(server, s22):
    """Test that update() with empty dict does nothing."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[0] = "First"
    vis.bookmarks.update({})

    assert len(vis.bookmarks) == 1
    assert vis.bookmarks[0] == "First"


@pytest.mark.parametrize(
    "update_data,expected_error",
    [
        ({999: "Out of bounds"}, IndexError),
        ({0: 123}, TypeError),
        ({0: ""}, ValueError),
    ],
    ids=["out_of_bounds", "non_string_label", "empty_label"],
)
def test_vis_bookmarks_update_errors(server, s22, update_data, expected_error):
    """Test that update() properly validates inputs."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(expected_error):
        vis.bookmarks.update(update_data)


# ==================== Clear Method Tests ====================


def test_vis_bookmarks_clear(server, s22):
    """Test clear() removes all bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 1: "Second", 2: "Third"})
    assert len(vis.bookmarks) == 3

    vis.bookmarks.clear()

    assert len(vis.bookmarks) == 0
    assert dict(vis.bookmarks) == {}


def test_vis_bookmarks_clear_empty(server, s22):
    """Test that clear() is safe on empty bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.clear()

    assert len(vis.bookmarks) == 0


# ==================== Dict Method Tests ====================


def test_vis_bookmarks_keys(server, s22):
    """Test keys() returns bookmark indices."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})

    assert set(vis.bookmarks.keys()) == {0, 2, 4}


def test_vis_bookmarks_values(server, s22):
    """Test values() returns bookmark labels."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})

    assert set(vis.bookmarks.values()) == {"First", "Second", "Third"}


def test_vis_bookmarks_items(server, s22):
    """Test items() returns index-label pairs."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})

    assert set(vis.bookmarks.items()) == {(0, "First"), (2, "Second"), (4, "Third")}


def test_vis_bookmarks_get_existing(server, s22):
    """Test get() returns value for existing bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[0] = "First"

    assert vis.bookmarks.get(0) == "First"


def test_vis_bookmarks_get_nonexistent_default(server, s22):
    """Test get() returns None for nonexistent bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    assert vis.bookmarks.get(0) is None


def test_vis_bookmarks_get_with_default(server, s22):
    """Test get() returns provided default for nonexistent bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    assert vis.bookmarks.get(0, "Default") == "Default"


def test_vis_bookmarks_contains(server, s22):
    """Test 'in' operator for bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[0] = "First"

    assert 0 in vis.bookmarks
    assert 1 not in vis.bookmarks


def test_vis_bookmarks_iteration(server, s22):
    """Test iterating over bookmarks yields indices."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})

    assert set(vis.bookmarks) == {0, 2, 4}


# ==================== Pop Method Tests ====================


def test_vis_bookmarks_pop_existing(server, s22):
    """Test pop() removes and returns existing bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    vis.bookmarks[0] = "First"

    label = vis.bookmarks.pop(0)

    assert label == "First"
    assert 0 not in vis.bookmarks


def test_vis_bookmarks_pop_with_default(server, s22):
    """Test pop() returns default for nonexistent bookmark."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    label = vis.bookmarks.pop(99, "Default")

    assert label == "Default"


def test_vis_bookmarks_pop_nonexistent_raises(server, s22):
    """Test pop() raises KeyError for nonexistent bookmark without default."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(KeyError):
        vis.bookmarks.pop(99)
