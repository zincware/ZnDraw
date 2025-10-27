import pytest

from zndraw.zndraw import ZnDraw


def test_vis_bookmarks(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis.bookmarks) == 0
    vis.extend(s22)

    # Use new MutableMapping interface
    vis.bookmarks[1] = "First Frame"
    vis.bookmarks[5] = "Middle Frame"
    vis.bookmarks[9] = "Last Frame"

    # Convert to dict for comparison
    assert dict(vis.bookmarks) == {1: "First Frame", 5: "Middle Frame", 9: "Last Frame"}

    del vis[2]
    assert dict(vis.bookmarks) == {1: "First Frame", 4: "Middle Frame", 8: "Last Frame"}
    del vis[1]
    assert dict(vis.bookmarks) == {3: "Middle Frame", 7: "Last Frame"}
    vis.insert(0, s22[0])
    assert dict(vis.bookmarks) == {4: "Middle Frame", 8: "Last Frame"}
    # replacing the frame removes the bookmark as well
    vis[4] = s22[3]
    assert dict(vis.bookmarks) == {8: "Last Frame"}

    # Clear all bookmarks by deleting them one by one
    for idx in list(vis.bookmarks):
        del vis.bookmarks[idx]
    assert len(vis.bookmarks) == 0


def test_vis_bookmarks_errors(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(IndexError):
        vis.bookmarks[999] = "Out of bounds"

    with pytest.raises(TypeError):
        vis.bookmarks[0] = 123  # Non-string label


def test_vis_bookmarks_multiple_clients(server, s22):
    """Test that bookmarks are room-wide and synced across multiple clients."""
    vis1 = ZnDraw(url=server, room="testroom", user="user1")
    vis2 = ZnDraw(url=server, room="testroom", user="user2")

    vis1.extend(s22)
    vis1.bookmarks[1] = "First Frame"
    vis1.bookmarks[5] = "Middle Frame"
    vis1.bookmarks[9] = "Last Frame"

    assert len(vis2) == len(s22)
    assert len(vis1) == len(s22)

    # Wait a moment for socket events to propagate
    import time

    time.sleep(0.1)

    # Both clients should see the same bookmarks (room-wide)
    assert (
        dict(vis2.bookmarks)
        == dict(vis1.bookmarks)
        == {
            1: "First Frame",
            5: "Middle Frame",
            9: "Last Frame",
        }
    )

    # vis2 adds new bookmarks and updates one
    vis2.bookmarks[2] = "Second Frame"
    vis2.bookmarks[6] = "Another Middle"
    vis2.bookmarks[9] = "Updated Last Frame"

    # Wait for propagation
    time.sleep(0.1)

    # Bookmarks are room-wide, so both clients see ALL bookmarks
    expected_bookmarks = {
        1: "First Frame",
        2: "Second Frame",
        5: "Middle Frame",
        6: "Another Middle",
        9: "Updated Last Frame",
    }
    assert dict(vis2.bookmarks) == expected_bookmarks
    assert dict(vis1.bookmarks) == expected_bookmarks


def test_vis_bookmarks_update(server, s22):
    """Test the update() method for batch bookmark updates."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Test update with dict
    vis.bookmarks.update({0: "First", 1: "Second", 2: "Third"})
    assert dict(vis.bookmarks) == {0: "First", 1: "Second", 2: "Third"}

    # Test update with iterable of tuples
    vis.bookmarks.update([(3, "Fourth"), (4, "Fifth")])
    assert dict(vis.bookmarks) == {
        0: "First",
        1: "Second",
        2: "Third",
        3: "Fourth",
        4: "Fifth",
    }

    # Test update overwriting existing bookmarks
    vis.bookmarks.update({0: "Updated First", 2: "Updated Third"})
    assert dict(vis.bookmarks) == {
        0: "Updated First",
        1: "Second",
        2: "Updated Third",
        3: "Fourth",
        4: "Fifth",
    }

    # Test update with empty dict (should do nothing)
    vis.bookmarks.update({})
    assert len(vis.bookmarks) == 5


def test_vis_bookmarks_clear(server, s22):
    """Test the clear() method to remove all bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Add several bookmarks
    vis.bookmarks.update({0: "First", 1: "Second", 2: "Third", 3: "Fourth"})
    assert len(vis.bookmarks) == 4

    # Clear all bookmarks
    vis.bookmarks.clear()
    assert len(vis.bookmarks) == 0
    assert dict(vis.bookmarks) == {}

    # Clear again (should be safe on empty)
    vis.bookmarks.clear()
    assert len(vis.bookmarks) == 0


def test_vis_bookmarks_update_errors(server, s22):
    """Test that update() properly validates inputs."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Test update with out-of-range index
    with pytest.raises(IndexError):
        vis.bookmarks.update({999: "Out of bounds"})

    # Test update with non-string label
    with pytest.raises(TypeError):
        vis.bookmarks.update({0: 123})

    # Test update with empty string label
    with pytest.raises(ValueError):
        vis.bookmarks.update({0: ""})


def test_vis_bookmarks_dict_methods(server, s22):
    """Test that bookmarks support standard dict methods."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Set up some bookmarks
    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})

    # Test keys()
    assert set(vis.bookmarks.keys()) == {0, 2, 4}

    # Test values()
    assert set(vis.bookmarks.values()) == {"First", "Second", "Third"}

    # Test items()
    assert set(vis.bookmarks.items()) == {(0, "First"), (2, "Second"), (4, "Third")}

    # Test get() with default
    assert vis.bookmarks.get(0) == "First"
    assert vis.bookmarks.get(1) is None
    assert vis.bookmarks.get(1, "Default") == "Default"

    # Test 'in' operator
    assert 0 in vis.bookmarks
    assert 1 not in vis.bookmarks


def test_vis_bookmarks_pop(server, s22):
    """Test the pop() method for removing and returning bookmarks."""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Set up some bookmarks
    vis.bookmarks.update({0: "First", 2: "Second", 4: "Third"})
    assert len(vis.bookmarks) == 3

    # Test pop() with existing bookmark
    label = vis.bookmarks.pop(2)
    assert label == "Second"
    assert len(vis.bookmarks) == 2
    assert 2 not in vis.bookmarks
    assert dict(vis.bookmarks) == {0: "First", 4: "Third"}

    # Test pop() with non-existing bookmark and default
    label = vis.bookmarks.pop(99, "Default")
    assert label == "Default"
    assert len(vis.bookmarks) == 2  # No change

    # Test pop() with non-existing bookmark without default (should raise KeyError)
    with pytest.raises(KeyError):
        vis.bookmarks.pop(99)

    # Test pop() removes bookmark from cache
    vis.bookmarks.pop(0)
    assert 0 not in vis.bookmarks
    assert len(vis.bookmarks) == 1
