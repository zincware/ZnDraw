"""Tests for deleting template frames from rooms.

Template frames are references to frames in another room. Deleting them
should remove the reference from the index, not the physical data.
"""

import uuid

import numpy as np

from zndraw import ZnDraw


def test_delete_copied_room_frames(server):
    """Deleting frames in a room copied from another should work."""
    # Create a source room with some frames
    source_room = f"source-{uuid.uuid4().hex}"
    source = ZnDraw(room=source_room, url=server)

    for i in range(5):
        with source.get_lock():
            source._append_frame({"index": np.array([i])})

    assert len(source) == 5
    source.disconnect()

    # Create a new room copied from source (uses copyFrom)
    user_room = f"user-{uuid.uuid4().hex}"
    user = ZnDraw(room=user_room, url=server, copy_from=source_room)

    # User room should have the copied frames
    assert len(user) == 5

    # Deleting frames should work - removes from index, not physical data
    del user[:3]

    assert len(user) == 2

    # Remaining frames should be the last two from source
    frame = user.get(0)
    assert np.array_equal(frame["index"], np.array([3]))

    frame = user.get(1)
    assert np.array_equal(frame["index"], np.array([4]))

    user.disconnect()

    # Source room should be unaffected
    source2 = ZnDraw(room=source_room, url=server)
    assert len(source2) == 5
    source2.disconnect()


def test_delete_single_copied_frame(server):
    """Deleting a single frame in a copied room should work."""
    source_room = f"source-{uuid.uuid4().hex}"
    source = ZnDraw(room=source_room, url=server)

    for i in range(3):
        with source.get_lock():
            source._append_frame({"index": np.array([i])})

    source.disconnect()

    user_room = f"user-{uuid.uuid4().hex}"
    user = ZnDraw(room=user_room, url=server, copy_from=source_room)

    assert len(user) == 3

    # Delete middle frame
    del user[1]

    assert len(user) == 2
    frame0 = user.get(0)
    frame1 = user.get(1)
    assert np.array_equal(frame0["index"], np.array([0]))
    assert np.array_equal(frame1["index"], np.array([2]))

    user.disconnect()


def test_delete_all_copied_frames(server):
    """Deleting all frames in a copied room should result in empty room."""
    source_room = f"source-{uuid.uuid4().hex}"
    source = ZnDraw(room=source_room, url=server)

    for i in range(3):
        with source.get_lock():
            source._append_frame({"index": np.array([i])})

    source.disconnect()

    user_room = f"user-{uuid.uuid4().hex}"
    user = ZnDraw(room=user_room, url=server, copy_from=source_room)

    assert len(user) == 3

    del user[:]

    assert len(user) == 0

    user.disconnect()


def test_delete_mixed_frames(server):
    """Deleting works when room has both copied and local frames."""
    source_room = f"source-{uuid.uuid4().hex}"
    source = ZnDraw(room=source_room, url=server)

    for i in range(2):
        with source.get_lock():
            source._append_frame({"index": np.array([i]), "source": "template"})

    source.disconnect()

    user_room = f"user-{uuid.uuid4().hex}"
    user = ZnDraw(room=user_room, url=server, copy_from=source_room)

    # Add local frames
    for i in range(2):
        with user.get_lock():
            user._append_frame({"index": np.array([100 + i]), "source": "local"})

    assert len(user) == 4

    # Delete first copied frame
    del user[0]
    assert len(user) == 3

    # Delete what was at index 2 (now at 1) - a local frame
    del user[1]
    assert len(user) == 2

    user.disconnect()
