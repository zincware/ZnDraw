from zndraw_communication import Client
import uuid
import numpy as np



def test_connection():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()
    assert client.sio.connected
    client.disconnect()
    assert not client.sio.connected

def test_len_frames_empty():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()
    assert client.len_frames() == 0
    client.disconnect()
    assert not client.sio.connected

def test_append_and_get_frame():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()
    for i in range(10):
        data = {
            "index": np.array([i]),
            "points": np.random.rand(5, 3),
            "colors": np.random.randint(0, 255, size=(5, 4), dtype=np.uint8)
        }
        client.append_frame(data)
        assert client.len_frames() == i + 1
        frame = client.get_frame(i)
        assert np.array_equal(frame["index"], data["index"])
        assert np.array_equal(frame["points"], data["points"])
        assert np.array_equal(frame["colors"], data["colors"])
    client.disconnect()
    assert not client.sio.connected

def test_delete_frame():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()
    for i in range(5):
        client.append_frame({"index": np.array([i])})
    assert client.len_frames() == 5
    client.delete_frame(0)
    assert client.len_frames() == 4
    frame = client.get_frame(0)
    assert np.array_equal(frame["index"], np.array([1]))
    client.delete_frame(1)
    assert client.len_frames() == 3
    frame = client.get_frame(1)
    assert np.array_equal(frame["index"], np.array([3]))
    client.disconnect()
    assert not client.sio.connected

def test_replace_frame():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()
    for i in range(10):
        client.append_frame({"index": np.array([i])})
    assert client.len_frames() == 10

    # Replace frame at index 5
    new_data = {
        "index": np.array([999]),
    }
    client.replace_frame(5, new_data)
    frame = client.get_frame(5)
    assert np.array_equal(frame["index"], new_data["index"])
    assert client.len_frames() == 10

    client.disconnect()
    assert not client.sio.connected

def test_extend_frames():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Start with a few frames
    for i in range(3):
        client.append_frame({"index": np.array([i])})
    assert client.len_frames() == 3

    # Extend with multiple frames at once
    extend_data = [
        {"index": np.array([10]), "points": np.random.rand(5, 3)},
        {"index": np.array([20]), "points": np.random.rand(5, 3)},
        {"index": np.array([30]), "points": np.random.rand(5, 3)},
    ]

    new_indices = client.extend_frames(extend_data)
    assert client.len_frames() == 6
    assert new_indices == [3, 4, 5]  # Should be appended at these logical positions

    # Verify the extended frames
    for i, expected_idx in enumerate([10, 20, 30]):
        frame = client.get_frame(3 + i)
        assert np.array_equal(frame["index"], np.array([expected_idx]))
        assert frame["points"].shape == (5, 3)

    client.disconnect()
    assert not client.sio.connected

def test_get_frames_with_indices():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Add some test frames
    test_data = []
    for i in range(10):
        data = {
            "index": np.array([i * 10]),
            "points": np.random.rand(3, 3)
        }
        client.append_frame(data)
        test_data.append(data)

    # Test fetching specific frames by indices
    indices = [0, 2, 5, 8]
    frames = client.get_frames(indices)

    assert len(frames) == len(indices)
    for i, frame_idx in enumerate(indices):
        assert np.array_equal(frames[i]["index"], test_data[frame_idx]["index"])
        assert np.array_equal(frames[i]["points"], test_data[frame_idx]["points"])

    client.disconnect()
    assert not client.sio.connected

def test_get_frames_with_slice():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Add some test frames
    test_data = []
    for i in range(20):
        data = {"index": np.array([i])}
        client.append_frame(data)
        test_data.append(data)

    # Test fetching frames with slice notation
    frames = client.get_frames(slice(5, 15, 2))  # Every 2nd frame from 5 to 15
    expected_indices = [5, 7, 9, 11, 13]

    assert len(frames) == len(expected_indices)
    for i, expected_idx in enumerate(expected_indices):
        assert np.array_equal(frames[i]["index"], test_data[expected_idx]["index"])

    # Test slice with just start
    frames = client.get_frames(slice(10, None))  # From 10 to end
    assert len(frames) == 10  # frames 10-19
    for i, frame in enumerate(frames):
        assert np.array_equal(frame["index"], test_data[10 + i]["index"])

    # Test slice with step only
    frames = client.get_frames(slice(None, None, 3))  # Every 3rd frame
    expected_count = (20 + 2) // 3  # Ceiling division
    assert len(frames) == expected_count

    # Test slice(None, None, None) - should get all frames
    frames = client.get_frames(slice(None, None, None))
    assert len(frames) == 20  # All frames
    for i, frame in enumerate(frames):
        assert np.array_equal(frame["index"], test_data[i]["index"])

    client.disconnect()
    assert not client.sio.connected

def test_get_frames_empty_result():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Add a few frames
    for i in range(5):
        client.append_frame({"index": np.array([i])})

    # Test slice that results in no frames
    frames = client.get_frames(slice(10, 20))  # Beyond available frames
    assert len(frames) == 0

    # Test empty indices list
    frames = client.get_frames([])
    assert len(frames) == 0

    client.disconnect()
    assert not client.sio.connected

def test_mutable_sequence_interface():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Test len()
    assert len(client) == 0

    # Test append()
    frame1 = {"index": np.array([10])}
    frame2 = {"index": np.array([20])}
    client.append(frame1)
    assert len(client) == 1
    client.append(frame2)
    assert len(client) == 2

    # Test __getitem__ with positive index
    retrieved = client[0]
    assert np.array_equal(retrieved["index"], frame1["index"])

    # Test __getitem__ with negative index
    retrieved = client[-1]
    assert np.array_equal(retrieved["index"], frame2["index"])

    # Test __getitem__ with slice
    frames = client[0:2]
    assert len(frames) == 2
    assert np.array_equal(frames[0]["index"], frame1["index"])
    assert np.array_equal(frames[1]["index"], frame2["index"])

    # Test __setitem__ (replace)
    new_frame = {"index": np.array([999])}
    client[1] = new_frame
    retrieved = client[1]
    assert np.array_equal(retrieved["index"], new_frame["index"])

    # Test extend()
    extend_data = [
        {"index": np.array([100])},
        {"index": np.array([200])}
    ]
    client.extend(extend_data)
    assert len(client) == 4
    assert np.array_equal(client[2]["index"], extend_data[0]["index"])
    assert np.array_equal(client[3]["index"], extend_data[1]["index"])

    # Test __delitem__
    client[1]  # Should be the replaced frame with index 999
    del client[1]
    assert len(client) == 3
    # After deletion, what was at index 2 should now be at index 1
    assert np.array_equal(client[1]["index"], extend_data[0]["index"])

    client.disconnect()
    assert not client.sio.connected

def test_insert_frame_functionality():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Start with some initial frames
    initial_frames = [
        {"index": np.array([0]), "value": np.array([100])},
        {"index": np.array([1]), "value": np.array([200])},
        {"index": np.array([2]), "value": np.array([300])},
    ]
    for frame in initial_frames:
        client.append_frame(frame)

    assert len(client) == 3

    # Test insert at beginning
    insert_frame_0 = {"index": np.array([999]), "value": np.array([999])}
    client.insert_frame(0, insert_frame_0)
    assert len(client) == 4

    # Verify the insertion shifted everything
    frame = client.get_frame(0)
    assert np.array_equal(frame["index"], insert_frame_0["index"])
    frame = client.get_frame(1)  # Should be the original frame 0
    assert np.array_equal(frame["index"], initial_frames[0]["index"])

    # Test insert in middle
    insert_frame_2 = {"index": np.array([888]), "value": np.array([888])}
    client.insert_frame(2, insert_frame_2)
    assert len(client) == 5

    # Verify the insertion
    frame = client.get_frame(2)
    assert np.array_equal(frame["index"], insert_frame_2["index"])
    # Original frame 1 should now be at position 3
    frame = client.get_frame(3)
    assert np.array_equal(frame["index"], initial_frames[1]["index"])

    # Test insert at end (equivalent to append)
    insert_frame_end = {"index": np.array([777]), "value": np.array([777])}
    client.insert_frame(len(client), insert_frame_end)
    assert len(client) == 6

    # Should be at the last position
    frame = client.get_frame(-1)
    assert np.array_equal(frame["index"], insert_frame_end["index"])

    client.disconnect()
    assert not client.sio.connected

def test_mutable_sequence_insert():
    client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
    client.connect()

    # Add initial frames using MutableSequence interface
    frames = [
        {"data": np.array([1, 2, 3])},
        {"data": np.array([4, 5, 6])},
        {"data": np.array([7, 8, 9])},
    ]

    for frame in frames:
        client.append(frame)

    assert len(client) == 3

    # Test MutableSequence insert at position 1
    insert_data = {"data": np.array([99, 99, 99])}
    client.insert(1, insert_data)
    assert len(client) == 4

    # Verify the order
    assert np.array_equal(client[0]["data"], frames[0]["data"])  # Original frame 0
    assert np.array_equal(client[1]["data"], insert_data["data"])  # Inserted frame
    assert np.array_equal(client[2]["data"], frames[1]["data"])   # Shifted from pos 1 to 2
    assert np.array_equal(client[3]["data"], frames[2]["data"])   # Shifted from pos 2 to 3

    client.disconnect()
    assert not client.sio.connected


# def test_replace_frame_additional_keys():
#     client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
#     client.connect()
#     for i in range(10):
#         client.append_frame({"index": np.array([i])})
#     assert client.len_frames() == 10

#     # Replace frame at index 5 with additional keys
#     new_data = {
#         "index": np.array([999]),
#         "points": np.random.rand(5, 3),
#         "colors": np.random.randint(0, 255, size=(5, 4), dtype=np.uint8)
#     }
#     client.replace_frame(5, new_data)
#     frame = client.get_frame(5)
#     assert np.array_equal(frame["index"], new_data["index"])
#     assert np.array_equal(frame["points"], new_data["points"])
#     assert np.array_equal(frame["colors"], new_data["colors"])
#     assert client.len_frames() == 10

#     client.disconnect()
#     assert not client.sio.connected

# replace with data of different shape
# def test_replace_frame_different_shape():
#     client = Client(room=uuid.uuid4().hex, url="http://localhost:5000")
#     client.connect()
#     for i in range(10):
#         client.append_frame({"index": np.array([i]), "points": np.random.rand(5, 3)})
#     assert client.len_frames() == 10

#     # Replace frame at index 5 with different shape
#     new_data = {
#         "index": np.array([999]),
#         "points": np.random.rand(10, 3)  # Different shape
#     }
#     client.replace_frame(5, new_data)
#     frame = client.get_frame(5)
#     assert np.array_equal(frame["index"], new_data["index"])
#     assert np.array_equal(frame["points"], new_data["points"])
#     assert client.len_frames() == 10

#     client.disconnect()
#     assert not client.sio.connected
