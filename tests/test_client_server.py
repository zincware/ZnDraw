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
