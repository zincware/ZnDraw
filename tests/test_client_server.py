import uuid

import numpy as np
import pytest

from zndraw import ZnDraw


def test_connection(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    assert client.socket.connected
    client.disconnect()
    assert not client.socket.connected


def test_len_frames_empty(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    assert len(client) == 0
    assert client.socket.connected


def test_append_and_get_frame(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    for i in range(10):
        data = {
            "index": np.array([i]),
            "points": np.random.rand(5, 3),
            "colors": np.random.randint(0, 255, size=(5, 4), dtype=np.uint8),
        }
        with client.get_lock():
            client._append_frame(data)
        assert len(client) == i + 1
        frame = client.get(i)
        assert np.array_equal(frame["index"], data["index"])
        assert np.array_equal(frame["points"], data["points"])
        assert np.array_equal(frame["colors"], data["colors"])

    client.disconnect()
    assert not client.socket.connected


def test_delete_frame(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    for i in range(5):
        with client.get_lock():
            client._append_frame({"index": np.array([i])})
    assert len(client) == 5
    del client[0]
    assert len(client) == 4
    frame = client.get(0)
    assert np.array_equal(frame["index"], np.array([1]))
    del client[1]
    assert len(client) == 3
    frame = client.get(1)
    assert np.array_equal(frame["index"], np.array([3]))
    client.disconnect()
    assert not client.socket.connected


def test_replace_frame(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    for i in range(10):
        with client.get_lock():
            client._append_frame({"index": np.array([i])})
    assert len(client) == 10

    # Replace frame at index 5
    new_data = {
        "index": np.array([999]),
    }
    with client.get_lock():
        client._replace_frame(5, new_data)
    frame = client.get(5)
    assert np.array_equal(frame["index"], new_data["index"])
    assert len(client) == 10

    client.disconnect()
    assert not client.socket.connected


def test_extend_frames(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Start with a few frames
    for i in range(3):
        with client.get_lock():
            client._append_frame({"index": np.array([i]), "points": np.random.rand(5, 3)})
    assert len(client) == 3

    # Extend with multiple frames at once
    extend_data = [
        {"index": np.array([10]), "points": np.random.rand(5, 3)},
        {"index": np.array([20]), "points": np.random.rand(5, 3)},
        {"index": np.array([30]), "points": np.random.rand(5, 3)},
    ]

    with client.get_lock():
        new_indices = client._extend_frames(extend_data)
    assert len(client) == 6
    assert new_indices == [3, 4, 5]  # Should be appended at these logical positions

    # Verify the extended frames
    for i, expected_idx in enumerate([10, 20, 30]):
        frame = client.get(3 + i)
        assert np.array_equal(frame["index"], np.array([expected_idx]))
        assert frame["points"].shape == (5, 3)

    client.disconnect()
    assert not client.socket.connected


def test_get_frames_with_indices(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Add some test frames
    test_data = []
    for i in range(10):
        data = {"index": np.array([i * 10]), "points": np.random.rand(3, 3)}
        with client.get_lock():
            client._append_frame(data)
        test_data.append(data)

    # Test fetching specific frames by indices
    indices = [0, 2, 5, 8]
    frames = client.get(indices)

    assert len(frames) == len(indices)
    for i, frame_idx in enumerate(indices):
        assert np.array_equal(frames[i]["index"], test_data[frame_idx]["index"])
        assert np.array_equal(frames[i]["points"], test_data[frame_idx]["points"])

    client.disconnect()
    assert not client.socket.connected


def test_get_frames_with_slice(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Add some test frames
    test_data = []
    for i in range(20):
        data = {"index": np.array([i])}
        with client.get_lock():
            client._append_frame(data)
        test_data.append(data)

    # Test fetching frames with slice notation
    frames = client.get(slice(5, 15, 2))  # Every 2nd frame from 5 to 15
    expected_indices = [5, 7, 9, 11, 13]

    assert len(frames) == len(expected_indices)
    for i, expected_idx in enumerate(expected_indices):
        assert np.array_equal(frames[i]["index"], test_data[expected_idx]["index"])

    # Test slice with just start
    frames = client.get(slice(10, None))  # From 10 to end
    assert len(frames) == 10  # frames 10-19
    for i, frame in enumerate(frames):
        assert np.array_equal(frame["index"], test_data[10 + i]["index"])

    # Test slice with step only
    frames = client.get(slice(None, None, 3))  # Every 3rd frame
    expected_count = (20 + 2) // 3  # Ceiling division
    assert len(frames) == expected_count

    # Test slice(None, None, None) - should get all frames
    frames = client.get(slice(None, None, None))
    assert len(frames) == 20  # All frames
    for i, frame in enumerate(frames):
        assert np.array_equal(frame["index"], test_data[i]["index"])

    client.disconnect()
    assert not client.socket.connected


def test_get_frames_empty_result(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Add a few frames
    with client.get_lock():
        for i in range(5):
            client._append_frame({"index": np.array([i])})

    # Test slice that results in no frames
    frames = client.get(slice(10, 20))  # Beyond available frames
    assert len(frames) == 0

    # Test empty indices list
    frames = client.get([])
    assert len(frames) == 0

    client.disconnect()
    assert not client.socket.connected


def test_mutable_sequence_interface(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Test len()
    assert len(client) == 0

    # Test append()
    frame1 = {"index": np.array([10])}
    frame2 = {"index": np.array([20])}
    with client.get_lock():
        client._append_frame(frame1)
    assert len(client) == 1
    with client.get_lock():
        client._append_frame(frame2)
    assert len(client) == 2

    # Test __getitem__ with positive index
    retrieved = client.get(0)
    assert np.array_equal(retrieved["index"], frame1["index"])

    # Test __getitem__ with negative index
    retrieved = client.get(-1)
    assert np.array_equal(retrieved["index"], frame2["index"])

    # Test __getitem__ with slice
    frames = client.get(slice(0, 2))
    assert len(frames) == 2
    assert np.array_equal(frames[0]["index"], frame1["index"])
    assert np.array_equal(frames[1]["index"], frame2["index"])

    # Test __setitem__ (replace)
    new_frame = {"index": np.array([999])}
    with client.get_lock():
        client._replace_frame(1, new_frame)
    retrieved = client.get(1)
    assert np.array_equal(retrieved["index"], new_frame["index"])

    # Test extend()
    extend_data = [{"index": np.array([100])}, {"index": np.array([200])}]
    with client.get_lock():
        client._extend_frames(extend_data)
    assert len(client) == 4
    assert np.array_equal(client.get(2)["index"], extend_data[0]["index"])
    assert np.array_equal(client.get(3)["index"], extend_data[1]["index"])

    # Test __delitem__
    client.get(1)  # Should be the replaced frame with index 999
    del client[1]
    assert len(client) == 3
    # After deletion, what was at index 2 should now be at index 1
    assert np.array_equal(client.get(1)["index"], extend_data[0]["index"])

    client.disconnect()
    assert not client.socket.connected


def test_insert_frame_functionality(server, s22):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.connect()

    # Start with some initial frames
    initial_frames = [
        {"index": np.array([0]), "value": np.array([100])},
        {"index": np.array([1]), "value": np.array([200])},
        {"index": np.array([2]), "value": np.array([300])},
    ]
    with client.get_lock():
        for frame in initial_frames:
            client._append_frame(frame)

    assert len(client) == 3

    # Test insert at beginning
    insert_frame_0 = {"index": np.array([999]), "value": np.array([999])}
    with client.get_lock():
        client._insert_frame(0, insert_frame_0)
    assert len(client) == 4

    # Verify the insertion shifted everything
    frame = client.get(0)
    assert np.array_equal(frame["index"], insert_frame_0["index"])
    frame = client.get(1)  # Should be the original frame 0
    assert np.array_equal(frame["index"], initial_frames[0]["index"])

    # Test insert in middle
    insert_frame_2 = {"index": np.array([888]), "value": np.array([888])}
    with client.get_lock():
        client._insert_frame(2, insert_frame_2)
    assert len(client) == 5

    # Verify the insertion
    frame = client.get(2)
    assert np.array_equal(frame["index"], insert_frame_2["index"])
    # Original frame 1 should now be at position 3
    frame = client.get(3)
    assert np.array_equal(frame["index"], initial_frames[1]["index"])

    # Test insert at end (equivalent to append)
    insert_frame_end = {"index": np.array([777]), "value": np.array([777])}
    with client.get_lock():
        client._insert_frame(len(client), insert_frame_end)
    assert len(client) == 6

    # Should be at the last position
    frame = client.get(-1)
    assert np.array_equal(frame["index"], insert_frame_end["index"])

    client.disconnect()
    assert not client.socket.sio.connected


def test_mutable_sequence_insert(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.connect()

    # Add initial frames using MutableSequence interface
    frames = [
        {"data": np.array([1, 2, 3])},
        {"data": np.array([4, 5, 6])},
        {"data": np.array([7, 8, 9])},
    ]
    with client.get_lock():
        for frame in frames:
            client._append_frame(frame)

    assert len(client) == 3

    # Test MutableSequence insert at position 1
    insert_data = {"data": np.array([99, 99, 99])}
    with client.get_lock():
        client._insert_frame(1, insert_data)
    assert len(client) == 4

    # Verify the order
    assert np.array_equal(client.get(0)["data"], frames[0]["data"])  # Original frame 0
    assert np.array_equal(client.get(1)["data"], insert_data["data"])  # Inserted frame
    assert np.array_equal(
        client.get(2)["data"], frames[1]["data"]
    )  # Shifted from pos 1 to 2
    assert np.array_equal(
        client.get(3)["data"], frames[2]["data"]
    )  # Shifted from pos 2 to 3

    client.disconnect()
    assert not client.socket.sio.connected


def test_slice_assignment_vs_python_list(server):
    """Test that slice assignment behavior matches Python list behavior exactly."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.connect()

    # Test cases covering various slice assignment scenarios
    test_cases = [
        # (initial_data, slice_assignment, values_to_assign)
        ([1, 2, 3, 4, 5], slice(1, 3), [10, 20]),  # data[1:3] = [10, 20]
        ([1, 2, 3], slice(3, None), [4, 5, 6]),  # data[3:] = [4, 5, 6] (extend)
        ([1, 2, 3, 4, 5], slice(None, 2), [10, 20]),  # data[:2] = [10, 20]
        ([1, 2, 3, 4, 5], slice(1, 4), [10]),  # data[1:4] = [10] (shrink)
        ([1, 2, 3], slice(1, 1), [10, 20]),  # data[1:1] = [10, 20] (insert)
        ([1, 2, 3, 4, 5], slice(2, 4), []),  # data[2:4] = [] (delete range)
    ]

    for i, (initial, slice_obj, values) in enumerate(test_cases):
        # Test with Python list
        python_list = initial.copy()
        python_list[slice_obj] = values

        # Test with client
        # Convert integers to numpy arrays for client
        initial_frames = [{"index": np.array([x])} for x in initial]
        value_frames = [{"index": np.array([x])} for x in values]

        # Clear client and set up initial data
        while len(client) > 0:
            del client[0]
        with client.get_lock():
            for frame in initial_frames:
                client._append_frame(frame)

        # Perform slice assignment
        client.set_frames(slice_obj, value_frames)

        # Compare results
        assert len(client) == len(python_list), (
            f"Test case {i}: Length mismatch. Client: {len(client)}, Python: {len(python_list)}"
        )

        for j in range(len(python_list)):
            client_value = client.get(j)["index"][0]
            python_value = python_list[j]
            assert client_value == python_value, (
                f"Test case {i}, index {j}: Client: {client_value}, Python: {python_value}"
            )

    client.disconnect()
    assert not client.socket.sio.connected


def test_extended_slice_assignment_vs_python_list(server):
    """Test extended slice assignment behavior matches Python list behavior."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.connect()

    # Test cases for extended slices (step != 1)
    test_cases = [
        # (initial_data, slice_assignment, values_to_assign)
        (
            [1, 2, 3, 4, 5, 6],
            slice(None, None, 2),
            [10, 20, 30],
        ),  # data[::2] = [10, 20, 30]
        (
            [1, 2, 3, 4, 5, 6],
            slice(1, None, 2),
            [10, 20, 30],
        ),  # data[1::2] = [10, 20, 30]
        (
            [1, 2, 3, 4, 5, 6, 7, 8],
            slice(2, 7, 2),
            [10, 20, 30],
        ),  # data[2:7:2] = [10, 20, 30]
    ]

    for i, (initial, slice_obj, values) in enumerate(test_cases):
        # Test with Python list
        python_list = initial.copy()
        python_list[slice_obj] = values

        # Test with client
        initial_frames = [{"index": np.array([x])} for x in initial]
        value_frames = [{"index": np.array([x])} for x in values]

        # Clear client and set up initial data
        while len(client) > 0:
            del client[0]

        with client.get_lock():
            for frame in initial_frames:
                client._append_frame(frame)

        # Perform extended slice assignment
        client.set_frames(slice_obj, value_frames)

        # Compare results
        assert len(client) == len(python_list), (
            f"Extended test case {i}: Length mismatch. Client: {len(client)}, Python: {len(python_list)}"
        )

        for j in range(len(python_list)):
            client_value = client.get(j)["index"][0]
            python_value = python_list[j]
            assert client_value == python_value, (
                f"Extended test case {i}, index {j}: Client: {client_value}, Python: {python_value}"
            )

    client.disconnect()
    assert not client.socket.sio.connected


def test_slice_assignment_error_conditions(server):
    """Test that slice assignment error conditions match Python behavior."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Set up initial data
    initial_frames = [{"index": np.array([i])} for i in range(1, 6)]  # [1, 2, 3, 4, 5]
    
    with client.get_lock():
        for frame in initial_frames:
            client._append_frame(frame)

    # Test extended slice with wrong number of values (should raise ValueError)
    value_frames = [{"index": np.array([10])}]  # Only 1 value for 3 positions

    # Test that Python raises ValueError for this case
    python_list = [1, 2, 3, 4, 5]
    with pytest.raises(
        ValueError,
        match="attempt to assign sequence of size 1 to extended slice of size 3",
    ):
        python_list[::2] = [10]  # Trying to assign 1 value to 3 positions

    # Test that client also raises ValueError with same message
    with pytest.raises(
        ValueError,
        match="attempt to assign sequence of size 1 to extended slice of size 3",
    ):
        client.set_frames(slice(0, None, 2), value_frames)

    client.disconnect()
    assert not client.socket.connected


def test_slice_assignment_connection_error(server):
    """Test that slice assignment fails when not connected."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.disconnect()
    assert not client.socket.connected
    # Deliberately not connecting

    value_frames = [{"index": np.array([10])}]

    # Test that both simple and extended slice assignment fail when not connected
    with pytest.raises(RuntimeError, match="Client is not connected"):
        client.set_frames(slice(1, 3), value_frames)

    with pytest.raises(RuntimeError, match="Client is not connected"):
        client.set_frames(slice(0, None, 2), value_frames)


def test_slice_assignment_edge_cases(server):
    """Test additional edge cases for slice assignment."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    client.connect()

    # Test cases for edge conditions
    edge_cases = [
        # Empty list operations
        ([], slice(0, 0), [1, 2, 3]),  # empty[:] = [1, 2, 3]
        ([], slice(None), [1, 2]),  # empty[:] = [1, 2]
        # Negative indices in slices
        ([1, 2, 3, 4, 5], slice(-3, -1), [10, 20]),  # data[-3:-1] = [10, 20]
        ([1, 2, 3, 4, 5], slice(-2, None), [10]),  # data[-2:] = [10]
        # Step = -1 (reverse slice assignment)
        (
            [1, 2, 3, 4, 5],
            slice(None, None, -1),
            [10, 20, 30, 40, 50],
        ),  # data[::-1] = [...]
        # Out of bounds operations
        ([1, 2, 3], slice(5, 10), [10, 20]),  # Beyond end of list
        ([1, 2, 3], slice(-10, 2), [10, 20]),  # Before start of list
    ]

    for i, (initial, slice_obj, values) in enumerate(edge_cases):
        # Test with Python list
        python_list = initial.copy()
        python_list[slice_obj] = values

        # Test with client
        initial_frames = [{"index": np.array([x])} for x in initial]
        value_frames = [{"index": np.array([x])} for x in values]

        # Clear client and set up initial data
        while len(client) > 0:
            del client[0]
        with client.get_lock():
            for frame in initial_frames:
                client._append_frame(frame)

        # Perform slice assignment
        client.set_frames(slice_obj, value_frames)

        # Compare results
        assert len(client) == len(python_list), (
            f"Edge case {i}: Length mismatch. Client: {len(client)}, Python: {len(python_list)}"
        )

        for j in range(len(python_list)):
            client_value = client.get(j)["index"][0]
            python_value = python_list[j]
            assert client_value == python_value, (
                f"Edge case {i}, index {j}: Client: {client_value}, Python: {python_value}"
            )

    client.disconnect()
    assert not client.socket.connected


def test_slice_assignment_single_value(server):
    """Test slice assignment with single value (not in list)."""
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Set up initial data
    initial_frames = [{"index": np.array([i])} for i in range(1, 4)]  # [1, 2, 3]
    with client.get_lock():
        for frame in initial_frames:
            client._append_frame(frame)

    # Test assigning single frame to slice (should be treated as [frame])
    single_frame = {"index": np.array([99])}

    # Python behavior: single value gets wrapped in list for slice assignment
    python_list = [1, 2, 3]
    python_list[1:3] = [99]  # Equivalent to our operation

    # Client operation
    client.set_frames(slice(1, 3), [single_frame])  # Should be equivalent

    # Compare results
    assert len(client) == len(python_list), (
        f"Single value: Length mismatch. Client: {len(client)}, Python: {len(python_list)}"
    )

    for j in range(len(python_list)):
        client_value = client.get(j)["index"][0]
        python_value = python_list[j]
        assert client_value == python_value, (
            f"Single value, index {j}: Client: {client_value}, Python: {python_value}"
        )

    client.disconnect()
    assert not client.socket.connected


def test_nested_dict_handling(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)

    # Append a frame with nested dictionaries
    nested_data = {
        "index": np.array([1, 2, 3]),
        "meta": {
            "author": "Test",
            "details": {"version": 1.0, "tags": ["test", "nested"]},
        },
    }
    
    with client.get_lock():
        client._append_frame(nested_data)
    assert len(client) == 1

    # Retrieve and verify the nested structure
    frame = client.get(0)
    assert np.array_equal(frame["index"], nested_data["index"])
    assert frame["meta"]["author"] == nested_data["meta"]["author"]
    assert (
        frame["meta"]["details"]["version"] == nested_data["meta"]["details"]["version"]
    )
    assert frame["meta"]["details"]["tags"] == nested_data["meta"]["details"]["tags"]

    client.disconnect()
    assert not client.socket.connected


def test_comprehensive_atom_dict(server):
    data = {
        "numbers": np.array([6, 1, 1, 1, 1]),
        "positions": np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
            ]
        ),
        "tags": np.array([1, 2, 3, 4, 5]),
        "initial_charges": np.array(
            [0.68062226, 0.77952437, 0.98666274, 0.12193437, 0.94664167]
        ),
        "momenta": np.array(
            [
                [0.11615927, 0.03696042, 0.46672605],
                [0.15866194, 0.11610305, 0.83340719],
                [0.59826544, 0.67848578, 0.74957987],
                [0.14949151, 0.68352953, 0.92811033],
                [0.40530313, 0.61190158, 0.43241552],
            ]
        ),
        "name": np.array(["Carbon", "Hydrogen", "Hydrogen", "Hydrogen", "Hydrogen"]),
        "info": {
            "string": "Lorem Ipsum",
            "float": 3.14,
            "list": [1, 2, 3],
            "array": np.array([1.0, 2.0, 3.0]),
            "dict": {"a": 1, "b": 2},
            "d2": {"a": [1, 2], "b": [3, 4]},
            "REPEATED_KEY": "info",
        },
        "REPEATED_KEY": np.array([0, 1, 2, 3, 4]),
        "cell": np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
        "pbc": np.array([True, False, True]),
        "celldisp": np.array([0.1, 0.2, 0.3]),
        "<SinglePointCalculator>": {
            "energy": 1.0,
            "forces": np.array(
                [
                    [0.22055347, 0.10920619, 0.27019844],
                    [0.72417325, 0.71259895, 0.77200743],
                    [0.07070797, 0.14246075, 0.31810531],
                    [0.34846727, 0.33212284, 0.09877173],
                    [0.46943827, 0.29961628, 0.43061745],
                ]
            ),
            "string-calc": "this is a string from calc",
        },
    }

    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    with client.get_lock():
        client._append_frame(data)
    assert len(client) == 1
    frame = client.get(0)
    assert np.array_equal(frame["numbers"], data["numbers"])
    assert np.array_equal(frame["positions"], data["positions"])
    assert np.array_equal(frame["tags"], data["tags"])
    assert np.allclose(frame["initial_charges"], data["initial_charges"])
    assert np.array_equal(frame["momenta"], data["momenta"])
    assert np.array_equal(frame["name"], data["name"])
    assert frame["info"]["string"] == data["info"]["string"]
    assert frame["info"]["float"] == data["info"]["float"]
    assert frame["info"]["list"] == data["info"]["list"]
    assert np.array_equal(frame["info"]["array"], data["info"]["array"])
    assert frame["info"]["dict"] == data["info"]["dict"]
    assert frame["info"]["d2"] == data["info"]["d2"]
    assert frame["info"]["REPEATED_KEY"] == data["info"]["REPEATED_KEY"]
    assert np.array_equal(frame["REPEATED_KEY"], data["REPEATED_KEY"])
    assert np.array_equal(frame["cell"], data["cell"])
    assert np.array_equal(frame["pbc"], data["pbc"])
    assert np.array_equal(frame["celldisp"], data["celldisp"])
    assert (
        frame["<SinglePointCalculator>"]["energy"]
        == data["<SinglePointCalculator>"]["energy"]
    )
    assert np.array_equal(
        frame["<SinglePointCalculator>"]["forces"],
        data["<SinglePointCalculator>"]["forces"],
    )
    assert (
        frame["<SinglePointCalculator>"]["string-calc"]
        == data["<SinglePointCalculator>"]["string-calc"]
    )


def test_partial_key_retrieval(server):
    client = ZnDraw(room=uuid.uuid4().hex, url=server)
    data = {
        "index": np.array([1, 2, 3]),
        "points": np.random.rand(5, 3),
        "colors": np.random.randint(0, 255, size=(5, 4), dtype=np.uint8),
        "extra": "This is some extra info",
    }
    with client.get_lock():
        client._append_frame(data)
    assert len(client) == 1

    # Retrieve only 'points' key
    frame = client.get(0, keys=["points"])
    assert frame.keys() == {"points"}
    assert np.array_equal(frame["points"], data["points"])

    # Retrieve 'index' and 'extra' keys
    frame = client.get(0, keys=["index", "extra"])
    assert frame.keys() == {"index", "extra"}
    assert np.array_equal(frame["index"], data["index"])
    assert frame["extra"] == data["extra"]

    # do the same with get using list
    frames = client.get([0], keys=["colors"])
    assert len(frames) == 1
    assert frames[0].keys() == {"colors"}
    assert np.array_equal(frames[0]["colors"], data["colors"])

    with pytest.raises(KeyError):
        client.get(0, keys=["nonexistent_key"])

    with pytest.raises(KeyError):
        client.get([0], keys=["nonexistent_key"])

    with pytest.raises(IndexError):
        client.get(5, keys=["points"])

    with pytest.raises(IndexError):
        client.get([0, 5], keys=["points"])

    client.disconnect()
    assert not client.socket.connected


# TODO: test keys that contain "." or `<SinglePointCalculator>`
# TODO: test replacing only single keys / nested keys via a.b (think atoms.arrats["plots"])


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
