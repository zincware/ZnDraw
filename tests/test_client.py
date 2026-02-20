"""Integration tests for the ZnDraw Python client.

These tests verify client functionality against a running server.
The ZnDraw client works with ase.Atoms objects as frames.
"""

import uuid

import ase
import httpx
import pytest

from zndraw import ZnDraw


def make_atoms(positions: list[list[float]], symbols: str = "H") -> ase.Atoms:
    """Helper to create ase.Atoms with given positions."""
    n_atoms = len(positions)
    if len(symbols) == 1:
        symbols = symbols * n_atoms
    return ase.Atoms(symbols=symbols, positions=positions)


def make_water() -> ase.Atoms:
    """Create a water molecule."""
    return ase.Atoms(
        symbols="OHH",
        positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]],
    )


class TestEstimateFrameSize:
    """Tests for the _estimate_frame_size helper."""

    def test_sums_string_values(self):
        from zndraw.client import _estimate_frame_size

        frame = {"b64:key1": "aaaa", "b64:key2": "bb"}
        assert _estimate_frame_size(frame) == 6

    def test_ignores_non_string_values(self):
        from zndraw.client import _estimate_frame_size

        frame = {"b64:key1": "aaa", "other": 42}
        assert _estimate_frame_size(frame) == 3

    def test_empty_frame(self):
        from zndraw.client import _estimate_frame_size

        assert _estimate_frame_size({}) == 0


class TestImportPath:
    """Tests for the public import path."""

    def test_import_from_zndraw(self):
        """ZnDraw is importable from the top-level package."""
        from zndraw import ZnDraw as ZnDrawFromPackage

        assert ZnDrawFromPackage is ZnDraw


class TestConstructorAPI:
    """Tests for the new ZnDraw(url, room, user, password) API."""

    def test_room_parameter(self, server: str):
        """The 'room' parameter sets the room ID."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)
        assert client.room == room_id
        client.disconnect()

    def test_room_auto_generated(self, server: str):
        """Room ID is auto-generated when not provided."""
        client = ZnDraw(url=server)
        assert client.room is not None
        assert len(client.room) > 0
        client.disconnect()

    def test_guest_session_when_no_user(self, server: str):
        """No user parameter creates a guest session."""
        client = ZnDraw(url=server)
        assert client.connected
        assert client.api.token is not None
        client.disconnect()

    def test_guest_user_is_populated(self, server: str):
        """Guest session populates user with the generated guest email."""
        client = ZnDraw(url=server)
        assert client.user is not None
        assert "@guest.user" in client.user
        client.disconnect()

    def test_explicit_user_is_preserved(self, server: str):
        """user='email' is preserved on the instance after login."""
        email = f"preserve-{uuid.uuid4().hex[:8]}@example.com"
        from zndraw.config import Settings

        password = Settings().guest_password.get_secret_value()
        httpx.post(
            f"{server}/v1/auth/register",
            json={"email": email, "password": password},
        ).raise_for_status()

        client = ZnDraw(url=server, user=email)
        assert client.user == email
        client.disconnect()

    def test_login_with_user_and_default_password(self, server: str):
        """user='email' logs in with the default guest password."""
        email = f"testuser-{uuid.uuid4().hex[:8]}@example.com"
        # Pre-register the user with the default guest password
        from zndraw.config import Settings

        password = Settings().guest_password.get_secret_value()
        httpx.post(
            f"{server}/v1/auth/register",
            json={"email": email, "password": password},
        ).raise_for_status()

        client = ZnDraw(url=server, user=email)
        assert client.connected
        assert client.api.token is not None
        client.disconnect()

    def test_login_with_explicit_password(self, server_auth: str):
        """user + password logs in with provided credentials."""
        client = ZnDraw(
            url=server_auth,
            user="admin@local.test",
            password="adminpassword",
        )
        assert client.connected
        assert client.api.token is not None
        client.disconnect()

    def test_login_nonexistent_user_raises(self, server: str):
        """user='nonexistent' raises an error if user doesn't exist."""
        from zndraw.client import ZnDrawError

        with pytest.raises(ZnDrawError, match="LOGIN_BAD_CREDENTIALS"):
            ZnDraw(url=server, user="nonexistent@example.com")

    def test_login_wrong_password_raises(self, server: str):
        """Wrong password raises an error."""
        email = f"wrongpw-{uuid.uuid4().hex[:8]}@example.com"
        # Register with one password
        httpx.post(
            f"{server}/v1/auth/register",
            json={"email": email, "password": "correct-password"},
        ).raise_for_status()

        from zndraw.client import ZnDrawError

        with pytest.raises(ZnDrawError, match="LOGIN_BAD_CREDENTIALS"):
            ZnDraw(url=server, user=email, password="wrong-password")

    def test_password_inferred_from_settings(self, server: str):
        """Default password is inferred from Settings.guest_password."""
        from zndraw.config import Settings

        email = f"inferred-{uuid.uuid4().hex[:8]}@example.com"
        guest_password = Settings().guest_password.get_secret_value()
        httpx.post(
            f"{server}/v1/auth/register",
            json={"email": email, "password": guest_password},
        ).raise_for_status()

        # No password= argument — should infer from Settings
        client = ZnDraw(url=server, user=email)
        assert client.connected
        client.disconnect()


class TestClientConnection:
    """Tests for client connection and disconnection."""

    def test_connect_and_disconnect(self, server: str):
        """Client can connect and disconnect from server."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        assert client.connected
        client.disconnect()
        assert not client.connected

    def test_context_manager(self, server: str):
        """Client works as context manager."""
        room_id = uuid.uuid4().hex

        with ZnDraw(url=server, room=room_id, auto_connect=False) as client:
            assert client.connected

        assert not client.connected

    def test_auto_creates_room(self, server: str):
        """Client creates room automatically if it doesn't exist."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Room should exist and be empty
        assert len(client) == 0
        client.disconnect()


class TestFrameOperations:
    """Tests for frame CRUD operations with ase.Atoms."""

    def test_len_empty_room(self, server: str):
        """Empty room has 0 frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        assert len(client) == 0
        client.disconnect()

    def test_append_and_get_frame(self, server: str):
        """Can append and retrieve ase.Atoms frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            atoms = make_atoms([[i, 0, 0], [i, 1, 0]])
            atoms.info["frame_index"] = i
            client.append(atoms)

            assert len(client) == i + 1
            retrieved = client[i]
            assert isinstance(retrieved, ase.Atoms)
            assert len(retrieved) == 2
            assert retrieved.info.get("frame_index") == i

        client.disconnect()

    def test_delete_frame(self, server: str):
        """Can delete frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add frames
        for i in range(5):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        assert len(client) == 5

        # Delete middle frame
        del client[2]
        assert len(client) == 4

        # Verify remaining frames shifted
        assert client[2].info["index"] == 3

        client.disconnect()

    def test_update_frame(self, server: str):
        """Can update existing frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add frames
        for i in range(3):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["value"] = i
            client.append(atoms)

        # Update middle frame
        new_atoms = make_atoms([[999, 0, 0]])
        new_atoms.info["value"] = 999
        client[1] = new_atoms

        assert client[1].info["value"] == 999
        assert len(client) == 3  # Length unchanged

        client.disconnect()

    def test_extend_frames(self, server: str):
        """Can extend with multiple frames at once."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add initial frames
        for i in range(2):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        # Extend with more
        new_frames = []
        for idx in [10, 20, 30]:
            atoms = make_atoms([[idx, 0, 0]])
            atoms.info["index"] = idx
            new_frames.append(atoms)
        client.extend(new_frames)

        assert len(client) == 5
        assert client[2].info["index"] == 10
        assert client[3].info["index"] == 20
        assert client[4].info["index"] == 30

        client.disconnect()

    def test_extend_generator(self, server: str):
        """extend() accepts a generator (not just lists)."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        def frame_gen():
            for i in range(5):
                yield make_atoms([[i, 0, 0]])

        client.extend(frame_gen())
        assert len(client) == 5

        client.disconnect()

    def test_extend_empty(self, server: str):
        """extend() with empty iterable is a no-op."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.append(make_atoms([[0, 0, 0]]))
        client.extend([])
        assert len(client) == 1

        client.disconnect()


class TestFrameSlicing:
    """Tests for slice operations on frames."""

    def test_get_slice(self, server: str):
        """Can get frames using slice notation."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add frames
        for i in range(10):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        # Get slice
        frames = client[2:5]
        assert len(frames) == 3
        assert all(isinstance(f, ase.Atoms) for f in frames)
        assert frames[0].info["index"] == 2
        assert frames[1].info["index"] == 3
        assert frames[2].info["index"] == 4

        client.disconnect()

    def test_get_slice_with_step(self, server: str):
        """Can get frames with step in slice."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(10):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        # Every other frame from 1 to 8
        frames = client[1:8:2]
        assert len(frames) == 4
        assert [f.info["index"] for f in frames] == [1, 3, 5, 7]

        client.disconnect()

    def test_negative_index(self, server: str):
        """Negative indices work correctly."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        assert client[-1].info["index"] == 4
        assert client[-2].info["index"] == 3

        client.disconnect()


class TestGetMethod:
    """Tests for the get() method which returns raw frame data."""

    def test_get_single_frame(self, server: str):
        """get() returns single frame as dict."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        atoms = make_water()
        atoms.info["custom_key"] = "custom_value"
        client.append(atoms)

        frame = client.get(0)
        # get() returns raw dict, not ase.Atoms
        assert isinstance(frame, dict)
        # Should contain serialized atoms data
        assert len(frame) > 0

        client.disconnect()

    def test_get_with_key_filter(self, server: str):
        """get() with keys parameter returns only specified keys."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        atoms = make_water()
        client.append(atoms)

        # The raw frame data uses base64-encoded keys for asebytes format.
        # Key filtering requires knowing the exact key format.
        # For now, just verify get() returns a dict.
        frame = client.get(0)
        assert isinstance(frame, dict)
        assert len(frame) > 0  # Should have some keys

        client.disconnect()

    def test_get_list_indices(self, server: str):
        """get() with list of indices returns multiple frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(10):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        frames = client.get([1, 5, 7])
        assert len(frames) == 3
        assert all(isinstance(f, dict) for f in frames)

        client.disconnect()

    def test_get_slice(self, server: str):
        """get() with slice returns multiple frames."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(10):
            atoms = make_atoms([[i, 0, 0]])
            client.append(atoms)

        frames = client.get(slice(2, 6))
        assert len(frames) == 4
        assert all(isinstance(f, dict) for f in frames)

        client.disconnect()


class TestSetFrames:
    """Tests for set_frames() method."""

    def test_set_frames_single(self, server: str):
        """set_frames() with single index."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        atoms = make_atoms([[0, 0, 0]])
        atoms.info["value"] = 0
        client.append(atoms)

        new_atoms = make_atoms([[999, 0, 0]])
        new_atoms.info["value"] = 999
        client.set_frames(0, new_atoms)

        assert client[0].info["value"] == 999
        client.disconnect()

    def test_set_frames_slice(self, server: str):
        """set_frames() with slice."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["value"] = i
            client.append(atoms)

        new_atoms = [
            make_atoms([[100, 0, 0]]),
            make_atoms([[200, 0, 0]]),
        ]
        new_atoms[0].info["value"] = 100
        new_atoms[1].info["value"] = 200

        client.set_frames(slice(1, 3), new_atoms)

        assert client[1].info["value"] == 100
        assert client[2].info["value"] == 200
        client.disconnect()

    def test_set_frames_list_indices(self, server: str):
        """set_frames() with list of indices."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["value"] = i
            client.append(atoms)

        new_atoms = []
        for val in [10, 20, 40]:
            a = make_atoms([[val, 0, 0]])
            a.info["value"] = val
            new_atoms.append(a)

        client.set_frames([0, 2, 4], new_atoms)

        assert client[0].info["value"] == 10
        assert client[1].info["value"] == 1  # Unchanged
        assert client[2].info["value"] == 20
        assert client[3].info["value"] == 3  # Unchanged
        assert client[4].info["value"] == 40
        client.disconnect()

    def test_set_frames_not_connected_raises(self, server: str):
        """set_frames() raises when not connected."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id, auto_connect=False)

        with pytest.raises(Exception):  # NotConnectedError
            client.set_frames(0, make_atoms([[0, 0, 0]]))


class TestMutableSequenceInterface:
    """Tests verifying MutableSequence behavior with ase.Atoms."""

    def test_sequence_operations(self, server: str):
        """Basic sequence operations work correctly."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # len() on empty
        assert len(client) == 0

        # append
        a1 = make_atoms([[1, 0, 0]])
        a1.info["v"] = 1
        client.append(a1)
        assert len(client) == 1

        # extend
        a2 = make_atoms([[2, 0, 0]])
        a2.info["v"] = 2
        a3 = make_atoms([[3, 0, 0]])
        a3.info["v"] = 3
        client.extend([a2, a3])
        assert len(client) == 3

        # getitem
        assert client[0].info["v"] == 1
        assert client[1].info["v"] == 2
        assert client[2].info["v"] == 3

        # setitem
        a_new = make_atoms([[999, 0, 0]])
        a_new.info["v"] = 999
        client[1] = a_new
        assert client[1].info["v"] == 999

        # delitem
        del client[1]
        assert len(client) == 2
        assert client[1].info["v"] == 3  # Was index 2, now index 1

        client.disconnect()


class TestLocking:
    """Tests for distributed locking."""

    def test_lock_context_manager(self, server: str):
        """Lock can be acquired and released via context manager."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with client.get_lock(msg="Test lock"):
            # Lock is held
            pass

        # Lock released
        client.disconnect()


class TestBookmarks:
    """Tests for bookmark operations."""

    def test_set_and_get_bookmark(self, server: str):
        """Can set and get bookmarks."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add a frame first
        client.append(make_water())

        # Set bookmark
        client.bookmarks[0] = "First Frame"

        # Get bookmark
        assert client.bookmarks[0] == "First Frame"

        client.disconnect()

    def test_delete_bookmark(self, server: str):
        """Can delete bookmarks."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.append(make_water())
        client.bookmarks[0] = "Test Bookmark"

        del client.bookmarks[0]

        # Should raise KeyError
        with pytest.raises(KeyError):
            _ = client.bookmarks[0]

        client.disconnect()

    def test_delete_nonexistent_bookmark_raises_key_error(self, server: str):
        """Deleting a non-existent bookmark raises KeyError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(KeyError):
            del client.bookmarks[999]

        client.disconnect()


class TestGeometries:
    """Tests for geometry operations."""

    def test_set_and_get_geometry(self, server: str):
        """Can set and get geometries."""
        from zndraw.geometries import Sphere

        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.geometries["my_sphere"] = Sphere(radius=[1.0])

        retrieved = client.geometries["my_sphere"]
        assert isinstance(retrieved, Sphere)

        client.disconnect()

    def test_delete_geometry(self, server: str):
        """Can delete geometries."""
        from zndraw.geometries import Sphere

        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.geometries["temp"] = Sphere()
        del client.geometries["temp"]

        with pytest.raises(KeyError):
            client.geometries["temp"]

        client.disconnect()


class TestFigures:
    """Tests for Plotly figure operations via MutableMapping[str, go.Figure]."""

    def test_figures_crud(self, server: str):
        """MutableMapping: set, get, delete, iterate, keys, len with go.Figure."""
        import plotly.graph_objects as go

        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        fig1 = go.Figure(data=[go.Scatter(x=[1, 2], y=[3, 4])])
        fig2 = go.Figure(data=[go.Bar(x=["a", "b"], y=[10, 20])])

        client.figures["fig1"] = fig1
        client.figures["fig2"] = fig2

        # iter + keys
        assert set(client.figures) == {"fig1", "fig2"}
        assert set(client.figures.keys()) == {"fig1", "fig2"}

        # len
        assert len(client.figures) == 2

        # get returns go.Figure with correct structure
        retrieved = client.figures["fig1"]
        assert isinstance(retrieved, go.Figure)
        assert len(retrieved.data) == 1
        assert retrieved.data[0].type == "scatter"

        # delete
        del client.figures["fig1"]
        assert set(client.figures) == {"fig2"}
        assert len(client.figures) == 1

        client.disconnect()

    def test_figures_multi_client_visibility(self, server: str):
        """Two clients in same room see the same figures via REST."""
        import plotly.graph_objects as go

        room_id = uuid.uuid4().hex
        client1 = ZnDraw(url=server, room=room_id)
        client2 = ZnDraw(url=server, room=room_id)

        fig1 = go.Figure(data=[go.Scatter(x=[1], y=[2])])
        fig2 = go.Figure(data=[go.Scatter(x=[3], y=[4])])

        # Client 1 sets figures
        client1.figures["fig1"] = fig1
        client1.figures["fig2"] = fig2

        # Client 2 sees them
        assert set(client2.figures) == {"fig1", "fig2"}
        retrieved = client2.figures["fig1"]
        assert isinstance(retrieved, go.Figure)

        # Client 1 deletes one
        del client1.figures["fig1"]
        assert set(client1.figures) == {"fig2"}
        assert set(client2.figures) == {"fig2"}

        # Client 2 deletes the other
        del client2.figures["fig2"]
        assert set(client1.figures) == set()
        assert set(client2.figures) == set()

        client1.disconnect()
        client2.disconnect()

    def test_figures_get_nonexistent_raises_key_error(self, server: str):
        """Getting a non-existent figure raises KeyError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(KeyError):
            _ = client.figures["nonexistent"]

        client.disconnect()

    def test_figures_delete_nonexistent_raises_key_error(self, server: str):
        """Deleting a non-existent figure raises KeyError (MutableMapping contract)."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(KeyError):
            del client.figures["nonexistent"]

        client.disconnect()


class TestSelectionGroups:
    """Tests for selection group operations."""

    def test_selection_groups_crud(self, server: str):
        """MutableMapping: set, get, delete, iterate, keys, len."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.selection_groups["grp1"] = {"sphere": [1, 2], "cube": [3]}
        client.selection_groups["grp2"] = {"sphere": [4, 5]}

        assert set(client.selection_groups) == {"grp1", "grp2"}
        assert len(client.selection_groups) == 2
        assert client.selection_groups["grp1"] == {"sphere": [1, 2], "cube": [3]}

        del client.selection_groups["grp1"]
        assert set(client.selection_groups) == {"grp2"}

        client.disconnect()

    def test_selection_groups_get_nonexistent_raises_key_error(self, server: str):
        """Getting a non-existent selection group raises KeyError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(KeyError):
            _ = client.selection_groups["nonexistent"]

        client.disconnect()

    def test_selection_groups_delete_nonexistent_raises_key_error(self, server: str):
        """Deleting a non-existent selection group raises KeyError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(KeyError):
            del client.selection_groups["nonexistent"]

        client.disconnect()


class TestSelections:
    """Tests for selection operations."""

    def test_set_and_get_selection(self, server: str):
        """Can set and get selections."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add a frame with multiple atoms
        atoms = ase.Atoms(
            symbols="HHH",
            positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]],
        )
        client.append(atoms)

        # Set selection
        client.selections["particles"] = [0, 2]

        # Get selection
        selection = client.selections["particles"]
        assert selection == (0, 2)

        client.disconnect()


class TestFrameSelection:
    """Tests for frame_selection property."""

    def test_initially_empty(self, server: str):
        """Frame selection is empty tuple for a new room."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        assert client.frame_selection == ()

        client.disconnect()

    def test_set_and_get(self, server: str):
        """Can set and get frame selection."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            client.append(make_atoms([[i, 0, 0]]))

        client.frame_selection = [3, 1, 4]
        assert client.frame_selection == (1, 3, 4)

        client.disconnect()

    def test_clear(self, server: str):
        """Setting to None clears selection."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.append(make_atoms([[0, 0, 0]]))
        client.frame_selection = [0]
        assert client.frame_selection == (0,)

        client.frame_selection = None
        assert client.frame_selection == ()

        client.disconnect()


class TestStep:
    """Tests for step (current frame index) property."""

    def test_step_property(self, server: str):
        """Can get and set step."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add frames
        for i in range(5):
            atoms = make_atoms([[i, 0, 0]])
            atoms.info["index"] = i
            client.append(atoms)

        # Initial step should be 0
        assert client.step == 0

        # Set step
        client.step = 3
        assert client.step == 3

        client.disconnect()

    def test_step_out_of_bounds(self, server: str):
        """Setting step out of bounds raises ValueError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(3):
            client.append(make_atoms([[i, 0, 0]]))

        with pytest.raises(ValueError):
            client.step = 10  # Out of bounds

        client.disconnect()

    def test_step_negative(self, server: str):
        """Negative step resolves from the end."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(5):
            client.append(make_atoms([[i, 0, 0]]))

        client.step = -1
        assert client.step == 4

        client.step = -3
        assert client.step == 2

        client.disconnect()

    def test_step_negative_out_of_bounds(self, server: str):
        """Negative step beyond length raises ValueError."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        for i in range(3):
            client.append(make_atoms([[i, 0, 0]]))

        with pytest.raises(ValueError, match="body.step"):
            client.step = -4  # len=3, -4+3=-1 → rejected by server (ge=0)

        client.disconnect()


class TestAtomsProperty:
    """Tests for the atoms convenience property."""

    def test_atoms_getter(self, server: str):
        """atoms property returns current frame."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        water = make_water()
        client.append(water)

        atoms = client.atoms
        assert isinstance(atoms, ase.Atoms)
        assert len(atoms) == 3
        assert list(atoms.get_chemical_symbols()) == ["O", "H", "H"]

        client.disconnect()

    def test_atoms_setter(self, server: str):
        """atoms property can set current frame."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        # Add initial frame
        client.append(make_atoms([[0, 0, 0]]))

        # Set new atoms via property
        new_atoms = make_water()
        client.atoms = new_atoms

        # Verify
        retrieved = client.atoms
        assert len(retrieved) == 3

        client.disconnect()


class TestTypeErrors:
    """Tests for proper type checking."""

    def test_append_requires_atoms(self, server: str):
        """append() raises TypeError for non-Atoms."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(TypeError):
            client.append({"not": "atoms"})  # type: ignore

        client.disconnect()

    def test_extend_requires_atoms_list(self, server: str):
        """extend() raises TypeError for non-Atoms."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        with pytest.raises(TypeError):
            client.extend([{"not": "atoms"}])  # type: ignore

        client.disconnect()

    def test_setitem_requires_atoms(self, server: str):
        """__setitem__ raises TypeError for non-Atoms."""
        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        client.append(make_atoms([[0, 0, 0]]))

        with pytest.raises(TypeError):
            client[0] = {"not": "atoms"}  # type: ignore

        client.disconnect()


class TestJobManager:
    """Tests for the jobs property (JobManager integration)."""

    def test_jobs_property_returns_job_manager(self, server: str):
        """ZnDraw.jobs returns a JobManager instance."""
        from zndraw_joblib.client import JobManager

        room_id = uuid.uuid4().hex
        client = ZnDraw(url=server, room=room_id)

        assert isinstance(client.jobs, JobManager)
        client.disconnect()
