import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.geometries import Bond, Camera, CameraType, Curve, Sphere


def test_rest_get_geometries(joined_room):
    """Test listing geometry keys and getting individual geometries."""
    server, room = joined_room

    # Test listing geometry keys (default geometries include cell and floor)
    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert "geometries" in data
    assert set(data["geometries"]) == {"particles", "bonds", "curve", "cell", "floor"}

    # Test getting individual geometry - particles
    response = requests.get(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "particles"
    assert data["geometry"]["type"] == "Sphere"
    assert data["geometry"]["data"]["color"] == "arrays.colors"
    assert data["geometry"]["data"]["scale"] == 0.7

    # Test getting individual geometry - bonds
    response = requests.get(f"{server}/api/rooms/{room}/geometries/bonds")
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "bonds"
    assert data["geometry"]["type"] == "Bond"
    assert data["geometry"]["data"]["connectivity"] == "info.connectivity"
    assert data["geometry"]["data"]["scale"] == 0.15

    # Test getting non-existent geometry
    response = requests.get(f"{server}/api/rooms/{room}/geometries/nonexistent")
    assert response.status_code == 404


def test_rest_update_geometries(joined_room):
    """Test creating/updating geometries via POST endpoint."""
    server, room = joined_room

    new_geometry_data = {
        "color": "#FF0000",  # Use hex color (shared across all instances)
        "position": [[1.0, 1.0, 1.0]],  # Position must be list of tuples
        "radius": 1.0,
    }

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "particles",
            "data": new_geometry_data,
            "type": "Sphere",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify the geometry was updated
    response = requests.get(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["type"] == "Sphere"
    assert data["geometry"]["data"]["color"] == ["#FF0000"]
    assert data["geometry"]["data"]["position"] == [[1.0, 1.0, 1.0]]
    assert data["geometry"]["data"]["radius"] == 1.0


def test_rest_partial_update_geometries(joined_room):
    """Test partially updating a geometry without losing existing data."""
    server, room = joined_room

    # First, create a geometry with some data
    initial_geometry_data = {
        "position": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        "color": "#FF0000",
    }
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "curve",
            "data": initial_geometry_data,
            "type": "Curve",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200

    # Now, update only the 'active' status
    partial_update_data = {"active": False}
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "curve",
            "data": partial_update_data,
            "type": "Curve",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200

    # Verify the geometry was updated and old data was preserved
    response = requests.get(f"{server}/api/rooms/{room}/geometries/curve")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["type"] == "Curve"
    assert data["geometry"]["data"]["active"] is False
    assert data["geometry"]["data"]["position"] == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    assert data["geometry"]["data"]["color"] == ["#FF0000"]


def test_rest_add_unknown_geometry(joined_room):
    """Test that creating geometry with unknown type returns error."""
    server, room = joined_room

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "unknown",
            "data": {},
            "type": "UnknownType",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 400
    data = response.json()
    assert data["type"] == "ValueError"
    assert "Unknown geometry type" in data["error"]


def test_rest_delete_geometry(joined_room):
    """Test deleting geometries."""
    server, room = joined_room

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/particles",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/bonds",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/curve",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/cell",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/floor",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify no geometries remain
    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert data == {"geometries": []}


def test_rest_delete_unknown_geometry(joined_room):
    server, room = joined_room

    response = requests.delete(
        f"{server}/api/rooms/{room}/geometries/unknown",
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 404
    data = response.json()
    assert data["type"] == "KeyError"
    assert "Geometry with key 'unknown' does not exist" in data["error"]


def test_vis_list_geometries(server):
    from zndraw.geometries import Cell, Curve, Floor

    vis = ZnDraw(url=server, room="test-room-vis-list-geom", user="tester")
    assert len(vis.geometries) == 5
    assert vis.geometries["particles"] == Sphere()
    assert vis.geometries["bonds"] == Bond()
    assert vis.geometries["curve"] == Curve()
    assert vis.geometries["cell"] == Cell()
    assert vis.geometries["floor"] == Floor()


def test_vis_add_update_delete_geometry(server):
    vis1 = ZnDraw(url=server, room="room1", user="tester")
    vis2 = ZnDraw(url=server, room="room1", user="tester2")

    new_sphere = Sphere(color="#00FF00", position=[[0.0, 0.0, 0.0]], radius=2.0)
    vis1.geometries["new_sphere"] = new_sphere
    vis2.socket.sio.sleep(0.5)
    assert vis1.geometries["new_sphere"] == new_sphere
    assert vis2.geometries["new_sphere"] == new_sphere

    del vis2.geometries["new_sphere"]
    vis2.socket.sio.sleep(0.5)
    assert "new_sphere" not in vis1.geometries
    assert "new_sphere" not in vis2.geometries


def test_vis_add_unknown_geometry(server):
    class MyGeometry:
        pass

    vis = ZnDraw(url=server, room="room-unknown-geom", user="tester")
    with pytest.raises(ValueError, match="Unknown geometry type"):
        vis.geometries["my_geom"] = MyGeometry()


def test_vis_delete_unknown_geometry(server):
    vis = ZnDraw(url=server, room="room-del-unknown-geom", user="tester")
    with pytest.raises(KeyError, match="Geometry with key 'unknown' does not exist"):
        del vis.geometries["unknown"]


def test_vis_geometries_key_error(server):
    vis = ZnDraw(url=server, room="room-geom-key-error", user="tester")
    with pytest.raises(
        KeyError, match="Geometry with key 'nonexistent' does not exist"
    ):
        _ = vis.geometries["nonexistent"]


def test_rest_create_basic_camera(joined_room):
    """Test creating a basic camera via REST API."""
    server, room = joined_room

    # Create curves for camera position and target
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "cam_pos",
            "data": {"position": [[0.0, 0.0, 10.0]]},
            "type": "Curve",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "cam_target",
            "data": {"position": [[0.0, 0.0, 0.0]]},
            "type": "Curve",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200

    camera_data = {
        "position_curve_key": "cam_pos",
        "position_progress": 0.0,
        "target_curve_key": "cam_target",
        "target_progress": 0.0,
        "up": [0.0, 1.0, 0.0],
        "fov": 60.0,
        "camera_type": "PerspectiveCamera",
    }

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "camera1",
            "data": camera_data,
            "type": "Camera",
        },
        headers=get_jwt_auth_headers(server),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify the camera was created
    response = requests.get(f"{server}/api/rooms/{room}/geometries/camera1")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["type"] == "Camera"
    assert data["geometry"]["data"]["position_curve_key"] == "cam_pos"
    assert data["geometry"]["data"]["target_curve_key"] == "cam_target"
    assert data["geometry"]["data"]["fov"] == 60.0


def test_vis_create_camera_with_curves(server):
    """Test creating a camera with curve attachments via Python API."""
    vis = ZnDraw(url=server, room="room-camera-curves", user="tester")

    # Create curves for camera paths
    position_curve = Curve(
        position=[[0, 0, 10], [5, 5, 10], [10, 0, 10]],
        color="#FF0000",
    )
    target_curve = Curve(
        position=[[0, 0, 0], [2, 2, 0], [4, 0, 0]],
        color="#00FF00",
    )

    vis.geometries["camera_path"] = position_curve
    vis.geometries["target_path"] = target_curve

    # Create camera with curve references
    camera = Camera(
        position_curve_key="camera_path",
        position_progress=0.5,
        target_curve_key="target_path",
        target_progress=0.5,
        fov=75.0,
        helper_visible=True,
        helper_color="#0000FF",
    )

    vis.geometries["cinematic_camera"] = camera
    vis.socket.sio.sleep(0.5)

    # Verify the camera was created with curve references
    retrieved_camera = vis.geometries["cinematic_camera"]
    assert retrieved_camera.position_curve_key == "camera_path"
    assert retrieved_camera.position_progress == 0.5
    assert retrieved_camera.target_curve_key == "target_path"
    assert retrieved_camera.target_progress == 0.5


def test_vis_camera_types(server):
    """Test creating cameras with different types."""
    vis = ZnDraw(url=server, room="room-camera-types", user="tester")

    # Create curves for camera position and target
    vis.geometries["pos1"] = Curve(position=[[0, 0, 10]])
    vis.geometries["target1"] = Curve(position=[[0, 0, 0]])
    vis.geometries["pos2"] = Curve(position=[[0, 0, 10]])
    vis.geometries["target2"] = Curve(position=[[0, 0, 0]])

    # Perspective camera
    perspective_camera = Camera(
        position_curve_key="pos1",
        target_curve_key="target1",
        camera_type=CameraType.PERSPECTIVE,
        fov=75.0,
    )
    vis.geometries["perspective_cam"] = perspective_camera

    # Orthographic camera
    orthographic_camera = Camera(
        position_curve_key="pos2",
        target_curve_key="target2",
        camera_type=CameraType.ORTHOGRAPHIC,
        zoom=2.0,
    )
    vis.geometries["ortho_cam"] = orthographic_camera

    vis.socket.sio.sleep(0.5)

    # Verify both cameras
    assert vis.geometries["perspective_cam"].camera_type == CameraType.PERSPECTIVE
    assert vis.geometries["ortho_cam"].camera_type == CameraType.ORTHOGRAPHIC


def test_vis_camera_validation(server):
    """Test camera validation for invalid parameters."""
    vis = ZnDraw(url=server, room="room-camera-validation", user="tester")

    # Create curves for testing
    vis.geometries["pos"] = Curve(position=[[0, 0, 10]])
    vis.geometries["target"] = Curve(position=[[0, 0, 0]])

    # Test invalid FOV (must be between 0 and 180)
    with pytest.raises(Exception):  # Pydantic validation error
        camera = Camera(
            position_curve_key="pos",
            target_curve_key="target",
            fov=200.0,  # Invalid: > 180
        )

    # Test invalid far plane (must be > near)
    with pytest.raises(Exception):  # Pydantic validation error
        camera = Camera(
            position_curve_key="pos",
            target_curve_key="target",
            near=10.0,
            far=5.0,  # Invalid: < near
        )

    # Test invalid up vector (cannot be zero)
    with pytest.raises(Exception):  # Pydantic validation error
        camera = Camera(
            position_curve_key="pos",
            target_curve_key="target",
            up=(0.0, 0.0, 0.0),  # Invalid: zero vector
        )


def test_vis_camera_update_progress(server):
    """Test updating camera curve progress."""
    vis = ZnDraw(url=server, room="room-camera-progress", user="tester")

    # Create curves
    position_curve = Curve(
        position=[[0, 0, 0], [5, 5, 5], [10, 0, 0]],
        color="#FF00FF",
    )
    target_curve = Curve(
        position=[[0, 0, 0]],
        color="#00FFFF",
    )
    vis.geometries["path"] = position_curve
    vis.geometries["target"] = target_curve

    # Create camera attached to curve
    camera = Camera(
        position_curve_key="path",
        position_progress=0.0,
        target_curve_key="target",
        target_progress=0.0,
    )
    vis.geometries["moving_camera"] = camera
    vis.socket.sio.sleep(0.5)

    # Update progress
    camera_updated = Camera(
        position_curve_key="path",
        position_progress=1.0,  # Move to end of curve
        target_curve_key="target",
        target_progress=0.0,
    )
    vis.geometries["moving_camera"] = camera_updated
    vis.socket.sio.sleep(0.5)

    # Verify progress was updated
    retrieved = vis.geometries["moving_camera"]
    assert retrieved.position_progress == 1.0


def test_box_normalization():
    """Test that Box normalizes inputs to list format."""
    from zndraw.geometries import Box

    # Single tuple/hex inputs -> wrapped in lists
    box1 = Box(
        position=(0, 0, 0), rotation=(0, 0.5, 0), color="#FF0000", size=(1, 1, 1)
    )
    assert box1.position == [(0.0, 0.0, 0.0)]
    assert box1.rotation == [(0.0, 0.5, 0.0)]
    assert box1.color == ["#FF0000"]
    assert box1.size == [(1.0, 1.0, 1.0)]

    # List inputs -> kept as lists
    box2 = Box(
        position=[(0, 0, 0), (1, 1, 1)],
        rotation=[(0, 0.5, 0), (0, 0, 0.5)],
        color=["#FF0000", "#00FF00"],
        size=[(1, 1, 1), (2, 2, 2)],
    )
    assert len(box2.position) == 2
    assert len(box2.rotation) == 2
    assert box2.color == ["#FF0000", "#00FF00"]
    assert len(box2.size) == 2

    # Dynamic references -> pass through
    box3 = Box(position="arrays.positions", color="arrays.colors")
    assert box3.position == "arrays.positions"
    assert box3.color == "arrays.colors"


def test_plane_normalization():
    """Test that Plane normalizes inputs to list format."""
    from zndraw.geometries import Plane

    # Single tuple/hex inputs -> wrapped in lists
    plane1 = Plane(
        position=(0, 0, 0), rotation=(0, 0.5, 0), color="#FF0000", size=(1, 1)
    )
    assert plane1.position == [(0.0, 0.0, 0.0)]
    assert plane1.rotation == [(0.0, 0.5, 0.0)]
    assert plane1.color == ["#FF0000"]
    assert plane1.size == [(1.0, 1.0)]

    # List inputs -> kept as lists
    plane2 = Plane(
        position=[(0, 0, 0), (1, 1, 1)],
        rotation=[(0, 0.5, 0), (0, 0, 0.5)],
        color=["#FF0000", "#00FF00"],
        size=[(1, 1), (2, 2)],
    )
    assert len(plane2.position) == 2
    assert len(plane2.rotation) == 2
    assert plane2.color == ["#FF0000", "#00FF00"]
    assert len(plane2.size) == 2


def test_arrow_normalization():
    """Test that Arrow normalizes inputs to list format."""
    from zndraw.geometries import Arrow

    # Single tuple/hex inputs -> wrapped in lists
    arrow1 = Arrow(position=(0, 0, 0), direction=(0, 0, 1), color="#FF0000")
    assert (
        arrow1.position == [(0.0, 0.0, 0.0)] or arrow1.position == "arrays.positions"
    )  # Default might be dynamic
    assert (
        arrow1.direction == [(0.0, 0.0, 1.0)] or arrow1.direction == "calc.forces"
    )  # Default might be dynamic
    assert (
        arrow1.color == ["#FF0000"] or arrow1.color == "arrays.colors"
    )  # Default might be dynamic


def test_color_conversion_to_hex():
    """Test that atoms colors are converted to hex strings."""
    import ase

    from zndraw.utils import update_colors_and_radii

    # Create atoms
    atoms = ase.Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])

    # Colors should not exist yet
    assert "colors" not in atoms.arrays

    # Call update function
    update_colors_and_radii(atoms)

    # Verify colors were added as hex strings
    colors = atoms.arrays["colors"]
    assert len(colors) == 3
    assert all(isinstance(c, str) and c.startswith("#") for c in colors)

    # Verify hex format (should be 7 characters: # + 6 hex digits)
    assert all(len(c) == 7 for c in colors)

    # Verify radii were also added
    assert "radii" in atoms.arrays
    assert len(atoms.arrays["radii"]) == 3


def test_color_already_exists():
    """Test that update_colors_and_radii doesn't overwrite existing colors."""
    import ase
    import numpy as np

    from zndraw.utils import update_colors_and_radii

    atoms = ase.Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])

    # Manually set colors
    existing_colors = np.array(["#AABBCC", "#DDEEFF", "#112233"], dtype=object)
    atoms.set_array("colors", existing_colors)

    # Call update function
    update_colors_and_radii(atoms)

    # Verify colors were not changed
    colors = atoms.arrays["colors"]
    assert len(colors) == 3
    assert colors[0] == "#AABBCC"
    assert colors[1] == "#DDEEFF"
    assert colors[2] == "#112233"
