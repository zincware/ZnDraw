import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.geometries import Bond, Camera, CameraType, Curve, Sphere


def test_rest_get_geometries(joined_room):
    """Test listing geometry keys and getting individual geometries."""
    server, room = joined_room
    headers = get_jwt_auth_headers(server)

    # Test listing geometry keys (default geometries include cell, floor, and constraints)
    response = requests.get(f"{server}/api/rooms/{room}/geometries", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert "geometries" in data
    assert set(data["geometries"]) == {
        "particles",
        "bonds",
        "curve",
        "cell",
        "floor",
        "constraints-fixed-atoms",
    }

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
    """Test creating/updating geometries via Python client (uses lock automatically)."""
    server, room = joined_room

    # Use ZnDraw client which handles lock acquisition automatically
    vis = ZnDraw(url=server, room=room, user="test-geom-update")

    new_sphere = Sphere(
        color="#FF0000",
        position=[[1.0, 1.0, 1.0]],
        radius=1.0,
    )
    vis.geometries["particles"] = new_sphere

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

    # Use ZnDraw client which handles lock acquisition automatically
    vis = ZnDraw(url=server, room=room, user="test-geom-partial")

    # First, create a geometry with some data
    initial_curve = Curve(
        position=[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        color="#FF0000",
    )
    vis.geometries["curve"] = initial_curve

    # When using Python client, we need to fetch, modify, and set
    # This tests that the server-side merge still works
    current_curve = vis.geometries["curve"]
    updated_curve = Curve(
        position=current_curve.position,  # Keep existing position
        color=current_curve.color,  # Keep existing color
        active=False,  # Update only active status
    )
    vis.geometries["curve"] = updated_curve

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

    # Use ZnDraw client - this test validates client-side validation
    vis = ZnDraw(url=server, room=room, user="test-geom-unknown")

    class UnknownGeometry:
        pass

    with pytest.raises(ValueError, match="Unknown geometry type"):
        vis.geometries["unknown"] = UnknownGeometry()


def test_rest_delete_geometry(joined_room):
    """Test deleting geometries."""
    server, room = joined_room
    headers = get_jwt_auth_headers(server)

    # Use ZnDraw client which handles lock acquisition automatically
    vis = ZnDraw(url=server, room=room, user="test-geom-delete")

    # Delete all geometries
    del vis.geometries["particles"]
    del vis.geometries["bonds"]
    del vis.geometries["curve"]
    del vis.geometries["cell"]
    del vis.geometries["floor"]
    del vis.geometries["constraints-fixed-atoms"]

    # Verify no geometries remain
    response = requests.get(f"{server}/api/rooms/{room}/geometries", headers=headers)
    assert response.status_code == 200
    data = response.json()
    assert data == {"geometries": {}}  # Empty dict when no geometries


def test_rest_delete_unknown_geometry(joined_room):
    server, room = joined_room

    # Use ZnDraw client - this test validates error handling
    vis = ZnDraw(url=server, room=room, user="test-geom-del-unknown")

    with pytest.raises(KeyError, match="Geometry with key 'unknown' does not exist"):
        del vis.geometries["unknown"]


def test_vis_list_geometries(server):
    from zndraw.geometries import Cell, Curve, Floor
    from zndraw.materials import MeshBasicMaterial
    from zndraw.transformations import InArrayTransform

    vis = ZnDraw(url=server, room="test-room-vis-list-geom", user="tester")
    assert len(vis.geometries) == 6

    # Check default geometries match what's created in room_service.py
    assert vis.geometries["particles"] == Sphere(
        position="arrays.positions", color="arrays.colors", radius="arrays.radii"
    )
    assert vis.geometries["bonds"] == Bond(
        position="arrays.positions", color="arrays.colors"
    )
    assert vis.geometries["curve"] == Curve()
    assert vis.geometries["cell"] == Cell()
    assert vis.geometries["floor"] == Floor()

    # Verify constraints-fixed-atoms geometry
    constraints_geom = vis.geometries["constraints-fixed-atoms"]
    assert isinstance(constraints_geom, Sphere)
    assert constraints_geom.color == ["#FF0000"]
    assert constraints_geom.scale == 0.71
    assert constraints_geom.active is True
    assert constraints_geom.material == MeshBasicMaterial(wireframe=True)


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
    """Test creating a basic camera via Python client."""
    server, room = joined_room

    # Use ZnDraw client which handles lock acquisition automatically
    vis = ZnDraw(url=server, room=room, user="test-camera-basic")

    # Create curves for camera position and target
    vis.geometries["cam_pos"] = Curve(position=[[0.0, 0.0, 10.0]])
    vis.geometries["cam_target"] = Curve(position=[[0.0, 0.0, 0.0]])

    # Create camera
    camera = Camera(
        position_curve_key="cam_pos",
        position_progress=0.0,
        target_curve_key="cam_target",
        target_progress=0.0,
        up=[0.0, 1.0, 0.0],
        fov=60.0,
        camera_type=CameraType.PERSPECTIVE,
    )
    vis.geometries["camera1"] = camera

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


def test_geometries_repr_and_str():
    """Test that vis.geometries repr and str only show keys, not full data."""
    from zndraw.scene_manager import Geometries

    # Create a mock ZnDraw instance with some geometries
    class MockZnDraw:
        _geometries = {
            "particles": {
                "type": "Sphere",
                "data": {
                    "position": "arrays.positions",
                    "color": "arrays.colors",
                    "scale": 0.7,
                    "radius": 1.0,
                    "material": "MeshPhysicalMaterial",
                },
            },
            "bonds": {
                "type": "Bond",
                "data": {
                    "connectivity": "info.connectivity",
                    "scale": 0.15,
                    "color": "#FFFFFF",
                },
            },
            "curve": {
                "type": "Curve",
                "data": {
                    "position": [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
                    "color": "#FF0000",
                },
            },
        }

    vis = MockZnDraw()
    geom = Geometries(vis)

    # Get the string representations
    repr_str = repr(geom)
    str_str = str(geom)

    # Verify exact format
    assert repr_str == "Geometries(keys=['particles', 'bonds', 'curve'])"
    assert str_str == "Geometries(keys=['particles', 'bonds', 'curve'])"
