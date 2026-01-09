"""Test script to update camera from Python and observe frontend changes."""

import time

from zndraw import ZnDraw
from zndraw.geometries import Camera

# Connect to the server
vis = ZnDraw(url="http://127.0.0.1:5000")

print(f"Connected to ZnDraw server")
# print(f"Available rooms: {vis.api.get_all_room_names()}")

# List frontend sessions
sessions = vis.api.list_frontend_sessions()
print(f"\nFrontend sessions: {sessions}")

if sessions:
    session_id = sessions[0]
    print(f"Using session: {session_id}")

    # Get current camera
    cam = vis.api.get_session_camera(session_id)
    print(f"\nCurrent camera state:")
    print(f"  Position: {cam.position}")
    print(f"  Target: {cam.target}")
    print(f"  Zoom: {cam.zoom}")

    # Wait a bit
    time.sleep(2)

    # Move camera to a dramatically different position
    print(f"\nUpdating camera to position (50, 50, 50)...")
    cam.position = (50.0, 50.0, 50.0)
    vis.api.set_session_camera(session_id, cam)
    print("Camera updated!")

    # Verify the update
    time.sleep(1)
    updated_cam = vis.api.get_session_camera(session_id)
    print(f"\nVerified updated camera:")
    print(f"  Position: {updated_cam.position}")
    print(f"  Target: {updated_cam.target}")
    print(f"  Zoom: {updated_cam.zoom}")
else:
    print("\nERROR: No frontend sessions found!")
    print("Make sure the browser is open and connected to the room")
