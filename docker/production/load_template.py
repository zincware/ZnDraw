#!/usr/bin/env python
"""Load template file into ZnDraw and configure it as a locked default room."""

import os
import sys

import ase.io
import requests

from zndraw import ZnDraw

server_url = os.environ.get("ZNDRAW_SERVER_URL", "http://nginx")
template_file = sys.argv[1] if len(sys.argv) > 1 else "/app/templates/ethanol.xyz"
room_name = os.environ.get("TEMPLATE_ROOM_NAME", "template")
admin_user = os.environ.get("ZNDRAW_ADMIN_USERNAME")
admin_pass = os.environ.get("ZNDRAW_ADMIN_PASSWORD")

print(f"Loading template from {template_file} to room '{room_name}'")

# Connect as admin to create and configure the template room
vis = ZnDraw(url=server_url, room=room_name, user=admin_user, password=admin_pass)
print(f"Connected to room: {room_name} as {vis.user}")

# Read and upload the template file
atoms_list = list(ase.io.iread(template_file))
print(f"Read {len(atoms_list)} frames from {template_file}")

vis.extend(atoms_list)
print(f"Uploaded {len(atoms_list)} frames")

# Use the JWT token from the API manager for authenticated requests
headers = {"Authorization": f"Bearer {vis.api.jwt_token}"}

# Lock the room
resp = requests.patch(
    f"{server_url}/api/rooms/{room_name}",
    json={"locked": True},
    headers=headers,
)
resp.raise_for_status()
print(f"Room '{room_name}' locked")

# Set as default template
resp = requests.put(
    f"{server_url}/api/rooms/default",
    json={"roomId": room_name},
    headers=headers,
)
resp.raise_for_status()
print(f"Room '{room_name}' set as default template")

print("Template loading complete!")
