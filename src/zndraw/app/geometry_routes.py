"""Geometry, figure, and selection routes.

Handles geometry objects, figures for visualization, and atom/frame selections.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio
from zndraw.auth import require_auth

from .constants import SocketEvents
from .redis_keys import RoomKeys
from .route_utils import requires_lock

log = logging.getLogger(__name__)

geometries = Blueprint("geometries", __name__)


@geometries.route("/api/rooms/<string:room_id>/geometries", methods=["POST"])
@requires_lock(target="trajectory:meta")
def create_geometry(room_id: str, session_id: str, user_id: str):
    """Create or update a geometry in the room.

    Requires trajectory:meta lock (enforced by @requires_lock decorator).

    Headers:
        X-Session-ID: Session ID from /join

    Request body:
        {
            "key": "geometry_name",
            "type": "Sphere" | "Arrow" | "Bond" | "Curve",
            "data": {...}  // geometry-specific data
        }
    """
    # session_id, user_id, lock_token are injected by @requires_lock decorator
    # No need to manually authenticate or validate lock

    data = request.get_json() or {}
    key = data.get("key")
    geometry_type = data.get("type")
    geometry_data = data.get("data")

    if not key or not geometry_type or geometry_data is None:
        return {
            "error": "'key', 'type', and 'data' are required",
            "type": "ValueError",
        }, 400

    from zndraw.geometries import geometries

    if geometry_type not in geometries:
        return {
            "error": f"Unknown geometry type '{geometry_type}'",
            "type": "ValueError",
        }, 400

    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    existing_geometry_json = r.hget(keys.geometries(), key)

    if existing_geometry_json:
        existing_geometry = json.loads(existing_geometry_json)
        # Make sure we're updating the same type of geometry
        if existing_geometry.get("type") == geometry_type:
            # Merge new data into existing data
            merged_data = existing_geometry.get("data", {}).copy()
            merged_data.update(geometry_data)
            geometry_data = merged_data
        else:
            log.warning(f"Geometry type mismatch for key '{key}'. Overwriting.")

    # Validate and apply defaults through Pydantic model
    try:
        geometry_class = geometries[geometry_type]
        validated_geometry = geometry_class(**geometry_data)
        value_to_store = json.dumps(
            {"type": geometry_type, "data": validated_geometry.model_dump()}
        )
    except Exception as e:
        return {
            "error": f"Invalid geometry data: {str(e)}",
            "type": "ValidationError",
        }, 400

    r.hset(keys.geometries(), key, value_to_store)
    socketio.emit(
        SocketEvents.INVALIDATE_GEOMETRY,
        {
            "key": key,
            "operation": "set",
        },
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/geometries/<string:key>", methods=["GET"]
)
def get_geometry(room_id: str, key: str):
    """Get a specific geometry by key.

    Returns:
        {
            "key": "geometry_name",
            "geometry": {
                "type": "Sphere",
                "data": {...}
            }
        }
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    geometry_data = r.hget(keys.geometries(), key)
    if not geometry_data:
        return {
            "error": f"Geometry with key '{key}' not found",
            "type": "KeyError",
        }, 404
    geometry = json.loads(geometry_data)
    return {"key": key, "geometry": geometry}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/geometries/<string:key>", methods=["DELETE"]
)
@requires_lock(target="trajectory:meta")
def delete_geometry(room_id: str, key: str, session_id: str, user_id: str):
    """Delete a geometry from the room.

    Requires trajectory:meta lock (enforced by @requires_lock decorator).

    Headers:
        X-Session-ID: Session ID from /join
    """
    # session_id and user_id are injected by @requires_lock decorator

    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    response = r.hdel(keys.geometries(), key)
    if response == 0:
        return {
            "error": f"Geometry with key '{key}' does not exist",
            "type": "KeyError",
        }, 404

    socketio.emit(
        SocketEvents.INVALIDATE_GEOMETRY,
        {
            "key": key,
            "operation": "delete",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@geometries.route("/api/rooms/<string:room_id>/geometries", methods=["GET"])
@require_auth
def list_geometries(room_id: str):
    """Get all geometries with their full data.

    Returns:
        {
            "geometries": {
                "particles": {"type": "Sphere", "data": {...}},
                "bonds": {"type": "Bond", "data": {...}},
                ...
            }
        }
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    geometries_raw = r.hgetall(keys.geometries())
    geometries = {k: json.loads(v) for k, v in geometries_raw.items()}
    return {"geometries": geometries}, 200


@geometries.route("/api/rooms/<string:room_id>/geometries/schemas", methods=["GET"])
def list_geometry_schemas(room_id: str):
    """Return JSON schemas for all geometry types for form generation."""
    from zndraw.geometries import geometries

    schemas = {name: model.model_json_schema() for name, model in geometries.items()}
    return {"schemas": schemas}, 200


# -------------#
### FIGURES ###
# -------------#


@geometries.route("/api/rooms/<string:room_id>/figures", methods=["POST"])
def create_figure(room_id: str):
    data = request.get_json() or {}
    key = data.get("key")
    figure = data.get("figure")

    if not key or not figure:
        return {
            "error": "Both 'key' and 'figure' are required",
            "type": "ValueError",
        }, 400

    # store in hash
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    r.hset(keys.figures(), key, json.dumps(figure))
    socketio.emit(
        SocketEvents.INVALIDATE_FIGURE,
        {
            "key": key,
            "operation": "set",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@geometries.route("/api/rooms/<string:room_id>/figures/<string:key>", methods=["GET"])
def get_figure(room_id: str, key: str):
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    figure_data = r.hget(keys.figures(), key)
    if not figure_data:
        return {"error": f"Figure with key '{key}' not found", "type": "KeyError"}, 404
    figure = json.loads(figure_data)
    return {"key": key, "figure": figure}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/figures/<string:key>", methods=["DELETE"]
)
def delete_figure(room_id: str, key: str):
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    response = r.hdel(keys.figures(), key)
    if response == 0:
        return {
            "error": f"Figure with key '{key}' does not exist",
            "type": "KeyError",
        }, 404
    socketio.emit(
        SocketEvents.INVALIDATE_FIGURE,
        {
            "key": key,
            "operation": "delete",
        },
        to=f"room:{room_id}",
    )
    return {"status": "success"}, 200


@geometries.route("/api/rooms/<string:room_id>/figures", methods=["GET"])
def list_figures(room_id: str):
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    all_keys = r.hkeys(keys.figures())
    return {"figures": list(all_keys)}, 200


# ============================================================================
# Bookmarks API Routes
# ============================================================================


@geometries.route("/api/rooms/<string:room_id>/selections", methods=["GET"])
@require_auth
def get_all_selections(room_id: str):
    """Get all current selections and groups.

    Returns:
        {
            "selections": {"particles": [1,2,3], "forces": [2,3]},
            "groups": {"group1": {"particles": [1,3], "forces": [1,3]}},
            "activeGroup": "group1" | null
        }
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Get current selections
    selections_raw = r.hgetall(keys.selections())
    selections = {k: json.loads(v) for k, v in selections_raw.items()}

    # Get selection groups
    groups_raw = r.hgetall(keys.selection_groups())
    groups = {k: json.loads(v) for k, v in groups_raw.items()}

    # Get active group
    active_group = r.get(keys.active_selection_group())

    return {
        "selections": selections,
        "groups": groups,
        "activeGroup": active_group,
    }, 200


@geometries.route("/api/rooms/<string:room_id>/frame-selection", methods=["GET"])
@require_auth
def get_frame_selection(room_id: str):
    """Get frame selection for the room.

    Returns:
        {
            "frameSelection": [0, 1, 5, 10] | null
        }
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    frame_selection = r.get(keys.frame_selection())

    return {
        "frameSelection": json.loads(frame_selection) if frame_selection else None
    }, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["GET"]
)
def get_selection(room_id: str, geometry: str):
    """Get selection for a specific geometry."""
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    selection = r.hget(keys.selections(), geometry)

    if selection is None:
        return {"selection": []}, 200

    return {"selection": json.loads(selection)}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["PUT"]
)
def update_selection(room_id: str, geometry: str):
    """Update selection for a specific geometry.

    Body: {"indices": [1, 2, 3]}
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    data = request.get_json()

    indices = data.get("indices", [])
    if not isinstance(indices, list):
        return {"error": "indices must be a list"}, 400

    if any(not isinstance(idx, int) or idx < 0 for idx in indices):
        return {"error": "All indices must be non-negative integers"}, 400

    # Store selection
    r.hset(keys.selections(), geometry, json.dumps(indices))

    # Clear active group (manual edit breaks group association)
    r.delete(keys.active_selection_group())

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"geometry": geometry},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["GET"],
)
def get_selection_group(room_id: str, group_name: str):
    """Get a specific selection group."""
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    group = r.hget(keys.selection_groups(), group_name)

    if group is None:
        return {"error": "Group not found"}, 404

    return {"group": json.loads(group)}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["PUT"],
)
def create_update_selection_group(room_id: str, group_name: str):
    """Create or update a selection group.

    Body: {"particles": [1, 3], "forces": [1, 3]}
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)
    data = request.get_json()

    # Validate data is a dict of geometry -> indices
    if not isinstance(data, dict):
        return {"error": "Group must be a dictionary"}, 400

    for geometry, indices in data.items():
        if not isinstance(indices, list):
            return {"error": f"Indices for '{geometry}' must be a list"}, 400
        if any(not isinstance(idx, int) or idx < 0 for idx in indices):
            return {"error": f"Invalid indices for '{geometry}'"}, 400

    # Store group
    r.hset(keys.selection_groups(), group_name, json.dumps(data))

    # Emit invalidation (groups list changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION_GROUPS,
        {"operation": "group_saved", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>",
    methods=["DELETE"],
)
def delete_selection_group(room_id: str, group_name: str):
    """Delete a selection group."""
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Check if group exists
    if not r.hexists(keys.selection_groups(), group_name):
        return {"error": "Group not found"}, 404

    # Delete group
    r.hdel(keys.selection_groups(), group_name)

    # Clear active group if it was this one
    active_group = r.get(keys.active_selection_group())
    if active_group == group_name:
        r.delete(keys.active_selection_group())

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION_GROUPS,
        {"operation": "group_deleted", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200


@geometries.route(
    "/api/rooms/<string:room_id>/selections/groups/<string:group_name>/load",
    methods=["POST"],
)
def load_selection_group(room_id: str, group_name: str):
    """Load a selection group (apply it to current selections).

    This sets the active group and updates all selections to match the group.
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Get group
    group_data = r.hget(keys.selection_groups(), group_name)
    if group_data is None:
        return {"error": "Group not found"}, 404

    group = json.loads(group_data)

    # Apply group to current selections
    for geometry, indices in group.items():
        r.hset(keys.selections(), geometry, json.dumps(indices))

    # Set as active group
    r.set(keys.active_selection_group(), group_name)

    # Emit invalidation (all selections changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"operation": "group_loaded", "group": group_name},
        to=f"room:{room_id}",
    )

    return {"success": True}, 200
