"""Session API routes.

REST endpoints for frontend session management.
Sessions represent browser windows/tabs with independent cameras and settings.

Note: Python clients do NOT appear in sessions - only frontend browsers.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.auth import require_auth
from zndraw.geometries import Camera
from zndraw.server import socketio
from zndraw.settings import RoomConfig

from .constants import SocketEvents
from .redis_keys import RoomKeys

log = logging.getLogger(__name__)

session_bp = Blueprint("sessions", __name__)


@session_bp.route("/api/rooms/<string:room_id>/sessions", methods=["GET"])
@require_auth
def list_sessions(room_id: str):
    """List all frontend sessions in a room.

    Parameters
    ----------
    room_id : str
        Room identifier.

    Returns
    -------
    dict
        {"sessions": [{"session_id": str}, ...]}
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Get all frontend session IDs
    session_ids = r.smembers(keys.frontend_sessions())

    sessions = [{"session_id": session_id} for session_id in session_ids]

    log.debug(f"list_sessions: room={room_id}, count={len(sessions)}")
    return {"sessions": sessions}, 200


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/camera",
    methods=["GET"],
)
@require_auth
def get_session_camera(room_id: str, session_id: str):
    """Get camera state for a session.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Returns
    -------
    dict
        Camera state as JSON.
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Verify session exists
    if not r.sismember(keys.frontend_sessions(), session_id):
        return {"error": f"Session '{session_id}' not found"}, 404

    # Get camera state
    camera_json = r.hget(keys.session_cameras(), session_id)
    if camera_json is None:
        # Return default camera if not yet set
        camera = Camera()
        return {"camera": camera.model_dump()}, 200

    camera_data = json.loads(camera_json)
    log.debug(f"get_session_camera: room={room_id}, session={session_id}")
    return {"camera": camera_data}, 200


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/camera",
    methods=["PUT"],
)
@require_auth
def set_session_camera(room_id: str, session_id: str):
    """Set camera state for a session (programmatic control from Python).

    Validates the camera data and emits socket event to update frontend.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Request Body
    ------------
    JSON object with camera properties (position, target, fov, etc.)

    Returns
    -------
    dict
        {"status": "success"}
    """
    r = current_app.extensions["redis"]
    keys = RoomKeys(room_id)

    # Verify session exists
    if not r.sismember(keys.frontend_sessions(), session_id):
        return {"error": f"Session '{session_id}' not found"}, 404

    json_data = request.json
    if json_data is None:
        return {"error": "Request body must be JSON"}, 400

    # Validate camera data through Pydantic
    try:
        camera = Camera.model_validate(json_data)
    except Exception as e:
        return {"error": f"Invalid camera data: {e}"}, 400

    # Store in Redis
    r.hset(keys.session_cameras(), session_id, camera.model_dump_json())

    # Emit to frontend to update camera (programmatic control)
    # Use session-specific room to target only that session
    socketio.emit(
        SocketEvents.CAMERA_CONTROL,
        {
            "sessionId": session_id,
            "camera": camera.model_dump(),
        },
        to=f"session:{session_id}",
    )

    log.debug(f"set_session_camera: room={room_id}, session={session_id}")
    return {"status": "success"}, 200


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/settings",
    methods=["GET"],
)
@require_auth
def get_session_settings(room_id: str, session_id: str):
    """Get settings for a session.

    Returns both the JSON schema and current data for all settings categories.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Returns
    -------
    dict
        {"schema": RoomConfig schema, "data": all settings data}
    """
    settings_service = current_app.extensions["settings_service"]

    data = settings_service.get_all(room_id, session_id)
    schema = RoomConfig.model_json_schema()

    log.debug(f"get_session_settings: room={room_id}, session={session_id}")
    return {"schema": schema, "data": data}, 200


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/settings",
    methods=["PUT"],
)
@require_auth
def set_session_settings(room_id: str, session_id: str):
    """Update settings for a session.

    Accepts partial updates - only provided categories are updated.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Request Body
    ------------
    JSON object with category keys and settings data values, e.g.:
    {"studio_lighting": {"key_light": 0.8}}

    Returns
    -------
    dict
        {"status": "success"}
    """
    json_data = request.json
    if json_data is None:
        return {"error": "Request body must be JSON"}, 400

    # Validate categories
    valid_categories = set(RoomConfig.model_fields.keys())
    provided_categories = set(json_data.keys())
    invalid = provided_categories - valid_categories
    if invalid:
        return {"error": f"Unknown settings categories: {invalid}"}, 400

    settings_service = current_app.extensions["settings_service"]
    settings_service.update_all(room_id, session_id, json_data)

    # Emit invalidate event to notify clients
    socketio.emit(
        SocketEvents.INVALIDATE,
        {
            "sessionId": session_id,
            "category": "settings",
            "roomId": room_id,
        },
        to=f"room:{room_id}",
    )

    log.debug(f"set_session_settings: room={room_id}, session={session_id}")
    return {"status": "success"}, 200
