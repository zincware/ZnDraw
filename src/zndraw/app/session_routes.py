"""Session API routes.

REST endpoints for frontend session management.
Sessions represent browser windows/tabs with independent settings.

Note: Python clients do NOT appear in sessions - only frontend browsers.
Note: Session cameras are handled as regular geometries (key: cam:session:{session_id}).
"""

import json
import logging
import time
import uuid

from flask import Blueprint, current_app, request
from pydantic import BaseModel, Field, ValidationError

from zndraw.auth import get_current_user, require_auth
from zndraw.server import socketio
from zndraw.settings import RoomConfig

from .constants import SocketEvents
from .redis_keys import RoomKeys, SessionKeys

log = logging.getLogger(__name__)


class ScreenshotRequest(BaseModel):
    """Request model for screenshot endpoint."""

    timeout: float = Field(default=10.0, ge=0.0, le=60.0)


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


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/active-camera",
    methods=["GET"],
)
@require_auth
def get_active_camera(room_id: str, session_id: str):
    """Get active camera key for a session.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Returns
    -------
    dict
        {"active_camera": str}
    """
    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)
    key = room_keys.session_active_camera(session_id)
    value = r.get(key)

    log.debug(f"get_active_camera: room={room_id}, session={session_id}, value={value}")
    return {"active_camera": value}, 200


@session_bp.route(
    "/api/rooms/<string:room_id>/sessions/<string:session_id>/active-camera",
    methods=["PUT"],
)
@require_auth
def set_active_camera(room_id: str, session_id: str):
    """Set active camera key for a session.

    Parameters
    ----------
    room_id : str
        Room identifier.
    session_id : str
        Session identifier.

    Request Body
    ------------
    JSON object with active_camera key, e.g.:
    {"active_camera": "my_camera"}

    Returns
    -------
    dict
        {"status": "success"}
    """
    json_data = request.json
    if json_data is None:
        return {"error": "Request body must be JSON"}, 400

    active_camera = json_data.get("active_camera")
    if active_camera is None:
        return {"error": "active_camera is required"}, 400

    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)
    key = room_keys.session_active_camera(session_id)
    r.set(key, active_camera)

    # Emit to ONLY this session
    socketio.emit(
        SocketEvents.ACTIVE_CAMERA_UPDATE,
        {"active_camera": active_camera},
        to=f"session:{session_id}",
    )

    log.debug(
        f"set_active_camera: room={room_id}, session={session_id}, camera={active_camera}"
    )
    return {"status": "success"}, 200


@session_bp.route(
    "/api/sessions/<string:session_id>/screenshot",
    methods=["POST"],
)
@require_auth
def request_screenshot(session_id: str):
    """Request a screenshot from a frontend session.

    Triggers a screenshot capture in the browser session via Socket.IO.
    Browser uploads the screenshot via REST, which marks the request as complete.
    This endpoint polls for completion and returns the screenshot info.

    Parameters
    ----------
    session_id : str
        Session identifier of the browser to capture.

    Request Body (optional)
    -----------------------
    JSON object with optional settings:
    - timeout: Maximum wait time in seconds (default: 10, max: 60)

    Returns
    -------
    dict
        {"screenshot_id": int, "room_id": str}

    HTTP Status Codes
    -----------------
    200: Screenshot captured successfully
    400: Invalid request (e.g., invalid timeout)
    404: Session not found
    408: Timeout waiting for screenshot
    """
    r = current_app.extensions["redis"]

    session_key = SessionKeys.session_data(session_id)
    session_data_str = r.get(session_key)
    if not session_data_str:
        return {"error": f"Session '{session_id}' not found"}, 404

    session_data = json.loads(session_data_str)
    room_id = session_data.get("roomId")
    if not room_id:
        return {"error": "Session has no associated room"}, 400

    room_keys = RoomKeys(room_id)
    current_user = get_current_user()
    if not r.sismember(room_keys.users(), current_user):
        return {"error": "You are not a member of this room"}, 403

    try:
        request_model = ScreenshotRequest.model_validate(
            request.get_json(silent=True) or {}
        )
    except ValidationError as e:
        return {"error": e.errors()}, 400
    timeout = request_model.timeout

    request_id = str(uuid.uuid4())
    request_key = room_keys.screenshot_request(request_id)
    request_data = {
        "session_id": session_id,
        "room_id": room_id,
        "status": "pending",
    }
    r.setex(request_key, 120, json.dumps(request_data))

    socketio.emit(
        "screenshot:request",
        {
            "requestId": request_id,
            "uploadUrl": f"/api/rooms/{room_id}/screenshots/upload",
        },
        to=f"session:{session_id}",
    )

    log.debug(f"request_screenshot: session={session_id}, request_id={request_id}")

    poll_interval = 0.1
    max_polls = int(timeout / poll_interval)

    for _ in range(max_polls):
        data_str = r.get(request_key)
        if data_str:
            data = json.loads(data_str)
            if data.get("status") == "completed":
                screenshot_id = data.get("screenshot_id")
                r.delete(request_key)
                return {
                    "screenshot_id": screenshot_id,
                    "room_id": room_id,
                }, 200
        time.sleep(poll_interval)

    r.delete(request_key)
    return {"error": "Timeout waiting for screenshot"}, 408
