"""Settings API routes.

Consolidated endpoints for user settings management.
Settings are per-user, per-room and stored via JWT authentication.
"""

import logging

from flask import Blueprint, current_app, request

from zndraw.auth import get_current_user, require_auth
from zndraw.server import socketio
from zndraw.settings import RoomConfig

from .constants import SocketEvents

log = logging.getLogger(__name__)

settings_bp = Blueprint("settings", __name__)


@settings_bp.route("/api/rooms/<string:room_id>/settings", methods=["GET"])
@require_auth
def get_settings(room_id: str):
    """Get all settings with schema for the authenticated user.

    Returns both the JSON schema and current data for all settings categories.

    Parameters
    ----------
    room_id : str
        Room identifier

    Returns
    -------
    dict
        {"schema": RoomConfig schema, "data": all settings data}
    """
    user_name = get_current_user()
    settings_service = current_app.extensions["settings_service"]

    data = settings_service.get_all(room_id, user_name)
    schema = RoomConfig.model_json_schema()

    log.debug(f"get_settings: room={room_id}, user={user_name}")
    return {"schema": schema, "data": data}, 200


@settings_bp.route("/api/rooms/<string:room_id>/settings", methods=["PUT"])
@require_auth
def update_settings(room_id: str):
    """Update settings categories for the authenticated user.

    Accepts partial updates - only provided categories are updated.

    Parameters
    ----------
    room_id : str
        Room identifier

    Request Body
    ------------
    JSON object with category keys and settings data values, e.g.:
    {"camera": {"near_plane": 0.5}, "studio_lighting": {"key_light": 0.8}}

    Returns
    -------
    dict
        {"status": "success"}
    """
    user_name = get_current_user()
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
    log.debug(f"update_settings received categories: {list(json_data.keys())}")
    settings_service.update_all(room_id, user_name, json_data)

    # Emit invalidate event to notify other clients (same user, same room)
    socketio.emit(
        SocketEvents.INVALIDATE,
        {
            "userName": user_name,
            "category": "settings",
            "roomId": room_id,
        },
        to=f"room:{room_id}",
    )

    log.debug(f"Updated settings for room {room_id}: {list(json_data.keys())}")

    return {"status": "success"}, 200
