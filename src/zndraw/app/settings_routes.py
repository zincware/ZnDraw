"""Settings API routes.

Dedicated endpoints for user settings management.
Settings are per-user, per-room and stored via JWT authentication.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.auth import get_current_user, require_auth
from zndraw.server import socketio
from zndraw.settings import settings as settings_registry

from .constants import SocketEvents

log = logging.getLogger(__name__)

settings_bp = Blueprint("settings", __name__)


@settings_bp.route("/api/rooms/<string:room_id>/settings/schema", methods=["GET"])
def get_settings_schema(room_id: str):
    """Get the JSON schema for all settings categories.

    Returns schemas for camera, studio_lighting, pathtracing, property_inspector.
    Settings are always per-user, per-room (never public).

    Parameters
    ----------
    room_id : str
        Room identifier (unused for schema, but kept for API consistency)

    Returns
    -------
    list
        List of settings schemas with metadata
    """
    schema_list = []

    for name, cls in settings_registry.items():
        schema_list.append(
            {
                "name": name,
                "schema": cls.model_json_schema(),
                "provider": "settings",
            }
        )

    return schema_list


@settings_bp.route(
    "/api/rooms/<string:room_id>/settings/<string:category>",
    methods=["GET"],
)
@require_auth
def get_setting(room_id: str, category: str):
    """Get a specific settings category for the authenticated user.

    Parameters
    ----------
    room_id : str
        Room identifier
    category : str
        Settings category (camera, studio_lighting, pathtracing, property_inspector)

    Returns
    -------
    dict
        {"data": settings_data} - always returns valid data (defaults if not stored)
    """
    user_name = get_current_user()

    if category not in settings_registry:
        return {"error": f"Unknown settings category: {category}"}, 400

    settings_service = current_app.extensions["settings_service"]
    data = settings_service.get_category(room_id, user_name, category)

    # If no data stored, return defaults from Pydantic model
    if data is None:
        settings_class = settings_registry[category]
        data = settings_class().model_dump()
        log.debug(
            f"get_setting: room={room_id}, user={user_name}, category={category}, "
            f"returning defaults"
        )
    else:
        log.debug(
            f"get_setting: room={room_id}, user={user_name}, category={category}, "
            f"returning stored data"
        )

    return {"data": data}, 200


@settings_bp.route(
    "/api/rooms/<string:room_id>/settings/<string:category>",
    methods=["PUT"],
)
@require_auth
def update_setting(room_id: str, category: str):
    """Update a specific settings category for the authenticated user.

    Parameters
    ----------
    room_id : str
        Room identifier
    category : str
        Settings category (camera, studio_lighting, pathtracing, property_inspector)

    Request Body
    ------------
    JSON object with settings data

    Returns
    -------
    dict
        {"status": "success", "message": "Settings updated"}
    """
    user_name = get_current_user()

    if category not in settings_registry:
        return {"error": f"Unknown settings category: {category}"}, 400

    json_data = request.json
    if json_data is None:
        return {"error": "Request body must be JSON"}, 400

    # Allow both {"data": {...}} and direct {...} formats
    data = json_data.get("data", json_data)

    log.info(
        f"Settings update: room={room_id}, user={user_name}, "
        f"category={category}, data={json.dumps(data)}"
    )

    settings_service = current_app.extensions["settings_service"]
    settings_service.update_category(room_id, user_name, category, data)

    # Emit invalidate event to notify other clients (same user, same room)
    socketio.emit(
        SocketEvents.INVALIDATE,
        {
            "userName": user_name,
            "category": "settings",
            "extension": category,
            "roomId": room_id,
        },
        to=f"room:{room_id}",
    )

    log.info(f"Updated settings for room {room_id}, user {user_name}: {category}")

    return {"status": "success", "message": "Settings updated"}, 200
