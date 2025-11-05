"""Bookmark management routes.

Handles CRUD operations for frame bookmarks.
"""

import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import RoomKeys

log = logging.getLogger(__name__)

bookmarks = Blueprint("bookmarks", __name__)


@bookmarks.route("/api/rooms/<string:room_id>/bookmarks", methods=["GET"])
def get_all_bookmarks(room_id: str):
    """Get all bookmarks for a room.

    Returns:
        {"bookmarks": {1: "First Frame", 5: "Middle Frame"}}
    """
    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)
    bookmarks_raw = r.hgetall(room_keys.bookmarks())
    # Convert byte keys to integers
    bookmarks = {int(k): v for k, v in bookmarks_raw.items()}
    return {"bookmarks": bookmarks}, 200


@bookmarks.route("/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["GET"])
def get_bookmark(room_id: str, index: int):
    """Get a specific bookmark by frame index."""
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    room_keys = RoomKeys(room_id)

    # Check if frame index is valid
    frame_count = room_service.get_frame_count(room_id)
    if index < 0 or index >= frame_count:
        return {
            "error": f"Bookmark index {index} out of range (0-{frame_count - 1})",
            "type": "IndexError",
        }, 404

    label = r.hget(room_keys.bookmarks(), str(index))
    if label is None:
        return {
            "error": f"Bookmark at index {index} does not exist",
            "type": "KeyError",
        }, 404

    return {"index": index, "label": label}, 200


@bookmarks.route("/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["PUT"])
def set_bookmark(room_id: str, index: int):
    """Set or update a bookmark at a specific frame index.

    Body: {"label": "Frame Label"}
    """
    r = current_app.extensions["redis"]
    room_service = current_app.extensions["room_service"]
    room_keys = RoomKeys(room_id)
    data = request.get_json() or {}
    label = data.get("label")

    if not label or not isinstance(label, str):
        return {
            "error": "Bookmark label must be a non-empty string",
            "type": "ValueError",
        }, 400

    # Check if frame index is valid
    frame_count = room_service.get_frame_count(room_id)
    if index < 0 or index >= frame_count:
        return {
            "error": f"Bookmark index {index} out of range (0-{frame_count - 1})",
            "type": "IndexError",
        }, 400

    # Set the bookmark
    r.hset(room_keys.bookmarks(), str(index), label)

    # Emit invalidate event
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"index": index, "operation": "set"},
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200


@bookmarks.route(
    "/api/rooms/<string:room_id>/bookmarks/<int:index>", methods=["DELETE"]
)
def delete_bookmark(room_id: str, index: int):
    """Delete a bookmark at a specific frame index."""
    r = current_app.extensions["redis"]
    room_keys = RoomKeys(room_id)

    response = r.hdel(room_keys.bookmarks(), str(index))
    if response == 0:
        return {
            "error": f"Bookmark at index {index} does not exist",
            "type": "KeyError",
        }, 404

    # Emit invalidate event
    socketio.emit(
        SocketEvents.INVALIDATE_BOOKMARK,
        {"index": index, "operation": "delete"},
        to=f"room:{room_id}",
    )

    return {"status": "success"}, 200
