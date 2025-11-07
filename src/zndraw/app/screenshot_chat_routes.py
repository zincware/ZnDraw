"""Screenshot and chat message routes.

Handles screenshot uploads, downloads, metadata, and chat message retrieval.
"""

import json
import logging

from flask import Blueprint, current_app, request, send_from_directory

from zndraw.screenshot_manager import ScreenshotManager
from zndraw.server import socketio

from .redis_keys import RoomKeys

log = logging.getLogger(__name__)

media = Blueprint("media", __name__)

# 10MB max file size for screenshots
MAX_SCREENSHOT_SIZE = 10 * 1024 * 1024


def _get_screenshot_manager(room_id: str) -> ScreenshotManager:
    """Helper to create ScreenshotManager instance."""
    config = current_app.extensions["config"]
    return ScreenshotManager(room_id, config.storage_path)


@media.route("/api/rooms/<string:room_id>/chat/messages", methods=["GET"])
def get_chat_messages(room_id: str):
    """
    Get paginated chat messages for a room.

    Query Parameters:
        - limit (int): Number of messages (default: 30, max: 100)
        - before (int): Get messages before this timestamp
        - after (int): Get messages after this timestamp

    Returns:
    {
        "messages": [Message],
        "metadata": {
            "hasMore": bool,
            "totalCount": int,
            "oldestTimestamp": int | null,
            "newestTimestamp": int | null
        }
    }
    """
    r = current_app.extensions["redis"]

    # Parse and validate query parameters
    try:
        limit = int(request.args.get("limit", 30))
        limit = max(1, min(limit, 100))  # Clamp between 1 and 100
    except (ValueError, TypeError):
        return {"error": "Invalid limit parameter"}, 400

    before = request.args.get("before")
    after = request.args.get("after")

    try:
        before = int(before) if before else None
        after = int(after) if after else None
    except (ValueError, TypeError):
        return {"error": "Invalid before/after parameter"}, 400

    room_keys = RoomKeys(room_id)
    index_key = room_keys.chat_index()
    data_key = room_keys.chat_data()

    # Get total count
    total_count = r.zcard(index_key)

    # Determine range query based on before/after
    if after is not None:
        # Get messages after timestamp (ascending order, then reverse)
        message_ids = r.zrangebyscore(
            index_key, f"({after}", "+inf", start=0, num=limit
        )
        # Reverse to maintain newest-first order
        message_ids = list(reversed(message_ids))
    elif before is not None:
        # Get messages before timestamp (descending order)
        message_ids = r.zrevrangebyscore(
            index_key, f"({before}", "-inf", start=0, num=limit
        )
    else:
        # Get latest messages (descending order)
        message_ids = r.zrevrangebyscore(index_key, "+inf", "-inf", start=0, num=limit)

    # Fetch message data
    messages = []
    if message_ids:
        with r.pipeline() as pipe:
            for msg_id in message_ids:
                pipe.hget(data_key, msg_id)
            message_data_list = pipe.execute()

        messages = [json.loads(msg_data) for msg_data in message_data_list if msg_data]

    # Calculate metadata
    has_more = False
    oldest_timestamp = None
    newest_timestamp = None

    if messages:
        oldest_timestamp = messages[-1]["createdAt"]
        newest_timestamp = messages[0]["createdAt"]

        # Check if there are more messages
        if after is not None:
            # Check if there are more recent messages
            count_after = r.zcount(index_key, f"({newest_timestamp}", "+inf")
            has_more = count_after > 0
        else:
            # Check if there are older messages
            count_before = r.zcount(index_key, "-inf", f"({oldest_timestamp}")
            has_more = count_before > 0

    return {
        "messages": messages,
        "metadata": {
            "hasMore": has_more,
            "totalCount": total_count,
            "oldestTimestamp": oldest_timestamp,
            "newestTimestamp": newest_timestamp,
        },
    }, 200


@media.route("/api/rooms/<string:room_id>/screenshots/upload", methods=["POST"])
def upload_screenshot(room_id: str):
    """Upload screenshot from frontend.

    Request:
        multipart/form-data with:
        - file: image file
        - format: png/jpeg/webp
        - width: optional image width
        - height: optional image height

    Response:
        JSON with screenshot metadata
    """
    if "file" not in request.files:
        return {"error": "No file provided"}, 400

    file = request.files["file"]
    if file.filename == "":
        return {"error": "Empty filename"}, 400

    format = request.form.get("format", "png")
    width = request.form.get("width", type=int)
    height = request.form.get("height", type=int)

    try:
        image_data = file.read()

        # Validate file size
        if len(image_data) > MAX_SCREENSHOT_SIZE:
            return {
                "error": f"File too large. Maximum size is {MAX_SCREENSHOT_SIZE // 1024 // 1024}MB"
            }, 400

        manager = _get_screenshot_manager(room_id)
        screenshot = manager.save(image_data, format, width, height)

        socketio.emit(
            "screenshot:created",
            {"id": screenshot.id},
            to=f"room:{room_id}",
        )

        return {
            "id": screenshot.id,
            "format": screenshot.format,
            "size": screenshot.size,
            "url": f"/api/rooms/{room_id}/screenshots/{screenshot.id}",
        }, 201
    except ValueError as e:
        return {"error": str(e)}, 400
    except Exception as e:
        return {"error": f"Failed to save screenshot: {str(e)}"}, 500


@media.route("/api/rooms/<string:room_id>/screenshots", methods=["GET"])
def list_screenshots(room_id: str):
    """List all screenshots for a room.

    Query params:
        limit: max results (default 20)
        offset: skip N results (default 0)

    Response:
        JSON with screenshots array and total count
    """
    limit = request.args.get("limit", 20, type=int)
    offset = request.args.get("offset", 0, type=int)

    if limit < 1 or limit > 100:
        return {"error": "Limit must be between 1 and 100"}, 400
    if offset < 0:
        return {"error": "Offset must be non-negative"}, 400

    try:
        manager = _get_screenshot_manager(room_id)
        screenshots = manager.list(limit, offset)
        total = manager.count()

        return {
            "screenshots": [
                {
                    "id": s.id,
                    "format": s.format,
                    "size": s.size,
                    "width": s.width,
                    "height": s.height,
                    "url": f"/api/rooms/{room_id}/screenshots/{s.id}",
                }
                for s in screenshots
            ],
            "total": total,
            "limit": limit,
            "offset": offset,
        }, 200
    except Exception as e:
        return {"error": f"Failed to list screenshots: {str(e)}"}, 500


@media.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>", methods=["GET"]
)
def get_screenshot(room_id: str, screenshot_id: int):
    """Download a specific screenshot.

    Returns the image file with appropriate Content-Type header.
    """
    try:
        manager = _get_screenshot_manager(room_id)
        result = manager.get(screenshot_id)

        if not result:
            return {"error": "Screenshot not found"}, 404

        filepath, metadata = result
        return send_from_directory(
            filepath.parent,
            filepath.name,
            mimetype=f"image/{metadata.format}",
        )
    except Exception as e:
        return {"error": f"Failed to get screenshot: {str(e)}"}, 500


@media.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>/metadata",
    methods=["GET"],
)
def get_screenshot_metadata(room_id: str, screenshot_id: int):
    """Get metadata for a specific screenshot."""
    try:
        manager = _get_screenshot_manager(room_id)
        result = manager.get(screenshot_id)

        if not result:
            return {"error": "Screenshot not found"}, 404

        _, metadata = result
        return {
            "id": metadata.id,
            "format": metadata.format,
            "size": metadata.size,
            "width": metadata.width,
            "height": metadata.height,
            "url": f"/api/rooms/{room_id}/screenshots/{metadata.id}",
        }, 200
    except Exception as e:
        return {"error": f"Failed to get screenshot metadata: {str(e)}"}, 500


@media.route(
    "/api/rooms/<string:room_id>/screenshots/<int:screenshot_id>", methods=["DELETE"]
)
def delete_screenshot(room_id: str, screenshot_id: int):
    """Delete a screenshot."""
    try:
        manager = _get_screenshot_manager(room_id)

        if manager.delete(screenshot_id):
            socketio.emit(
                "screenshot:deleted",
                {"id": screenshot_id},
                to=f"room:{room_id}",
            )
            return {"success": True}, 200

        return {"error": "Screenshot not found"}, 404
    except Exception as e:
        return {"error": f"Failed to delete screenshot: {str(e)}"}, 500
