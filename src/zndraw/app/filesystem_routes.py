"""REST API endpoints for filesystem operations.

These endpoints use Socket.IO RPC to communicate with client workers
that have registered filesystems.
"""
import logging
import uuid

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user

from .redis_keys import FilesystemKeys

log = logging.getLogger(__name__)

filesystem_bp = Blueprint("filesystem", __name__)

# Configuration constants
SOCKET_IO_TIMEOUT = 30  # seconds


def _require_authentication():
    """Require JWT authentication for filesystem endpoints.

    Returns
    -------
    str
        Username of authenticated user

    Raises
    ------
    AuthError
        If authentication fails (handled by caller)
    """
    return get_current_user()


def _list_filesystems_impl(room_id: str | None):
    """Shared implementation for listing filesystems.

    Parameters
    ----------
    room_id : str | None
        Room ID to get filesystems for, or None for global only

    Returns
    -------
    list[dict]
        List of filesystem metadata dicts
    """
    r = current_app.extensions["redis"]
    filesystems = []

    # Get room-scoped filesystems if room_id provided
    if room_id:
        pattern = f"room:{room_id}:filesystems:*"
        for key in r.scan_iter(pattern):
            if key.endswith(":worker"):
                continue
            if r.type(key) != "hash":
                continue

            metadata = r.hgetall(key)
            if metadata:
                filesystems.append({
                    "name": metadata.get("name"),
                    "fsType": metadata.get("fsType"),
                    "userName": metadata.get("userName"),
                    "public": metadata.get("public") == "true",
                    "workerId": metadata.get("workerId"),
                })

    # Get global filesystems
    pattern = "global:filesystems:*"
    for key in r.scan_iter(pattern):
        if key.endswith(":worker"):
            continue
        if r.type(key) != "hash":
            continue

        metadata = r.hgetall(key)
        if metadata:
            filesystems.append({
                "name": metadata.get("name"),
                "fsType": metadata.get("fsType"),
                "userName": metadata.get("userName"),
                "public": metadata.get("public") == "true",
                "workerId": metadata.get("workerId"),
            })

    return filesystems


def _list_files_impl(room_id: str | None, fs_name: str):
    """Shared implementation for listing files from a filesystem.

    Parameters
    ----------
    room_id : str | None
        Room ID for room-scoped filesystem, or None for global
    fs_name : str
        Filesystem name

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from zndraw.server import socketio

    r = current_app.extensions["redis"]

    # Find the worker for this filesystem
    if room_id:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        worker_id = r.get(keys.worker)
        # Fallback to global if not found in room
        if not worker_id:
            keys = FilesystemKeys.for_global_filesystem(fs_name)
            worker_id = r.get(keys.worker)
    else:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        worker_id = r.get(keys.worker)

    if not worker_id:
        fs_type = "global" if not room_id else "room-scoped or global"
        return jsonify({"error": f"{fs_type.capitalize()} filesystem '{fs_name}' not found"}), 404

    # Generate request ID
    request_id = str(uuid.uuid4())

    # Parse query parameters
    path = request.args.get("path", "")
    recursive = request.args.get("recursive", "false").lower() == "true"
    filter_extensions_str = request.args.get("filterExtensions")
    filter_extensions = (
        filter_extensions_str.split(",") if filter_extensions_str else None
    )
    search = request.args.get("search")

    # Send request to client via Socket.IO with timeout
    try:
        response = socketio.call(
            "filesystem:list",
            {
                "requestId": request_id,
                "fsName": fs_name,
                "path": path,
                "recursive": recursive,
                "filterExtensions": filter_extensions,
                "search": search,
            },
            to=worker_id,
            timeout=SOCKET_IO_TIMEOUT,
        )

        if not response:
            return jsonify({"error": "No response from worker"}), 500

        if not response.get("success"):
            return jsonify({"error": response.get("error", "Unknown error")}), 500

        return jsonify({"files": response.get("files", [])})

    except Exception as e:
        log.error(f"Error listing files from filesystem '{fs_name}': {e}")
        return jsonify({"error": str(e)}), 500


def _load_file_impl(room_id: str | None, fs_name: str, require_target_room: bool = False):
    """Shared implementation for loading files from a filesystem.

    Parameters
    ----------
    room_id : str | None
        Room ID for room-scoped filesystem, or None for global
    fs_name : str
        Filesystem name
    require_target_room : bool
        Whether targetRoom is required in request body

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from zndraw.server import socketio

    r = current_app.extensions["redis"]

    data = request.get_json()
    if not data:
        return jsonify({"error": "Missing request body"}), 400

    path = data.get("path")
    if not path:
        return jsonify({"error": "Missing 'path' in request body"}), 400

    # Handle target room
    target_room = data.get("targetRoom")
    if require_target_room and not target_room:
        return jsonify({"error": "Missing 'targetRoom' in request body"}), 400
    if not target_room:
        target_room = room_id

    batch_size = data.get("batchSize", 10)
    start = data.get("start")
    stop = data.get("stop")
    step = data.get("step")

    # Find the worker for this filesystem
    if room_id:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        worker_id = r.get(keys.worker)
        # Fallback to global if not found in room
        if not worker_id:
            keys = FilesystemKeys.for_global_filesystem(fs_name)
            worker_id = r.get(keys.worker)
    else:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        worker_id = r.get(keys.worker)

    if not worker_id:
        fs_type = "global" if not room_id else "room-scoped or global"
        return jsonify({"error": f"{fs_type.capitalize()} filesystem '{fs_name}' not found"}), 404

    # Generate request ID
    request_id = str(uuid.uuid4())

    # Send request to client via Socket.IO with timeout
    try:
        response = socketio.call(
            "filesystem:load",
            {
                "requestId": request_id,
                "fsName": fs_name,
                "path": path,
                "room": target_room,
                "batchSize": batch_size,
                "start": start,
                "stop": stop,
                "step": step,
            },
            to=worker_id,
            timeout=SOCKET_IO_TIMEOUT,
        )

        if not response:
            return jsonify({"error": "No response from worker"}), 500

        if not response.get("success"):
            return jsonify({"error": response.get("error", "Unknown error")}), 500

        return jsonify({
            "success": True,
            "frameCount": response.get("frameCount", 0),
        })

    except Exception as e:
        log.error(f"Error loading file from filesystem '{fs_name}': {e}")
        return jsonify({"error": str(e)}), 500


@filesystem_bp.route("/api/rooms/<room_id>/filesystems", methods=["GET"])
def list_filesystems(room_id: str):
    """List all registered filesystems for a room.

    Returns both room-scoped and global filesystems.

    Requires JWT authentication.

    Returns
    -------
    {
        "filesystems": [
            {
                "name": str,
                "fsType": str,
                "userName": str,
                "public": bool,
                "workerId": str
            }
        ]
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    filesystems = _list_filesystems_impl(room_id)
    return jsonify({"filesystems": filesystems})


@filesystem_bp.route("/api/rooms/<room_id>/filesystems/<fs_name>/list", methods=["GET"])
def list_files(room_id: str, fs_name: str):
    """List files from a registered filesystem.

    Requires JWT authentication.

    Query Parameters
    ----------------
    path : str, optional
        Path to list files from (default: "")
    recursive : bool, optional
        Whether to list recursively (default: False)
    filterExtensions : str, optional
        Comma-separated list of extensions to filter by (e.g., "xyz,h5")
    search : str, optional
        Regex pattern to filter file names

    Returns
    -------
    {
        "files": [
            {
                "name": str,
                "path": str,
                "size": int,
                "type": str,
                "modified": float | None
            }
        ]
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    return _list_files_impl(room_id, fs_name)


@filesystem_bp.route(
    "/api/rooms/<room_id>/filesystems/<fs_name>/metadata", methods=["GET"]
)
def get_file_metadata(room_id: str, fs_name: str):
    """Get metadata for a specific file.

    Requires JWT authentication.

    Query Parameters
    ----------------
    path : str, required
        Path to the file

    Returns
    -------
    {
        "metadata": {
            "name": str,
            "path": str,
            "size": int,
            "type": str,
            "modified": float | None,
            "created": float | None
        }
    }
    """
    from zndraw.server import socketio

    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    r = current_app.extensions["redis"]

    path = request.args.get("path")
    if not path:
        return jsonify({"error": "Missing 'path' query parameter"}), 400

    # Find the worker for this filesystem
    keys = FilesystemKeys.for_filesystem(room_id, fs_name)
    worker_id = r.get(keys.worker)

    if not worker_id:
        # Try global filesystem
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        worker_id = r.get(keys.worker)

    if not worker_id:
        return jsonify({"error": f"Filesystem '{fs_name}' not found"}), 404

    # Generate request ID
    request_id = str(uuid.uuid4())

    # Send request to client via Socket.IO with timeout
    try:
        response = socketio.call(
            "filesystem:metadata",
            {"requestId": request_id, "fsName": fs_name, "path": path},
            to=worker_id,
            timeout=SOCKET_IO_TIMEOUT,
        )

        if not response:
            return jsonify({"error": "No response from worker"}), 500

        if not response.get("success"):
            return jsonify({"error": response.get("error", "Unknown error")}), 500

        return jsonify({"metadata": response.get("metadata")})

    except Exception as e:
        log.error(f"Error getting metadata from filesystem '{fs_name}': {e}")
        return jsonify({"error": str(e)}), 500


@filesystem_bp.route("/api/rooms/<room_id>/filesystems/<fs_name>/load", methods=["POST"])
def load_file(room_id: str, fs_name: str):
    """Load a file from a filesystem and upload to the target room.

    Requires JWT authentication.

    Request Body
    ------------
    {
        "path": str,  # Required
        "targetRoom": str,  # Optional, defaults to current room
        "batchSize": int,  # Optional, default 10
        "start": int,  # Optional slice parameter
        "stop": int,  # Optional slice parameter
        "step": int  # Optional slice parameter
    }

    Returns
    -------
    {
        "success": true,
        "frameCount": int
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    return _load_file_impl(room_id, fs_name, require_target_room=False)


# Global filesystem endpoints (no room_id required)


@filesystem_bp.route("/api/filesystems", methods=["GET"])
def list_global_filesystems():
    """List all registered global/public filesystems.

    Requires JWT authentication.

    Returns
    -------
    {
        "filesystems": [
            {
                "name": str,
                "fsType": str,
                "userName": str,
                "public": bool,
                "workerId": str
            }
        ]
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    filesystems = _list_filesystems_impl(None)
    return jsonify({"filesystems": filesystems})


@filesystem_bp.route("/api/filesystems/<fs_name>/list", methods=["GET"])
def list_global_files(fs_name: str):
    """List files from a global filesystem.

    Requires JWT authentication.

    Query Parameters
    ----------------
    path : str, optional
        Path to list files from (default: "")
    recursive : bool, optional
        Whether to list recursively (default: False)
    filterExtensions : str, optional
        Comma-separated list of extensions to filter by (e.g., "xyz,h5")
    search : str, optional
        Regex pattern to filter file names

    Returns
    -------
    {
        "files": [
            {
                "name": str,
                "path": str,
                "size": int,
                "type": str,
                "modified": float | None
            }
        ]
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    return _list_files_impl(None, fs_name)


@filesystem_bp.route("/api/filesystems/<fs_name>/load", methods=["POST"])
def load_global_file(fs_name: str):
    """Load a file from a global filesystem and upload to a room.

    Requires JWT authentication.

    Request Body
    ------------
    {
        "path": str,  # Required
        "targetRoom": str,  # Required
        "batchSize": int,  # Optional, default 10
        "start": int,  # Optional slice parameter
        "stop": int,  # Optional slice parameter
        "step": int  # Optional slice parameter
    }

    Returns
    -------
    {
        "success": true,
        "frameCount": int
    }
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    return _load_file_impl(None, fs_name, require_target_room=True)
