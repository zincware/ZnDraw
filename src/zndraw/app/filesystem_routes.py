"""REST API endpoints for filesystem operations.

These endpoints use Socket.IO RPC to communicate with client workers
that have registered filesystems.
"""
import logging
import uuid

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import AuthError, get_current_user

from .redis_keys import FilesystemKeys, RoomKeys

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
        room_keys = RoomKeys(room_id)
        pattern = room_keys.filesystems_pattern()
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
                    "public": metadata.get("public") == "true",
                    "sessionId": metadata.get("sessionId"),
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
                "public": metadata.get("public") == "true",
                "sessionId": metadata.get("sessionId"),
            })

    return filesystems


def _list_files_impl(fs_name: str, public: bool, room_id: str | None = None):
    """Shared implementation for listing files from a filesystem.

    Parameters
    ----------
    fs_name : str
        Filesystem name
    public : bool
        Whether this is a global (public=True) or room-scoped (public=False) filesystem
    room_id : str | None
        Room ID, required if public=False

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from zndraw.server import socketio

    r = current_app.extensions["redis"]

    # Validate parameters
    if not public and not room_id:
        return jsonify({"error": "room_id required for room-scoped filesystem"}), 400

    # Find the worker for this filesystem
    if public:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        fs_type = "global"
    else:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        fs_type = "room-scoped"

    worker_id = r.get(keys.worker)
    if not worker_id:
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
                "public": public,
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


def _load_file_impl(fs_name: str, public: bool, room_id: str | None = None, require_target_room: bool = False):
    """Shared implementation for loading files from a filesystem.

    Parameters
    ----------
    fs_name : str
        Filesystem name
    public : bool
        Whether this is a global (public=True) or room-scoped (public=False) filesystem
    room_id : str | None
        Room ID, required if public=False
    require_target_room : bool
        Whether targetRoom is required in request body

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from zndraw.server import socketio

    r = current_app.extensions["redis"]

    # Validate parameters
    if not public and not room_id:
        return jsonify({"error": "room_id required for room-scoped filesystem"}), 400

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
    if public:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        fs_type = "global"
    else:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        fs_type = "room-scoped"

    worker_id = r.get(keys.worker)
    if not worker_id:
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
                "public": public,
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
    [
        {
            "name": str,
            "fsType": str,
            "public": bool,
            "sessionId": str
        }
    ]
    """
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    filesystems = _list_filesystems_impl(room_id)
    return jsonify(filesystems)


@filesystem_bp.route("/api/rooms/<room_id>/filesystems/<fs_name>/list", methods=["GET"])
def list_files(room_id: str, fs_name: str):
    """List files from a room-scoped filesystem.

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

    return _list_files_impl(fs_name=fs_name, public=False, room_id=room_id)


def _get_file_metadata_impl(fs_name: str, public: bool, room_id: str | None = None):
    """Shared implementation for getting file metadata.

    Parameters
    ----------
    fs_name : str
        Filesystem name
    public : bool
        Whether this is a global (public=True) or room-scoped (public=False) filesystem
    room_id : str | None
        Room ID, required if public=False

    Returns
    -------
    tuple[dict, int]
        JSON response and HTTP status code
    """
    from zndraw.server import socketio

    r = current_app.extensions["redis"]

    # Validate parameters
    if not public and not room_id:
        return jsonify({"error": "room_id required for room-scoped filesystem"}), 400

    path = request.args.get("path")
    if not path:
        return jsonify({"error": "Missing 'path' query parameter"}), 400

    # Find the worker for this filesystem
    if public:
        keys = FilesystemKeys.for_global_filesystem(fs_name)
        fs_type = "global"
    else:
        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        fs_type = "room-scoped"

    worker_id = r.get(keys.worker)
    if not worker_id:
        return jsonify({"error": f"{fs_type.capitalize()} filesystem '{fs_name}' not found"}), 404

    # Generate request ID
    request_id = str(uuid.uuid4())

    # Send request to client via Socket.IO with timeout
    try:
        response = socketio.call(
            "filesystem:metadata",
            {"requestId": request_id, "fsName": fs_name, "public": public, "path": path},
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


@filesystem_bp.route(
    "/api/rooms/<room_id>/filesystems/<fs_name>/metadata", methods=["GET"]
)
def get_file_metadata(room_id: str, fs_name: str):
    """Get metadata for a specific file from a room-scoped filesystem.

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
    try:
        _require_authentication()
    except AuthError as e:
        return jsonify({"error": e.message}), e.status_code

    return _get_file_metadata_impl(fs_name=fs_name, public=False, room_id=room_id)


@filesystem_bp.route("/api/rooms/<room_id>/filesystems/<fs_name>/load", methods=["POST"])
def load_file(room_id: str, fs_name: str):
    """Load a file from a room-scoped filesystem and upload to the target room.

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

    return _load_file_impl(fs_name=fs_name, public=False, room_id=room_id, require_target_room=False)


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

    return _list_files_impl(fs_name=fs_name, public=True)


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

    return _load_file_impl(fs_name=fs_name, public=True, require_target_room=True)
