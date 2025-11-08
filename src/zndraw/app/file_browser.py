"""File browser endpoints for browsing and loading local files."""

import logging
import os
from datetime import datetime, timezone
from pathlib import Path

from flask import Blueprint, current_app, jsonify, request

from zndraw.auth import require_auth

from .redis_keys import RoomKeys

log = logging.getLogger(__name__)

file_browser = Blueprint("file_browser", __name__, url_prefix="/api/file-browser")


def _get_all_rooms_metadata(redis_client) -> dict[str, dict]:
    """Get metadata for all rooms.

    Parameters
    ----------
    redis_client
        Redis client instance.

    Returns
    -------
    dict[str, dict]
        Dictionary mapping room_id to room metadata and description.
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    rooms = {}
    for key in redis_client.scan_iter(match="room:*:metadata"):
        room_id = key.split(":")[1]
        manager = RoomMetadataManager(redis_client, room_id)
        metadata = manager.get_all()

        # Get room description
        room_keys = RoomKeys(room_id)
        description = redis_client.get(room_keys.description())

        rooms[room_id] = {
            "room_id": room_id,
            "metadata": metadata,
            "description": description,
        }

    return rooms


def _find_matching_room(
    file_path: Path, root_path: str, rooms_metadata: dict[str, dict]
) -> dict | None:
    """Find a room that contains this exact file.

    Parameters
    ----------
    file_path : Path
        The file path to search for.
    root_path : str
        The root path for file browser.
    rooms_metadata : dict[str, dict]
        Pre-loaded metadata for all rooms.

    Returns
    -------
    dict | None
        Room info if match found, None otherwise.
    """
    # Compute relative path from root
    try:
        relative_path = str(file_path.relative_to(Path(root_path).resolve()))
    except ValueError:
        # If file is not relative to root, use absolute path
        relative_path = str(file_path)

    # Get current file stats
    try:
        current_stat = file_path.stat()
        current_size = str(current_stat.st_size)
        current_mtime = datetime.fromtimestamp(
            current_stat.st_mtime, tz=timezone.utc
        ).isoformat()
    except OSError:
        return None

    # Check each room's metadata
    for room_id, room_info in rooms_metadata.items():
        metadata = room_info["metadata"]

        stored_path = metadata.get("relative_file_path")
        stored_size = metadata.get("file_size")
        stored_mtime = metadata.get("file_modified_timestamp")

        # Check exact match: same path, size, and mtime
        if (
            stored_path == relative_path
            and stored_size == current_size
            and stored_mtime == current_mtime
        ):
            return room_info

    return None


def _find_room_with_exact_file(
    file_path: Path, redis_client, root_path: str
) -> str | None:
    """Find a room that has this exact file (matching path, size, and mtime).

    Parameters
    ----------
    file_path : Path
        The absolute file path to search for.
    redis_client
        Redis client instance.
    root_path : str
        The root path for file browser to compute relative path.

    Returns
    -------
    str | None
        Room ID if exact match found, None otherwise.
    """
    from zndraw.app.metadata_manager import RoomMetadataManager

    # Compute relative path from root
    try:
        relative_path = str(file_path.relative_to(Path(root_path).resolve()))
    except ValueError:
        # If file is not relative to root, use absolute path
        relative_path = str(file_path)

    # Get current file stats
    current_stat = file_path.stat()
    current_size = str(current_stat.st_size)
    current_mtime = datetime.fromtimestamp(
        current_stat.st_mtime, tz=timezone.utc
    ).isoformat()

    log.info(
        f"Looking for duplicate file: relative_path={relative_path}, size={current_size}, mtime={current_mtime}"
    )

    # Scan all rooms for matching metadata
    for key in redis_client.scan_iter(match="room:*:metadata"):
        room_id = key.split(":")[1]
        manager = RoomMetadataManager(redis_client, room_id)
        metadata = manager.get_all()

        stored_path = metadata.get("relative_file_path")
        stored_size = metadata.get("file_size")
        stored_mtime = metadata.get("file_modified_timestamp")

        log.debug(
            f"Checking room {room_id}: stored_path={stored_path}, stored_size={stored_size}, stored_mtime={stored_mtime}"
        )

        # Check exact match: same path, size, and mtime
        if (
            stored_path == relative_path
            and stored_size == current_size
            and stored_mtime == current_mtime
        ):
            log.info(f"Found duplicate in room: {room_id}")
            return room_id

    log.info("No duplicate found")
    return None


# File format to backend/library mapping
# Maps file extensions (without dot) to the backend libraries that handle them
FORMAT_BACKENDS = {
    "h5": ["ZnH5MD"],
    "h5md": ["ZnH5MD"],
    "hdf5": ["ZnH5MD"],
    "xyz": ["ASE"],
    "extxyz": ["ASE"],
    "pdb": ["ASE"],
    "traj": ["ASE"],
    "nc": ["ASE"],
    "cif": ["ASE"],
    "mol": ["ASE"],
    "sdf": ["ASE"],
    "db": ["ASE-DB"],
    "aselmdb": ["ASE-DB"],
    "gro": ["ASE"],
    "car": ["ASE"],
    "xsf": ["ASE"],
    "cube": ["ASE"],
    "json": ["ASE-DB"],
    "vasp": ["ASE"],
    "poscar": ["ASE"],
    "contcar": ["ASE"],
    "xdatcar": ["ASE"],
    "outcar": ["ASE"],
    "xml": ["ASE"],
    "pwi": ["ASE"],
    "pwo": ["ASE"],
    "out": ["ASE"],
    "castep": ["ASE"],
    "cell": ["ASE"],
    "geom": ["ASE"],
    "md": ["ASE"],
    "gjf": ["ASE"],
    "com": ["ASE"],
    "log": ["ASE"],
    "arc": ["ASE"],
    "dmol": ["ASE"],
    "gen": ["ASE"],
    "g96": ["ASE"],
    "xtl": ["ASE"],
    "rmc6f": ["ASE"],
    "shelx": ["ASE"],
    "res": ["ASE"],
    "vtu": ["ASE"],
    "vti": ["ASE"],
    "x3d": ["ASE"],
    "xsd": ["ASE"],
    "xtd": ["ASE"],
    "dcd": ["ASE"],
    "restart": ["ASE"],
    "dat": ["ASE"],
    "config": ["ASE"],
    "phonon": ["ASE"],
    "cfg": ["ASE"],
    "cjson": ["ASE"],
    "f34": ["ASE"],
    "con": ["ASE"],
    "nwi": ["ASE"],
    "nwo": ["ASE"],
}


def get_format_info(extension: str) -> str | None:
    """
    Get backend information for a file extension.

    Parameters
    ----------
    extension : str
        File extension (with or without leading dot).

    Returns
    -------
    str | None
        Backend information string or None if not in known formats.
    """
    ext = extension.lstrip(".")
    backends = FORMAT_BACKENDS.get(ext)
    if backends:
        return f"Supported by {', '.join(backends)}"
    return None


def validate_path(requested_path: str, root: str) -> Path | None:
    """
    Validate that requested path doesn't escape root directory.

    Parameters
    ----------
    requested_path : str
        The requested relative path.
    root : str
        The root directory.

    Returns
    -------
    Path | None
        The resolved path if valid, None if path escapes root.
    """
    try:
        root_path = Path(root).resolve()
        target_path = (root_path / requested_path).resolve()

        # Check if target is within root
        if not target_path.is_relative_to(root_path):
            return None

        return target_path
    except (ValueError, OSError) as e:
        log.warning(f"Path validation error: {e}")
        return None


def is_supported_file(filepath: Path) -> bool:
    """
    Check if file has a known supported extension.

    Note: This is informational only. The backend will attempt to read
    unknown formats via ASE, which may succeed for additional formats.

    Parameters
    ----------
    filepath : Path
        The file path to check.

    Returns
    -------
    bool
        True if file extension is in known supported formats.
    """
    ext = filepath.suffix.lstrip(".").lower()
    return ext in FORMAT_BACKENDS


def check_feature_enabled() -> tuple[dict[str, str], int] | None:
    """
    Check if file browser feature is enabled.

    Returns
    -------
    tuple[dict[str, str], int] | None
        Error response if disabled, None if enabled.
    """
    config = current_app.extensions.get("config")
    if not config or not config.file_browser_enabled:
        return {"error": "File browser feature is not enabled"}, 403
    return None


@file_browser.route("/list", methods=["GET"])
@require_auth
def list_directory():
    """
    List contents of a directory.

    Query Parameters
    ----------------
    path : str, optional
        Relative path from browser root. Defaults to root.
    search : str, optional
        Regex pattern to filter file/directory names.

    Returns
    -------
    Response
        JSON response with directory contents or error.
    """
    import re

    # Check if feature is enabled
    error = check_feature_enabled()
    if error:
        return jsonify(error[0]), error[1]

    # Get requested path and search pattern
    requested_path = request.args.get("path", "")
    search_pattern = request.args.get("search")
    config = current_app.extensions["config"]
    root = config.file_browser_root

    # Validate path
    target_path = validate_path(requested_path, root)
    if target_path is None:
        return jsonify({"error": "Invalid path or path outside root directory"}), 400

    # Check if path exists
    if not target_path.exists():
        return jsonify({"error": "Path does not exist"}), 404

    # Check if path is a directory
    if not target_path.is_dir():
        return jsonify({"error": "Path is not a directory"}), 400

    # Get redis client and preload all room metadata for file matching
    redis_client = current_app.extensions["redis"]
    rooms_metadata = _get_all_rooms_metadata(redis_client)

    # List directory contents
    try:
        items = []
        for item in sorted(
            target_path.iterdir(), key=lambda x: (not x.is_dir(), x.name)
        ):
            # Skip hidden files
            if item.name.startswith("."):
                continue

            # Apply search filter if provided
            if search_pattern:
                try:
                    pattern = re.compile(search_pattern, re.IGNORECASE)
                    if not pattern.search(item.name):
                        continue
                except re.error:
                    # Invalid regex, include all items
                    pass

            item_info = {
                "name": item.name,
                "type": "directory" if item.is_dir() else "file",
                "size": item.stat().st_size if item.is_file() else None,
                "modified": datetime.fromtimestamp(
                    item.stat().st_mtime, tz=timezone.utc
                ).isoformat(),
            }

            # Mark if file is supported and add format info
            if item.is_file():
                item_info["supported"] = is_supported_file(item)
                # Add format information for all files
                format_info = get_format_info(item.suffix)
                if format_info:
                    item_info["format_info"] = format_info

                # Check if this file is already loaded in any room
                if item_info["supported"]:
                    loaded_room = _find_matching_room(item, root, rooms_metadata)
                    if loaded_room:
                        item_info["alreadyLoaded"] = {
                            "room": loaded_room["room_id"],
                            "description": loaded_room.get("description"),
                        }

            items.append(item_info)

        # Get parent path
        parent = None
        if target_path != Path(root).resolve():
            parent_path = target_path.parent
            parent = str(parent_path.relative_to(Path(root).resolve()))

        return jsonify(
            {
                "current_path": str(target_path.relative_to(Path(root).resolve())),
                "items": items,
                "parent": parent,
            }
        )

    except PermissionError:
        return jsonify({"error": "Permission denied"}), 403
    except Exception as e:
        log.error(f"Error listing directory: {e}")
        return jsonify({"error": "Internal server error"}), 500


@file_browser.route("/load", methods=["POST"])
@require_auth
def load_file():
    """
    Load a file into ZnDraw.

    Request Body
    ------------
    path : str
        Relative path to the file.
    room : str, optional
        Room name. If not provided, generated from filename.
    start : int, optional
        Start frame for trajectory files.
    stop : int, optional
        Stop frame for trajectory files.
    step : int, optional
        Step frame for trajectory files.
    make_default : bool, optional
        Whether to make this the default room. Defaults to False.

    Returns
    -------
    Response
        JSON response with task information or error.
    """
    # Check if feature is enabled
    error = check_feature_enabled()
    if error:
        return jsonify(error[0]), error[1]

    # Get request data
    data = request.get_json()
    if not data or "path" not in data:
        return jsonify({"error": "Missing required field: path"}), 400

    requested_path = data["path"]
    config = current_app.extensions["config"]
    root = config.file_browser_root

    # Validate path
    target_path = validate_path(requested_path, root)
    if target_path is None:
        return jsonify({"error": "Invalid path or path outside root directory"}), 400

    # Check if file exists
    if not target_path.exists():
        return jsonify({"error": "File does not exist"}), 404

    # Check if path is a file
    if not target_path.is_file():
        return jsonify({"error": "Path is not a file"}), 400

    # Note: We don't check file type here - we let the backend reader handle it
    # Known formats use their assigned backend, unknown formats are tried via ASE

    # Check for duplicate file unless force_upload is true
    force_upload = data.get("force_upload", False)
    redis_client = current_app.extensions["redis"]

    if not force_upload:
        existing_room = _find_room_with_exact_file(target_path, redis_client, root)

        if existing_room:
            return jsonify(
                {
                    "status": "file_already_loaded",
                    "existingRoom": existing_room,
                    "message": f"This file is already loaded in room '{existing_room}'",
                    "filePath": str(target_path),
                    "options": {
                        "openExisting": f"Open room '{existing_room}'",
                        "createNew": "Create new room (reuse storage)",
                        "forceUpload": "Upload anyway (ignore existing)",
                    },
                }
            ), 200

    # Get or generate room name
    room = data.get("room")
    if not room:
        from zndraw.utils import generate_room_name

        # Use filename as base, not full path
        base_name = target_path.name
        room = generate_room_name(base_name, redis_client)

    # Get optional parameters
    start = data.get("start")
    stop = data.get("stop")
    step = data.get("step")
    make_default = data.get("make_default", False)

    # Queue file loading task
    from zndraw.app.tasks import read_file

    config = current_app.extensions["config"]
    server_url = config.server_url

    # Create description indicating file browser load
    description = f"{target_path.name} (loaded from file browser)"

    try:
        task = read_file.delay(
            file=str(target_path),
            room=room,
            server_url=server_url,
            start=start,
            stop=stop,
            step=step,
            make_default=make_default,
            root_path=root,
            description=description,
        )

        return jsonify(
            {
                "status": "queued",
                "room": room,
                "message": "File loading queued",
                "task_id": task.id,
            }
        )

    except Exception as e:
        log.error(f"Error queuing file load task: {e}")
        return jsonify({"error": "Failed to queue file loading task"}), 500


@file_browser.route("/upload", methods=["POST"])
def upload_file():
    """Upload file content for loading into ZnDraw.

    Note: This endpoint is always available (does not require file browser to be enabled)
    since drag/drop upload is a core feature.

    Form Data
    ---------
    file : FileStorage
        File to upload (multipart/form-data).
    room : str, optional
        Room name. If not provided, generated from filename.
    start : int, optional
        Start frame for trajectory files.
    stop : int, optional
        Stop frame for trajectory files.
    step : int, optional
        Step frame for trajectory files.
    make_default : bool, optional
        Whether to make this the default room. Defaults to False.

    Returns
    -------
    Response
        JSON response with task information or error.
    """
    import uuid

    from werkzeug.utils import secure_filename

    # Check if file was provided
    if "file" not in request.files:
        return jsonify({"error": "No file provided"}), 400

    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No file selected"}), 400

    # Validate file type
    ext = Path(file.filename).suffix.lstrip(".").lower()
    if ext not in FORMAT_BACKENDS:
        return jsonify({"error": f"Unsupported format: {ext}"}), 400

    # Create temporary directory for uploads
    config = current_app.extensions["config"]
    temp_dir = config.upload_temp
    os.makedirs(temp_dir, exist_ok=True)

    # Generate unique filename to avoid collisions
    unique_id = str(uuid.uuid4())
    temp_filename = f"{unique_id}_{secure_filename(file.filename)}"
    temp_path = os.path.join(temp_dir, temp_filename)

    # Save uploaded file
    try:
        file.save(temp_path)
    except Exception as e:
        log.error(f"Failed to save uploaded file: {e}")
        return jsonify({"error": "Failed to save uploaded file"}), 500

    # Get or generate room name
    room = request.form.get("room")
    if not room:
        from zndraw.utils import generate_room_name

        redis_client = current_app.extensions["redis"]
        room = generate_room_name(file.filename, redis_client)

    # Get optional parameters
    start = request.form.get("start", type=int)
    stop = request.form.get("stop", type=int)
    step = request.form.get("step", type=int)
    make_default = request.form.get("make_default", type=bool, default=False)

    # Queue Celery task with temp file path
    from zndraw.app.tasks import read_file

    config = current_app.extensions["config"]
    server_url = config.server_url

    # Create description indicating drag/drop upload
    description = f"{file.filename} (uploaded via drag & drop)"

    try:
        task = read_file.delay(
            file=temp_path,
            room=room,
            server_url=server_url,
            start=start,
            stop=stop,
            step=step,
            make_default=make_default,
            cleanup_after=True,
            description=description,
        )

        return jsonify(
            {
                "status": "queued",
                "room": room,
                "task_id": task.id,
                "message": "File uploaded and loading queued",
            }
        )

    except Exception as e:
        # Clean up temp file if task queuing failed
        try:
            os.remove(temp_path)
        except Exception:
            pass
        log.error(f"Error queuing file upload task: {e}")
        return jsonify({"error": "Failed to queue file loading task"}), 500


@file_browser.route("/create-room-from-file", methods=["POST"])
@require_auth
def create_room_from_existing_file():
    """Create a new room reusing physical storage from an existing room.

    Request Body
    ------------
    sourceRoom : str
        The room ID to copy storage from.
    newRoom : str, optional
        Name for the new room. If not provided, will be auto-generated.
    description : str, optional
        Description for the new room.

    Returns
    -------
    Response
        JSON response with new room information or error.

    Notes
    -----
    Creates a new room with identity mapping {i: i for i in range(frame_count_original)}
    WITHOUT re-uploading the file. Reuses the same physical storage backend.
    """
    # Check if feature is enabled
    error = check_feature_enabled()
    if error:
        return jsonify(error[0]), error[1]

    data = request.get_json()
    source_room = data.get("sourceRoom")
    new_room = data.get("newRoom")
    description = data.get("description")

    if not source_room:
        return jsonify({"error": "sourceRoom is required"}), 400

    redis_client = current_app.extensions["redis"]

    # Get source room metadata
    from zndraw.app.metadata_manager import RoomMetadataManager

    source_metadata_manager = RoomMetadataManager(redis_client, source_room)
    source_metadata = source_metadata_manager.get_all()

    if not source_metadata:
        return jsonify({"error": "Source room has no metadata"}), 400

    # Check if source room exists
    source_room_keys = RoomKeys(source_room)
    frame_count_int = redis_client.zcard(source_room_keys.trajectory_indices())

    if frame_count_int == 0:
        return jsonify({"error": "Source room not found or has no frames"}), 404

    # Generate new room name if not provided
    if not new_room:
        from pathlib import Path

        from zndraw.utils import generate_room_name

        # Use just the filename (not the full path) to generate consistent name
        file_path_str = source_metadata.get("relative_file_path") or source_room
        # Extract just the filename, not the directory path
        base_name = Path(file_path_str).name
        new_room = generate_room_name(base_name, redis_client=redis_client)

    # Create identity mapping for new room
    new_room_keys = RoomKeys(new_room)

    # Copy all frame references with identity mapping
    # Get all members from source room
    source_frames = redis_client.zrange(
        source_room_keys.trajectory_indices(), 0, -1, withscores=True
    )

    # Create mapping in new room
    for frame_key, score in source_frames:
        redis_client.zadd(new_room_keys.trajectory_indices(), {frame_key: score})

    # Copy metadata to new room
    new_metadata_manager = RoomMetadataManager(redis_client, new_room)
    new_metadata_manager.update(source_metadata)

    # Set room description and properties
    if description:
        redis_client.set(new_room_keys.description(), description)
    redis_client.set(new_room_keys.locked(), "0")
    redis_client.set(new_room_keys.hidden(), "0")

    log.info(
        f"Created room '{new_room}' from '{source_room}' with {frame_count_int} frames"
    )

    return jsonify(
        {
            "status": "success",
            "roomId": new_room,
            "sourceRoom": source_room,
            "frameCount": frame_count_int,
            "message": f"Room '{new_room}' created from '{source_room}' without re-uploading",
        }
    ), 201


@file_browser.route("/supported-types", methods=["GET"])
def supported_types():
    """
    Get list of supported file types.

    Note: This endpoint is always available (does not require file browser to be enabled)
    since it's needed for drag/drop upload validation which is a core feature.

    Returns
    -------
    Response
        JSON response with supported extensions, descriptions, and backend handlers.
    """
    # Generate extensions with dot prefix and descriptions from FORMAT_BACKENDS
    extensions = [f".{ext}" for ext in FORMAT_BACKENDS.keys()]
    descriptions = {
        f".{ext}": f"Supported by {', '.join(backends)}"
        for ext, backends in FORMAT_BACKENDS.items()
    }

    return jsonify(
        {
            "extensions": extensions,
            "descriptions": descriptions,
            "backends": FORMAT_BACKENDS,
        }
    )
