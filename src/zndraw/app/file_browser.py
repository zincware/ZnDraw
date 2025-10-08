"""File browser endpoints for browsing and loading local files."""

import logging
from datetime import datetime, timezone
from pathlib import Path

from flask import Blueprint, current_app, jsonify, request

from zndraw.cli import path_to_room

log = logging.getLogger(__name__)

file_browser = Blueprint("file_browser", __name__, url_prefix="/api/file-browser")

# Supported file extensions for molecular structure files
SUPPORTED_EXTENSIONS = {
    ".xyz": "XYZ coordinate file",
    ".pdb": "Protein Data Bank file",
    ".h5md": "H5MD molecular data",
    ".h5": "H5MD molecular data",
    ".hdf5": "HDF5 file",
    ".traj": "Trajectory file",
    ".nc": "NetCDF trajectory file",
    ".extxyz": "Extended XYZ file",
    ".cif": "Crystallographic Information File",
    ".mol": "MDL Molfile",
    ".sdf": "Structure Data File",
}


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
    Check if file is supported by ASE or znh5md.

    Parameters
    ----------
    filepath : Path
        The file path to check.

    Returns
    -------
    bool
        True if file extension is supported.
    """
    return filepath.suffix.lower() in SUPPORTED_EXTENSIONS


def check_feature_enabled() -> tuple[dict[str, str], int] | None:
    """
    Check if file browser feature is enabled.

    Returns
    -------
    tuple[dict[str, str], int] | None
        Error response if disabled, None if enabled.
    """
    if not current_app.config.get("FILE_BROWSER_ENABLED", False):
        return {"error": "File browser feature is not enabled"}, 403
    return None


@file_browser.route("/list", methods=["GET"])
def list_directory():
    """
    List contents of a directory.

    Query Parameters
    ----------------
    path : str, optional
        Relative path from browser root. Defaults to root.

    Returns
    -------
    Response
        JSON response with directory contents or error.
    """
    # Check if feature is enabled
    error = check_feature_enabled()
    if error:
        return jsonify(error[0]), error[1]

    # Get requested path
    requested_path = request.args.get("path", "")
    root = current_app.config.get("FILE_BROWSER_ROOT", ".")

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

    # List directory contents
    try:
        items = []
        for item in sorted(target_path.iterdir(), key=lambda x: (not x.is_dir(), x.name)):
            # Skip hidden files
            if item.name.startswith("."):
                continue

            item_info = {
                "name": item.name,
                "type": "directory" if item.is_dir() else "file",
                "size": item.stat().st_size if item.is_file() else None,
                "modified": datetime.fromtimestamp(
                    item.stat().st_mtime, tz=timezone.utc
                ).isoformat(),
            }

            # Mark if file is supported
            if item.is_file():
                item_info["supported"] = is_supported_file(item)

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
    root = current_app.config.get("FILE_BROWSER_ROOT", ".")

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

    # Check if file type is supported
    if not is_supported_file(target_path):
        return (
            jsonify(
                {
                    "error": f"Unsupported file type: {target_path.suffix}",
                    "supported_extensions": list(SUPPORTED_EXTENSIONS.keys()),
                }
            ),
            400,
        )

    # Get or generate room name
    room = data.get("room")
    if not room:
        room = path_to_room(str(target_path))

    # Get optional parameters
    start = data.get("start")
    stop = data.get("stop")
    step = data.get("step")
    make_default = data.get("make_default", False)

    # Queue file loading task
    from zndraw.app.tasks import read_file

    server_url = current_app.config.get("SERVER_URL", "http://localhost:5000")

    try:
        task = read_file.delay(
            file=str(target_path),
            room=room,
            server_url=server_url,
            start=start,
            stop=stop,
            step=step,
            make_default=make_default,
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


@file_browser.route("/supported-types", methods=["GET"])
def supported_types():
    """
    Get list of supported file types.

    Returns
    -------
    Response
        JSON response with supported extensions and descriptions.
    """
    # Check if feature is enabled
    error = check_feature_enabled()
    if error:
        return jsonify(error[0]), error[1]

    return jsonify(
        {
            "extensions": list(SUPPORTED_EXTENSIONS.keys()),
            "descriptions": SUPPORTED_EXTENSIONS,
        }
    )
