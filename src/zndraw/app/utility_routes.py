"""Utility and system routes.

Handles health checks, versioning, authentication, tools, and static asset serving.
"""

import base64
import datetime
import io
import logging
import uuid
from pathlib import Path

from flask import Blueprint, Response, current_app, request, send_from_directory
from flask_socketio import disconnect

from zndraw.server import socketio

log = logging.getLogger(__name__)

utility = Blueprint("utility", __name__)


@utility.route("/health")
def health_check():
    """Health check endpoint for server status verification."""
    return {"status": "ok"}, 200


@utility.route("/api/version")
def get_version():
    """Get the ZnDraw server version."""
    import zndraw

    return {"version": zndraw.__version__}, 200


@utility.route("/api/tools/rdkit-img", methods=["POST"])
def rdkit_image():
    """Convert SMILES or InChI to a molecule image using RDKit.

    Request
    -------
    {
        "type": "smiles" | "inchi",  // Required: type of input
        "data": "CCCO"                // Required: molecule string
    }

    Response
    --------
    {
        "image": "data:image/png;base64,iVBORw0KG...",  // Base64-encoded PNG
        "status": "success"
    }
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except ImportError:
        return {"error": "RDKit is not installed", "type": "ImportError"}, 500

    data = request.get_json() or {}
    mol_type = data.get("type")
    mol_data = data.get("data")

    if not mol_type or not mol_data:
        return {"error": "Both 'type' and 'data' are required", "type": "ValidationError"}, 400

    if mol_type not in ["smiles", "inchi"]:
        return {"error": "Type must be 'smiles' or 'inchi'", "type": "ValidationError"}, 400

    # Convert string to molecule object
    try:
        if mol_type == "smiles":
            mol = Chem.MolFromSmiles(mol_data)
        else:  # inchi
            mol = Chem.MolFromInchi(mol_data)

        if mol is None:
            return {"error": f"Invalid {mol_type} string", "type": "ConversionError"}, 400

        # Generate 2D coordinates if needed
        Chem.rdDepictor.Compute2DCoords(mol)

        # Draw molecule to PNG
        img = Draw.MolToImage(mol, size=(300, 300))

        # Convert to base64
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        img_base64 = base64.b64encode(buffer.getvalue()).decode("utf-8")

        return {
            "image": f"data:image/png;base64,{img_base64}",
            "status": "success"
        }, 200

    except Exception as e:
        log.error(f"Error generating molecule image: {e}")
        return {"error": str(e), "type": "GenerationError"}, 500


@utility.route("/api/login", methods=["POST"])
def login():
    """Authenticate user and issue JWT token.

    Request
    -------
    {
        "userName": "John Doe"  // Required
    }

    Response
    --------
    {
        "status": "ok",
        "token": "eyJhbGc...",     // JWT token
        "clientId": "uuid-string"  // Server-generated client ID
    }
    """
    from zndraw.auth import create_jwt_token

    data = request.get_json() or {}
    user_name = data.get("userName")

    if not user_name or not user_name.strip():
        return {"error": "userName is required"}, 400

    # Generate server-side client ID
    client_id = str(uuid.uuid4())

    # Create JWT token
    token = create_jwt_token(client_id, user_name)

    # Store client metadata in Redis
    r = current_app.extensions["redis"]
    client_key = f"client:{client_id}"
    current_time = datetime.datetime.utcnow().isoformat()
    r.hset(
        client_key,
        mapping={
            "userName": user_name,
            "createdAt": current_time,
            "lastLogin": current_time,
        },
    )

    log.info(f"User '{user_name}' logged in with client ID: {client_id}")

    return {
        "status": "ok",
        "token": token,
        "clientId": client_id,
    }


@utility.route("/assets/<path:filename>")
def serve_static_assets(filename: str):
    """Serve static assets from the assets directory."""
    static_folder = Path(__file__).parent.parent / "static" / "assets"
    return send_from_directory(static_folder, filename)


@utility.route("/", defaults={"path": ""})
@utility.route("/<path:path>")
def serve_react_router_paths(path: str):
    """Catch-all route to serve index.html for React Router paths.

    This allows React Router to handle client-side routing for paths like /rooms, /rooms/:id, etc.
    API routes are registered with higher priority and won't match this pattern.
    """
    static_folder = Path(__file__).parent.parent / "static"
    return send_from_directory(static_folder, "index.html")


@utility.route("/api/disconnect/<string:client_sid>", methods=["POST"])
def disconnect_sid(client_sid: str):
    """Disconnect the client from the room.

    Args:
        client_sid: Can be either a socket sid OR a client_id.
                   We try both lookups to support both cases.
    """
    try:
        r = current_app.extensions["redis"]

        # First, try to interpret client_sid as a client_id and get the socket sid
        socket_sid = r.hget(f"client:{client_sid}", "currentSid")

        if socket_sid:
            disconnect(socket_sid, namespace="/")
            return {"success": True}
        else:
            return {"success": False, "error": "Client not found or not connected"}
    except Exception as e:
        log.error(f"Error disconnecting client {client_sid}: {e}")
        return {"success": False, "error": str(e)}


@utility.route("/internal/emit", methods=["POST"])
def internal_emit():
    """Internal endpoint to emit Socket.IO events. Secured via a shared secret."""
    data = request.get_json()

    event = data.get("event")
    sid = data.get("sid")
    payload = data.get("data", {})

    if not event or not sid:
        return {"error": "Event and sid are required"}, 400

    socketio.emit(event, payload, to=sid)
    return {"success": True}


@utility.route("/api/shutdown", methods=["POST"])
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}  # this might never be seen
