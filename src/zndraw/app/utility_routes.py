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

from zndraw.auth import get_current_user, require_admin, require_auth
from zndraw.server import socketio

from .redis_keys import SessionKeys, UserKeys

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


@utility.route("/api/config/global-settings")
def get_global_settings():
    """Get global settings for headbar features.

    Returns
    -------
    dict
        Global settings including feature flags for SiMGen and other features.
    """
    config = current_app.extensions["config"]
    return {
        "simgen": {
            "enabled": config.simgen_enabled,
        }
    }, 200


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
        return {
            "error": "Both 'type' and 'data' are required",
            "type": "ValidationError",
        }, 400

    if mol_type not in ["smiles", "inchi"]:
        return {
            "error": "Type must be 'smiles' or 'inchi'",
            "type": "ValidationError",
        }, 400

    # Convert string to molecule object
    try:
        if mol_type == "smiles":
            mol = Chem.MolFromSmiles(mol_data)
        else:  # inchi
            mol = Chem.MolFromInchi(mol_data)

        if mol is None:
            return {
                "error": f"Invalid {mol_type} string",
                "type": "ConversionError",
            }, 400

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
            "status": "success",
        }, 200

    except Exception as e:
        log.error(f"Error generating molecule image: {e}")
        return {"error": str(e), "type": "GenerationError"}, 500


@utility.route("/api/login", methods=["POST"])
def login():
    """Authenticate user and issue JWT token.

    Supports three modes:
    1. Guest login (no password): Creates new anonymous user
    2. User login (username + password): Authenticates existing registered user
    3. Admin login (admin credentials): Authenticates as admin

    Request
    -------
    {
        "userName": "JohnDoe",  // Optional, if not provided, generates guest username
        "password": "secret"    // Optional, if provided, authenticates existing user
    }

    Response
    --------
    {
        "status": "ok",
        "token": "eyJhbGc...",  // JWT token
        "userName": "JohnDoe",  // Actual username (may be generated)
        "role": "guest"         // User role: "guest", "user", or "admin"
    }
    """
    from zndraw.auth import create_jwt_token

    data = request.get_json() or {}
    user_name = data.get("userName")
    password = data.get("password")

    admin_service = current_app.extensions["admin_service"]
    user_service = current_app.extensions["user_service"]

    # CASE 1: Password provided - Authenticate existing user
    if password:
        if not user_name or not user_name.strip():
            return {"error": "userName is required when password is provided"}, 400

        user_name = user_name.strip()

        # Check if admin credentials
        is_admin_creds = admin_service.validate_admin_credentials(user_name, password)

        if is_admin_creds:
            # Admin login via env vars
            # Ensure admin user exists
            if not user_service.username_exists(user_name):
                user_service.create_user(user_name)
                user_service.register_user(user_name, user_name, password)

            # Grant admin privileges
            admin_service.grant_admin(user_name)
            user_service.update_last_login(user_name)

            token = create_jwt_token(user_name, role="admin")
            log.info(f"Admin '{user_name}' logged in")

            return {
                "status": "ok",
                "token": token,
                "userName": user_name,
                "role": "admin",
            }

        # Not admin credentials - check regular user
        if not user_service.username_exists(user_name):
            return {"error": "Invalid username or password"}, 401

        if not user_service.verify_password(user_name, password):
            return {"error": "Invalid username or password"}, 401

        # Valid user login
        user_service.update_last_login(user_name)
        is_admin = admin_service.is_admin(user_name)
        role = user_service.get_user_role(user_name)

        # In local mode, all users are admin; in deployment mode, check is_admin flag
        final_role = "admin" if is_admin else role.value

        token = create_jwt_token(user_name, role=final_role)
        log.info(f"User '{user_name}' logged in")

        return {
            "status": "ok",
            "token": token,
            "userName": user_name,
            "role": final_role,
        }

    # CASE 2: No password - Guest login or existing user re-auth
    if user_name and user_name.strip():
        user_name = user_name.strip()

        # Check if user exists
        if user_service.username_exists(user_name):
            # Existing user trying to login without password
            # Only allow if they're a guest (no password set)
            if user_service.is_registered(user_name):
                return {"error": "Password required for registered users"}, 401

            # Guest re-login
            user_service.update_last_login(user_name)
            is_admin = admin_service.is_admin(user_name)
            role = user_service.get_user_role(user_name)

            # In local mode, all users are admin
            final_role = "admin" if is_admin else role.value

            token = create_jwt_token(user_name, role=final_role)
            log.info(f"Guest '{user_name}' logged in")

            return {
                "status": "ok",
                "token": token,
                "userName": user_name,
                "role": final_role,
            }
        else:
            # New user with chosen username (still guest)
            try:
                user_service.create_user(user_name)
                is_admin = admin_service.is_admin(user_name)
                role = user_service.get_user_role(user_name)

                # In local mode, all users are admin
                final_role = "admin" if is_admin else role.value

                token = create_jwt_token(user_name, role=final_role)
                log.info(f"New guest '{user_name}' created")

                return {
                    "status": "ok",
                    "token": token,
                    "userName": user_name,
                    "role": final_role,
                }
            except ValueError as e:
                return {"error": str(e)}, 400

    # CASE 3: No username, no password - Generate anonymous guest
    # Generate unique guest username
    import secrets

    for _ in range(10):  # Try up to 10 times
        guest_name = f"user-{secrets.token_hex(4)}"
        if not user_service.username_exists(guest_name):
            user_service.create_user(guest_name)
            is_admin = admin_service.is_admin(guest_name)
            role = user_service.get_user_role(guest_name)

            # In local mode, all users are admin
            final_role = "admin" if is_admin else role.value

            token = create_jwt_token(guest_name, role=final_role)
            log.info(f"Anonymous guest '{guest_name}' created")

            return {
                "status": "ok",
                "token": token,
                "userName": guest_name,
                "role": final_role,
            }

    return {"error": "Failed to generate unique username"}, 500


@utility.route("/api/user/register", methods=["POST"])
@require_auth
def register_user():
    """Register a guest user with chosen username and password.

    Allows user to choose permanent username and set password.
    Old guest username will be deleted.

    Request
    -------
    {
        "userName": "MyUsername",  // Required, desired username
        "password": "mypassword"   // Required, no requirements
    }

    Response
    --------
    {
        "status": "ok",
        "token": "new-jwt-token",  // New token with new username
        "userName": "MyUsername",
        "role": "user"
    }
    """
    from zndraw.auth import create_jwt_token

    try:
        old_user_name = get_current_user()

        data = request.get_json() or {}
        new_user_name = data.get("userName")
        password = data.get("password")

        if not new_user_name:
            return {"error": "userName is required"}, 400

        if not password:
            return {"error": "password is required"}, 400

        user_service = current_app.extensions["user_service"]

        # Register the user (may change username)
        user_service.register_user(old_user_name, new_user_name.strip(), password)

        # Get role after registration (will be 'user' now)
        role = user_service.get_user_role(new_user_name.strip())

        # Create new JWT token with new username and role
        token = create_jwt_token(new_user_name.strip(), role=role.value)

        return {
            "status": "ok",
            "token": token,
            "userName": new_user_name.strip(),
            "role": "user",
        }

    except ValueError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except Exception as e:
        log.error(f"Error registering user: {e}")
        return {"error": "Registration failed", "type": "ServerError"}, 500


@utility.route("/api/user/change-password", methods=["POST"])
@require_auth
def change_password():
    """Change a user's password.

    User must provide current password for verification.

    Request
    -------
    {
        "oldPassword": "current",  // Required
        "newPassword": "new"       // Required
    }

    Response
    --------
    {
        "status": "ok"
    }
    """
    try:
        user_name = get_current_user()

        data = request.get_json() or {}
        old_password = data.get("oldPassword")
        new_password = data.get("newPassword")

        if not old_password or not new_password:
            return {"error": "oldPassword and newPassword are required"}, 400

        user_service = current_app.extensions["user_service"]
        user_service.change_password(user_name, old_password, new_password)

        return {"status": "ok"}

    except ValueError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except Exception as e:
        log.error(f"Error changing password: {e}")
        return {"error": "Password change failed", "type": "ServerError"}, 500


@utility.route("/api/user/role", methods=["GET"])
@require_auth
def get_user_role():
    """Get the current user's role.

    Response
    --------
    {
        "userName": "Name",
        "role": "guest" | "user" | "admin"
    }
    """
    try:
        user_name = get_current_user()

        user_service = current_app.extensions["user_service"]
        role = user_service.get_user_role(user_name)

        return {
            "userName": user_name,
            "role": role.value,
        }

    except Exception as e:
        log.error(f"Error getting user role: {e}")
        return {"error": "Failed to get user role", "type": "ServerError"}, 500


@utility.route("/api/rooms/<string:room_id>/settings", methods=["GET"])
@require_auth
def get_user_room_settings(room_id: str):
    """Get authenticated user's settings for a specific room.

    Settings are per-user, this endpoint returns the current user's
    settings scoped to the specified room.

    Parameters
    ----------
    room_id : str
        The room identifier

    Returns
    -------
    dict
        {
            "settings": {
                "category1": {...},
                "category2": {...}
            }
        }
    """
    try:
        user_name = get_current_user()  # From JWT
        settings_service = current_app.extensions["settings_service"]
        settings = settings_service.get_all(room_id, user_name)

        return {"settings": settings}, 200

    except Exception as e:
        log.error(f"Error getting user settings for room {room_id}: {e}")
        return {"error": "Failed to get user settings", "type": "ServerError"}, 500


@utility.route("/api/admin/users", methods=["GET"])
@require_admin
def list_users():
    """List all users (admin only).

    Response
    --------
    {
        "users": [
            {
                "userName": "Name",
                "role": "guest" | "user" | "admin",
                "createdAt": "2025-01-01T00:00:00",
                "lastLogin": "2025-01-01T00:00:00"
            },
            ...
        ]
    }
    """
    try:
        user_service = current_app.extensions["user_service"]
        users = user_service.list_all_users()

        return {"users": users}

    except Exception as e:
        log.error(f"Error listing users: {e}")
        return {"error": "Failed to list users", "type": "ServerError"}, 500


@utility.route("/api/admin/users/<string:user_name>/promote", methods=["POST"])
@require_admin
def promote_user(user_name: str):
    """Promote a user to admin (admin only).

    Response
    --------
    {
        "status": "ok",
        "role": "admin"
    }
    """
    try:
        admin_service = current_app.extensions["admin_service"]
        admin_service.grant_admin(user_name)

        return {
            "status": "ok",
            "role": "admin",
        }

    except Exception as e:
        log.error(f"Error promoting user: {e}")
        return {"error": "Failed to promote user", "type": "ServerError"}, 500


@utility.route("/api/admin/users/<string:user_name>/demote", methods=["POST"])
@require_admin
def demote_user(user_name: str):
    """Demote an admin to user (admin only).

    Response
    --------
    {
        "status": "ok",
        "role": "user"
    }
    """
    try:
        admin_service = current_app.extensions["admin_service"]
        admin_service.revoke_admin(user_name)

        return {
            "status": "ok",
            "role": "user",
        }

    except Exception as e:
        log.error(f"Error demoting user: {e}")
        return {"error": "Failed to demote user", "type": "ServerError"}, 500


@utility.route("/api/admin/users/<string:user_name>/reset-password", methods=["POST"])
@require_admin
def admin_reset_password(user_name: str):
    """Admin reset a user's password (admin only).

    Request
    -------
    {
        "newPassword": "newpass"  // Required
    }

    Response
    --------
    {
        "status": "ok"
    }
    """
    try:
        data = request.get_json() or {}
        new_password = data.get("newPassword")

        if not new_password:
            return {"error": "newPassword is required"}, 400

        user_service = current_app.extensions["user_service"]
        user_service.reset_password(user_name, new_password)

        return {"status": "ok"}

    except ValueError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except Exception as e:
        log.error(f"Error resetting password: {e}")
        return {"error": "Password reset failed", "type": "ServerError"}, 500


@utility.route("/api/admin/users/<string:user_name>", methods=["DELETE"])
@require_admin
def delete_user(user_name: str):
    """Delete a user (admin only).

    Hard deletion - removes all user data.

    Response
    --------
    {
        "status": "ok"
    }
    """
    try:
        user_service = current_app.extensions["user_service"]
        user_service.delete_user(user_name)

        return {"status": "ok"}

    except Exception as e:
        log.error(f"Error deleting user: {e}")
        return {"error": "Failed to delete user", "type": "ServerError"}, 500


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
        client_sid: Can be either a socket sid OR a userName.
                   We try both lookups to support both cases.
    """
    try:
        r = current_app.extensions["redis"]

        # First, try to interpret client_sid as a userName and get the socket sid
        user_keys = UserKeys(client_sid)
        socket_sid = r.hget(user_keys.hash_key(), "currentSid")

        if socket_sid:
            disconnect(socket_sid, namespace="/")
            return {"success": True}

        # Fallback: try to use client_sid directly as a socket sid
        # Check if this sid exists (has a userName mapping)
        session_keys = SessionKeys(client_sid)
        if r.exists(session_keys.username()):
            disconnect(client_sid, namespace="/")
            return {"success": True}

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
@require_admin
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}  # this might never be seen
