"""Utility and system routes.

Handles health checks, versioning, authentication, tools, and static asset serving.
"""

import base64
import logging
import secrets
from pathlib import Path

from flask import Blueprint, current_app, request, send_from_directory
from flask_socketio import disconnect

from zndraw.auth import (
    AuthError,
    create_jwt_token,
    get_current_user,
    require_admin,
    require_auth,
)
from zndraw.server import socketio
from zndraw.services.user_service import (
    PasswordValidationError,
    UsernameValidationError,
    validate_username,
)

from .redis_keys import SessionKeys, UserKeys

log = logging.getLogger(__name__)

# Generic authentication error message to prevent username enumeration
AUTH_FAILED_MSG = "Authentication failed"

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
        "data": "CCCO",              // Required: molecule string
        "dark": true                 // Optional: dark mode (default: false)
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
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return {"error": "RDKit is not installed", "type": "ImportError"}, 500

    data = request.get_json() or {}
    mol_type = data.get("type")
    mol_data = data.get("data")
    dark_mode = data.get("dark", False)

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

        # Draw molecule to PNG using MolDraw2DCairo for full control
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
        opts = drawer.drawOptions()

        # Transparent background - container handles visual background
        opts.clearBackground = False

        if dark_mode:
            # Dark mode: light bonds/atoms for visibility on dark containers
            opts.setSymbolColour((0.9, 0.9, 0.9))
            opts.updateAtomPalette({6: (0.9, 0.9, 0.9)})

        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        png_data = drawer.GetDrawingText()

        # Convert to base64
        img_base64 = base64.b64encode(png_data).decode("utf-8")

        return {
            "image": f"data:image/png;base64,{img_base64}",
            "status": "success",
        }, 200

    except Exception as e:
        log.error(f"Error generating molecule image: {e}")
        return {"error": str(e), "type": "GenerationError"}, 500


@utility.route("/api/login", methods=["POST"])
def login():
    """Authenticate user and issue JWT. NEVER creates users.

    User must exist (from /api/user/register) before login.
    Guest users: no password required.
    Registered users: password validation against Redis.
    Admin users: validate against configured credentials.

    Request
    -------
    {
        "userName": "JohnDoe",  // Required
        "password": "secret"    // Optional, required for registered users
    }

    Response
    --------
    {
        "status": "ok",
        "token": "eyJhbGc...",
        "userName": "JohnDoe",
        "role": "guest"
    }
    """
    data = request.get_json() or {}
    user_name = (data.get("userName") or "").strip()
    password = data.get("password")

    if not user_name:
        return {"error": "userName required"}, 400

    # Validate username format
    try:
        validate_username(user_name)
    except UsernameValidationError as e:
        return {"error": str(e), "type": "ValidationError"}, 400

    admin_service = current_app.extensions["admin_service"]
    user_service = current_app.extensions["user_service"]

    # Password provided - authenticate existing user or admin
    if password:
        return _handle_password_login(user_name, password, user_service, admin_service)

    # No password - guest login (user must exist from /api/user/register)
    return _handle_guest_login(user_name, user_service, admin_service)


def _handle_password_login(user_name, password, user_service, admin_service):
    """Handle login with password (registered user or admin)."""
    # Check admin credentials first
    if admin_service.validate_admin_credentials(user_name, password):
        # Admin login - ensure user exists in Redis for tracking
        user_service.ensure_user_exists(user_name)
        admin_service.grant_admin(user_name)
        user_service.update_last_login(user_name)
        return _issue_token(user_name, "admin")

    # If username matches configured admin but password is wrong
    if admin_service.is_admin_username(user_name):
        return {"error": AUTH_FAILED_MSG}, 401

    # User must exist (from /api/user/register)
    if not user_service.username_exists(user_name):
        return {"error": AUTH_FAILED_MSG}, 401

    if not user_service.verify_password(user_name, password):
        return {"error": AUTH_FAILED_MSG}, 401

    user_service.update_last_login(user_name)
    role = _determine_role(user_name, user_service, admin_service)
    return _issue_token(user_name, role)


def _handle_guest_login(user_name, user_service, admin_service):
    """Handle guest login (no password) - user must exist from /api/user/register."""
    # User must exist (from /api/user/register)
    if not user_service.username_exists(user_name):
        return {"error": AUTH_FAILED_MSG}, 401

    # If existing registered user, require password
    if user_service.is_registered(user_name):
        return {"error": AUTH_FAILED_MSG}, 401

    user_service.update_last_login(user_name)
    role = _determine_role(user_name, user_service, admin_service)
    log.debug(f"Issued JWT for guest '{user_name}'")
    return _issue_token(user_name, role)


def _determine_role(user_name, user_service, admin_service):
    """Determine user role based on admin status."""
    if admin_service.is_admin(user_name):
        return "admin"
    if user_service.username_exists(user_name):
        return user_service.get_user_role(user_name).value
    # In local mode (no admin configured), all users are admin
    if not admin_service.is_deployment_mode():
        return "admin"
    return "guest"


def _issue_token(user_name, role):
    """Issue JWT token response."""
    token = create_jwt_token(user_name, role=role)
    return {"status": "ok", "token": token, "userName": user_name, "role": role}


@utility.route("/api/user/register", methods=["POST"])
def register_user():
    """Register a new user (guest or with password).

    This is the ONLY endpoint that creates users. Call this before /api/login.
    - No password: Creates guest user (can be upgraded later)
    - With password: Creates registered user with credentials

    Request
    -------
    {
        "userName": "MyUsername",  // Optional, generates random if not provided
        "password": "mypassword"   // Optional, creates guest if not provided
    }

    Response
    --------
    {
        "status": "ok",
        "userName": "MyUsername"
    }
    """
    data = request.get_json() or {}
    user_name = (data.get("userName") or "").strip()
    password = data.get("password")

    user_service = current_app.extensions["user_service"]

    # Generate username if not provided (for guests)
    if not user_name:
        user_name = f"user-{secrets.token_hex(4)}"

    # Validate username format
    try:
        validate_username(user_name)
    except UsernameValidationError as e:
        return {"error": str(e), "type": "ValidationError"}, 400

    # Check if user already exists
    if user_service.username_exists(user_name):
        return {"error": f"Username '{user_name}' already exists"}, 409

    try:
        if password:
            # Register with password (full registration)
            user_service.register_user(user_name, user_name, password)
            log.debug(f"Registered new user '{user_name}' with password")
        else:
            # Guest registration (no password)
            user_service.ensure_user_exists(user_name)
            log.debug(f"Created guest user '{user_name}'")

        return {"status": "ok", "userName": user_name}, 201

    except PasswordValidationError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except ValueError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except Exception:
        log.exception("Error registering user")
        return {"error": "Registration failed", "type": "ServerError"}, 500


@utility.route("/api/user/upgrade", methods=["POST"])
@require_auth
def upgrade_user():
    """Upgrade a guest user to registered with chosen username and password.

    Allows guest user to choose permanent username and set password.
    Old guest username will be deleted if different.

    Request
    -------
    {
        "userName": "MyUsername",  // Required, desired username
        "password": "mypassword"   // Required, must meet password requirements
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

    except PasswordValidationError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except ValueError as e:
        return {"error": str(e), "type": "ValidationError"}, 400
    except Exception:
        log.exception("Error upgrading user")
        return {"error": "Upgrade failed", "type": "ServerError"}, 500


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
def exit_app():
    """Endpoint to gracefully shut down the server.

    Requires either:
    1. Shutdown token from PID file (for CLI shutdown from same user)
    2. Admin authentication (for admin users via UI)

    Security
    --------
    - Shutdown token is generated at server start and stored in PID file
    - Only the user who started the server has access to the PID file
    - Token is sent in X-Shutdown-Token header
    - Fallback to admin authentication for UI-based shutdown
    """
    # Check for shutdown token (from CLI)
    shutdown_token = request.headers.get("X-Shutdown-Token")
    expected_token = current_app.config.get("SHUTDOWN_TOKEN")

    if shutdown_token and expected_token and shutdown_token == expected_token:
        # Valid shutdown token from CLI
        log.debug("Server shutdown via CLI with valid token")
        socketio.stop()
        return {"success": True}

    # Fallback to admin authentication check
    try:
        user_name = get_current_user()
        admin_service = current_app.extensions["admin_service"]

        if admin_service.is_admin(user_name):
            log.debug(f"Server shutdown by admin user: {user_name}")
            socketio.stop()
            return {"success": True}
        else:
            # User is authenticated but not an admin
            log.warning(f"Non-admin user '{user_name}' attempted shutdown")
            return {
                "error": "Admin access required for shutdown",
                "type": "AdminAccessError",
            }, 403
    except AuthError as e:
        # Not authenticated
        log.warning("Unauthorized shutdown attempt - no valid token")
        return {"error": e.message, "type": "AuthError"}, 401

    log.warning("Unauthorized shutdown attempt")
    return {
        "error": "Unauthorized - admin access or valid shutdown token required",
        "type": "AuthError",
    }, 401
