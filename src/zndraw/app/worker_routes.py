"""Worker management routes.

Handles worker registration via REST endpoints.
"""

import json
import logging

from flask import Blueprint, current_app, request

from zndraw.server import socketio

from .constants import SocketEvents
from .redis_keys import ExtensionKeys, FilesystemKeys, SessionKeys

log = logging.getLogger(__name__)

workers = Blueprint("workers", __name__)


@workers.route("/api/workers/register", methods=["POST"])
def register_worker():
    """Register a worker extension via REST.

    Request Payload
    ---------------
    {
        "sessionId": str,
        "roomId": str,
        "name": str,
        "category": str,
        "schema": dict,
        "public": bool (optional, defaults to False)
    }

    Returns
    -------
    200 OK:
        {
            "success": true,
            "workerId": str,
            "message": str
        }
    400 Bad Request:
        {
            "success": false,
            "error": str,
            "code": "MISSING_FIELD"
        }
    403 Forbidden:
        {
            "success": false,
            "error": str,
            "code": "ADMIN_REQUIRED"
        }
    409 Conflict:
        {
            "success": false,
            "error": str,
            "code": "SCHEMA_MISMATCH"
        }
    423 Locked:
        {
            "success": false,
            "error": str,
            "code": "RESERVED_NAME"
        }
    """
    from zndraw.auth import AuthError, get_current_user

    r = current_app.extensions["redis"]

    # Authenticate and get user from JWT token
    try:
        user_name = get_current_user()
    except AuthError as e:
        return {"success": False, "error": e.message}, e.status_code

    # Parse request payload
    data = request.json
    if not data:
        return {
            "success": False,
            "error": "Request body required",
            "code": "MISSING_FIELD",
        }, 400

    try:
        session_id = data["sessionId"]
        room_id = data["roomId"]
        name = data["name"]
        category = data["category"]
        schema = data["schema"]
        public = data.get("public", False)  # Default to False if not specified
    except KeyError as e:
        return {
            "success": False,
            "error": f"Missing required field: {e}",
            "code": "MISSING_FIELD",
        }, 400

    # Resolve sessionId to socket sid (worker_id)
    sid = r.get(SessionKeys.session_to_sid(session_id))
    if not sid:
        return {
            "success": False,
            "error": "Session not found. Client must connect via socket before registering.",
            "code": "SESSION_NOT_FOUND",
        }, 400

    worker_id = sid

    # Check admin privileges for public extensions
    if public:
        try:
            # Get role from Redis (stored during socket connection)
            session_keys = SessionKeys(sid)
            role = r.get(session_keys.role())
            if role is None:
                log.error(f"Role not found for sid {sid}")
                return {
                    "success": False,
                    "error": "Failed to verify admin privileges",
                    "code": "AUTH_ERROR",
                }, 500

            if role != "admin":
                log.warning(
                    f"User {user_name} (role: {role}) attempted to register public extension '{name}'"
                )
                return {
                    "success": False,
                    "error": "Only admin users can register public extensions",
                    "code": "ADMIN_REQUIRED",
                }, 403
        except Exception as e:
            log.error(f"Failed to check user role: {e}")
            return {
                "success": False,
                "error": "Failed to verify admin privileges",
                "code": "AUTH_ERROR",
            }, 500

    # Prevent overriding server-side extensions (only for global/public extensions)
    # Room-scoped extensions are allowed to have same names as server-side extensions
    if public:
        from zndraw.extensions.analysis import analysis
        from zndraw.extensions.modifiers import modifiers
        from zndraw.extensions.selections import selections
        from zndraw.settings import settings

        category_map = {
            "selections": selections,
            "modifiers": modifiers,
            "settings": settings,
            "analysis": analysis,
        }

        if category in category_map and name in category_map[category]:
            log.warning(
                f"Blocked attempt to register global extension '{name}' in category '{category}' "
                f"- name conflicts with server-side extension (security violation)"
            )
            return {
                "success": False,
                "error": f"Cannot register extension '{name}': name is reserved for server-side extensions",
                "code": "RESERVED_NAME",
            }, 423

    scope = "global" if public else f"room {room_id}"
    log.info(
        f"Registering {'global' if public else 'room-scoped'} extension: name={name}, category={category}, worker_id={worker_id}, scope={scope}"
    )

    # Use global or room-scoped keys based on public flag
    if public:
        keys = ExtensionKeys.for_global_extension(category, name)
        worker_extensions_key = ExtensionKeys.global_user_extensions_key(
            category, worker_id
        )
    else:
        keys = ExtensionKeys.for_extension(room_id, category, name)
        worker_extensions_key = ExtensionKeys.user_extensions_key(
            room_id, category, worker_id
        )

    existing_schema = r.hget(keys.schema, name)

    if existing_schema is not None:
        existing_schema = json.loads(existing_schema)
        if existing_schema != schema:
            return {
                "success": False,
                "error": "Extension with this name already exists with a different schema",
                "code": "SCHEMA_MISMATCH",
            }, 409

        # Re-registration with same schema
        from datetime import datetime

        # Add to extension registry with timestamp
        r.hset(keys.workers, worker_id, datetime.utcnow().timestamp())

        # Initialize worker capacity if not exists
        capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
        if not r.exists(capacity_key):
            r.set(capacity_key, 1)  # Default: can handle 1 job

        r.sadd(worker_extensions_key, name)

        log.info(
            f"Worker {worker_id} re-registered for extension '{name}' "
            f"in category '{category}', invalidating schema"
        )

        # Assign any pending jobs to this newly idle worker
        from .job_dispatcher import assign_pending_jobs_for_extension

        ext_room_id = None if public else room_id
        assigned = assign_pending_jobs_for_extension(
            r, socketio, ext_room_id, category, name, worker_id
        )
        if assigned > 0:
            log.info(
                f"Assigned {assigned} pending job(s) to re-registered worker {worker_id}"
            )

        # Emit schema invalidation
        if public:
            # Broadcast to all connected clients (global extension)
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
            )
        else:
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_id}",
            )

        return {
            "success": True,
            "workerId": worker_id,
            "message": "Extension already registered with same schema. Worker marked as idle.",
        }, 200
    else:
        # Brand new extension
        from datetime import datetime

        with r.pipeline() as pipe:
            pipe.hset(keys.schema, name, json.dumps(schema))
            # Add to extension registry with timestamp
            pipe.hset(keys.workers, worker_id, datetime.utcnow().timestamp())
            pipe.sadd(worker_extensions_key, name)
            pipe.execute()

        # Initialize worker capacity if not exists (outside pipeline to check existence)
        capacity_key = ExtensionKeys.worker_capacity_key(worker_id)
        if not r.exists(capacity_key):
            r.set(capacity_key, 1)  # Default: can handle 1 job

        # Emit schema invalidation
        if public:
            # Broadcast to all connected clients (global extension)
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
            )
        else:
            socketio.emit(
                SocketEvents.INVALIDATE_SCHEMA,
                {"category": category},
                to=f"room:{room_id}",
            )

        log.info(f"Successfully registered new extension '{name}' for worker {worker_id}")

        # Assign any pending jobs to this newly idle worker
        from .job_dispatcher import assign_pending_jobs_for_extension

        ext_room_id = None if public else room_id
        assigned = assign_pending_jobs_for_extension(
            r, socketio, ext_room_id, category, name, worker_id
        )
        if assigned > 0:
            log.info(
                f"Assigned {assigned} pending job(s) to newly registered worker {worker_id}"
            )

        return {
            "success": True,
            "workerId": worker_id,
            "message": "Extension registered successfully",
        }, 200


@workers.route("/api/workers/filesystem/register", methods=["POST"])
def register_filesystem():
    """Register a worker filesystem via REST.

    Request Payload
    ---------------
    {
        "sessionId": str,
        "roomId": str,
        "name": str,
        "fsType": str,
        "public": bool (optional, defaults to False)
    }

    Returns
    -------
    200 OK:
        {
            "success": true,
            "workerId": str,
            "message": str
        }
    400 Bad Request:
        {
            "success": false,
            "error": str,
            "code": "MISSING_FIELD"
        }
    403 Forbidden:
        {
            "success": false,
            "error": str,
            "code": "ADMIN_REQUIRED"
        }
    409 Conflict:
        {
            "success": false,
            "error": str,
            "code": "ALREADY_REGISTERED"
        }
    """
    from zndraw.auth import AuthError, get_current_user

    r = current_app.extensions["redis"]

    # Authenticate and get user from JWT token
    try:
        user_name = get_current_user()
    except AuthError as e:
        return {"success": False, "error": e.message}, e.status_code

    # Parse request payload
    data = request.json
    if not data:
        return {
            "success": False,
            "error": "Request body required",
            "code": "MISSING_FIELD",
        }, 400

    try:
        session_id = data["sessionId"]
        room_id = data["roomId"]
        name = data["name"]
        fs_type = data["fsType"]
        public = data.get("public", False)
    except KeyError as e:
        return {
            "success": False,
            "error": f"Missing required field: {e}",
            "code": "MISSING_FIELD",
        }, 400

    # Resolve sessionId to socket sid (worker_id)
    sid = r.get(SessionKeys.session_to_sid(session_id))
    if not sid:
        return {
            "success": False,
            "error": "Session not found. Client must connect via socket before registering.",
            "code": "SESSION_NOT_FOUND",
        }, 400

    worker_id = sid

    # Check admin privileges for public filesystems
    if public:
        try:
            session_keys = SessionKeys(sid)
            role = r.get(session_keys.role())
            if role is None:
                log.error(f"Role not found for sid {sid}")
                return {
                    "success": False,
                    "error": "Failed to verify admin privileges",
                    "code": "AUTH_ERROR",
                }, 500

            if role != "admin":
                log.warning(
                    f"User {user_name} (role: {role}) attempted to register public filesystem '{name}'"
                )
                return {
                    "success": False,
                    "error": "Only admin users can register public filesystems",
                    "code": "ADMIN_REQUIRED",
                }, 403
        except Exception as e:
            log.error(f"Failed to check user role: {e}")
            return {
                "success": False,
                "error": "Failed to verify admin privileges",
                "code": "AUTH_ERROR",
            }, 500

    scope = "global" if public else f"room {room_id}"
    log.info(
        f"Registering {'global' if public else 'room-scoped'} filesystem: name={name}, type={fs_type}, worker_id={worker_id}, scope={scope}"
    )

    # Use global or room-scoped keys based on public flag
    if public:
        keys = FilesystemKeys.for_global_filesystem(name)
        worker_filesystems_key = FilesystemKeys.global_user_filesystems_key(worker_id)
    else:
        keys = FilesystemKeys.for_filesystem(room_id, name)
        worker_filesystems_key = FilesystemKeys.user_filesystems_key(room_id, worker_id)

    # Check if filesystem already exists
    existing_worker = r.get(keys.worker)

    if existing_worker is not None:
        if existing_worker != worker_id:
            # Different worker already registered this filesystem
            return {
                "success": False,
                "error": f"Filesystem '{name}' is already registered by another worker",
                "code": "ALREADY_REGISTERED",
            }, 409

        # Same worker re-registering - update metadata
        log.info(
            f"Worker {worker_id} re-registered filesystem '{name}', updating metadata"
        )

    # Store filesystem metadata
    fs_metadata = {
        "name": name,
        "fsType": fs_type,
        "sessionId": session_id,
        "userName": user_name,
        "public": "true" if public else "false",
    }

    with r.pipeline() as pipe:
        # Store filesystem metadata as a hash
        pipe.hset(keys.metadata, mapping=fs_metadata)
        # Store worker ID
        pipe.set(keys.worker, worker_id)
        # Add filesystem name to worker's filesystem set
        pipe.sadd(worker_filesystems_key, name)
        pipe.execute()

    # Emit filesystem update event
    if public:
        # Global filesystem - notify all clients
        socketio.emit(SocketEvents.FILESYSTEMS_UPDATE, {"scope": "global"})
    else:
        # Room-scoped filesystem - notify clients in that room
        socketio.emit(
            SocketEvents.FILESYSTEMS_UPDATE,
            {"scope": "room"},
            to=f"room:{room_id}",
        )

    message = (
        "Filesystem registered successfully"
        if existing_worker is None
        else "Filesystem metadata updated successfully"
    )
    log.info(f"Filesystem '{name}' registered successfully by worker {worker_id}")

    return {
        "success": True,
        "workerId": worker_id,
        "message": message,
    }, 200
