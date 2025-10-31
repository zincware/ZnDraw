"""JWT authentication utilities for ZnDraw server."""

import logging
import typing as t
import uuid
from functools import wraps

import jwt
from flask import current_app, request

log = logging.getLogger(__name__)


class AuthError(Exception):
    """Authentication error.

    Parameters
    ----------
    message : str
        Error message
    status_code : int
        HTTP status code (default: 401)
    """

    def __init__(self, message: str, status_code: int = 401):
        self.message = message
        self.status_code = status_code
        super().__init__(self.message)


class AdminAccessError(Exception):
    """Admin access required error.

    Parameters
    ----------
    message : str
        Error message
    status_code : int
        HTTP status code (default: 403)
    """

    def __init__(self, message: str = "Admin access required", status_code: int = 403):
        self.message = message
        self.status_code = status_code
        super().__init__(self.message)


def create_jwt_token(client_id: str, user_name: str) -> str:
    """Create JWT token for authenticated client.

    Parameters
    ----------
    client_id : str
        Unique client identifier (UUID)
    user_name : str
        Display name of the user

    Returns
    -------
    str
        Encoded JWT token
    """
    secret_key = current_app.config["SECRET_KEY"]
    algorithm = current_app.config.get("JWT_ALGORITHM", "HS256")

    payload = {
        "sub": client_id,  # Subject: client ID
        "userName": user_name,  # Display name
        "jti": str(uuid.uuid4()),  # JWT ID for revocation (future use)
    }

    token = jwt.encode(payload, secret_key, algorithm=algorithm)
    log.info(f"Created JWT for client {client_id}")

    return token


def decode_jwt_token(token: str) -> dict:
    """Decode and validate JWT token.

    Parameters
    ----------
    token : str
        JWT token string

    Returns
    -------
    dict
        Decoded payload with claims

    Raises
    ------
    AuthError
        If token is invalid, expired, or malformed
    """
    secret_key = current_app.config["SECRET_KEY"]
    algorithm = current_app.config.get("JWT_ALGORITHM", "HS256")

    try:
        payload = jwt.decode(token, secret_key, algorithms=[algorithm])
        return payload
    except jwt.ExpiredSignatureError:
        raise AuthError("Token expired", 401)
    except jwt.InvalidTokenError:
        raise AuthError("Invalid token", 401)


def extract_token_from_request() -> str | None:
    """Extract JWT token from request Authorization header.

    Expects: Authorization: Bearer <token>

    Returns
    -------
    str | None
        Token string if found, None otherwise
    """
    auth_header = request.headers.get("Authorization")
    if auth_header and auth_header.startswith("Bearer "):
        return auth_header[7:]  # Remove "Bearer " prefix
    return None


def get_current_client() -> dict:
    """Get current authenticated client from request.

    Returns
    -------
    dict
        JWT payload with clientId and userName

    Raises
    ------
    AuthError
        If no token found or token is invalid
    """
    token = extract_token_from_request()
    if not token:
        raise AuthError("No authentication token provided", 401)

    payload = decode_jwt_token(token)
    return {
        "clientId": payload["sub"],
        "userName": payload["userName"],
    }


def require_auth(f):
    """Decorator to require JWT authentication for route.

    Usage
    -----
    @app.route("/api/protected")
    @require_auth
    def protected_route():
        client = get_current_client()
        return {"clientId": client["clientId"]}
    """

    @wraps(f)
    def decorated_function(*args, **kwargs):
        try:
            get_current_client()  # Validate token
            return f(*args, **kwargs)
        except AuthError as e:
            return {"error": e.message}, e.status_code

    return decorated_function


def require_admin(f):
    """Decorator to require admin privileges for route.

    This decorator combines authentication and admin authorization.
    First validates JWT token, then checks admin status via AdminService.

    In local mode (no admin credentials configured), all authenticated users
    are considered admins. In deployment mode (admin credentials set), only
    users who logged in with correct admin credentials have access.

    Usage
    -----
    @app.route("/api/admin-only")
    @require_admin
    def admin_route():
        # Only admins can access this
        return {"status": "ok"}

    Raises
    ------
    Returns 401 if not authenticated
    Returns 403 if authenticated but not admin
    """

    @wraps(f)
    def decorated_function(*args, **kwargs):
        from flask import current_app

        try:
            # First, validate authentication
            client = get_current_client()
            client_id = client["clientId"]

            # Then, check admin status
            admin_service = current_app.extensions.get("admin_service")
            if not admin_service:
                log.error("AdminService not initialized")
                return {"error": "Server configuration error", "type": "ServerError"}, 500

            if not admin_service.is_admin(client_id):
                raise AdminAccessError()

            return f(*args, **kwargs)

        except AuthError as e:
            return {"error": e.message, "type": "AuthError"}, e.status_code
        except AdminAccessError as e:
            return {"error": e.message, "type": "AdminAccessError"}, e.status_code

    return decorated_function
