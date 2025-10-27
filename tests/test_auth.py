"""Unit tests for JWT authentication utilities."""

import jwt as pyjwt
import pytest

from zndraw.auth import (
    AuthError,
    create_jwt_token,
    decode_jwt_token,
    extract_token_from_request,
    get_current_client,
)


def test_create_jwt_token(app):
    """Test JWT token creation."""
    with app.app_context():
        token = create_jwt_token("client-123", "TestUser")
        assert token is not None
        assert isinstance(token, str)


def test_decode_jwt_token(app):
    """Test JWT token decoding."""
    with app.app_context():
        token = create_jwt_token("client-123", "TestUser")
        payload = decode_jwt_token(token)

        assert payload["sub"] == "client-123"
        assert payload["userName"] == "TestUser"
        assert "jti" in payload  # JWT ID present


def test_decode_invalid_token_raises_error(app):
    """Test that invalid tokens raise AuthError."""
    with app.app_context():
        with pytest.raises(AuthError, match="Invalid token"):
            decode_jwt_token("invalid.token.here")


def test_decode_token_with_wrong_secret_raises_error(app):
    """Test that tokens signed with wrong secret raise AuthError."""
    with app.app_context():
        # Create token with different secret
        wrong_token = pyjwt.encode(
            {"sub": "client-123", "userName": "TestUser"},
            "wrong-secret-key",
            algorithm="HS256",
        )

        with pytest.raises(AuthError, match="Invalid token"):
            decode_jwt_token(wrong_token)


def test_extract_token_from_request_with_bearer(app):
    """Test extracting JWT token from Authorization header."""
    with app.test_request_context(headers={"Authorization": "Bearer test-token-123"}):
        token = extract_token_from_request()
        assert token == "test-token-123"


def test_extract_token_from_request_without_bearer(app):
    """Test that non-Bearer auth headers return None."""
    with app.test_request_context(headers={"Authorization": "Basic xyz"}):
        token = extract_token_from_request()
        assert token is None


def test_extract_token_from_request_no_header(app):
    """Test that missing Authorization header returns None."""
    with app.test_request_context():
        token = extract_token_from_request()
        assert token is None


def test_get_current_client_with_valid_token(app):
    """Test getting current client with valid JWT token."""
    with app.app_context():
        token = create_jwt_token("client-456", "JohnDoe")

        with app.test_request_context(headers={"Authorization": f"Bearer {token}"}):
            client = get_current_client()

            assert client["clientId"] == "client-456"
            assert client["userName"] == "JohnDoe"


def test_get_current_client_without_token_raises_error(app):
    """Test that missing token raises AuthError."""
    with app.app_context():
        with app.test_request_context():
            with pytest.raises(AuthError, match="No authentication token provided"):
                get_current_client()


def test_get_current_client_with_invalid_token_raises_error(app):
    """Test that invalid token raises AuthError."""
    with app.app_context():
        with app.test_request_context(
            headers={"Authorization": "Bearer invalid-token"}
        ):
            with pytest.raises(AuthError, match="Invalid token"):
                get_current_client()
