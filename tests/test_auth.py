"""Unit tests for JWT authentication utilities."""

import datetime

import jwt as pyjwt
import pytest

from zndraw.auth import (
    AuthError,
    create_jwt_token,
    decode_jwt_token,
    extract_token_from_request,
    get_current_user,
)


def test_create_jwt_token(app):
    """Test JWT token creation."""
    with app.app_context():
        token = create_jwt_token("TestUser", role="guest")
        assert token is not None
        assert isinstance(token, str)


def test_decode_jwt_token(app):
    """Test JWT token decoding."""
    with app.app_context():
        token = create_jwt_token("TestUser", role="guest")
        payload = decode_jwt_token(token)

        assert payload["sub"] == "TestUser"
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
            {"sub": "TestUser"},
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


def test_get_current_user_with_valid_token(app):
    """Test getting current user with valid JWT token."""
    with app.app_context():
        token = create_jwt_token("JohnDoe", role="guest")

        with app.test_request_context(headers={"Authorization": f"Bearer {token}"}):
            user_name = get_current_user()

            assert user_name == "JohnDoe"


def test_get_current_user_without_token_raises_error(app):
    """Test that missing token raises AuthError."""
    with app.app_context():
        with app.test_request_context():
            with pytest.raises(AuthError, match="No authentication token provided"):
                get_current_user()


def test_get_current_user_with_invalid_token_raises_error(app):
    """Test that invalid token raises AuthError."""
    with app.app_context():
        with app.test_request_context(
            headers={"Authorization": "Bearer invalid-token"}
        ):
            with pytest.raises(AuthError, match="Invalid token"):
                get_current_user()


def test_jwt_token_has_expiration(app):
    """Test that JWT tokens have a 7-day expiration claim."""
    with app.app_context():
        token = create_jwt_token("TestUser", role="guest")
        payload = decode_jwt_token(token)

        # Token should have exp claim
        assert "exp" in payload
        assert "iat" in payload

        # exp should be ~7 days after iat
        exp_dt = datetime.datetime.fromtimestamp(
            payload["exp"], tz=datetime.timezone.utc
        )
        iat_dt = datetime.datetime.fromtimestamp(
            payload["iat"], tz=datetime.timezone.utc
        )

        # Check that expiration is approximately 7 days (within 1 minute tolerance)
        expected_delta = datetime.timedelta(days=7)
        actual_delta = exp_dt - iat_dt
        assert abs((actual_delta - expected_delta).total_seconds()) < 60


def test_expired_jwt_token_raises_error(app):
    """Test that expired tokens raise AuthError."""
    from zndraw.utils.time import utc_now

    with app.app_context():
        # Create a token with past expiration
        secret_key = app.config["SECRET_KEY"]
        algorithm = app.config.get("JWT_ALGORITHM", "HS256")

        now = utc_now()
        expired_payload = {
            "sub": "TestUser",
            "role": "guest",
            "jti": "test-jti",
            "iat": now - datetime.timedelta(days=8),
            "exp": now - datetime.timedelta(days=1),
        }

        expired_token = pyjwt.encode(expired_payload, secret_key, algorithm=algorithm)

        with pytest.raises(AuthError, match="Token expired"):
            decode_jwt_token(expired_token)
