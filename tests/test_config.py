"""Tests for ZnDrawConfig configuration module."""

import os

import pytest

from zndraw.config import ZnDrawConfig, get_config, reload_config


@pytest.fixture
def clean_env():
    """Save and restore environment variables after test."""
    original_env = os.environ.copy()
    yield
    # Restore original environment
    os.environ.clear()
    os.environ.update(original_env)
    # Force config reload after each test
    from zndraw import config as config_module

    config_module._config = None


def test_config_defaults(clean_env):
    """Test that config loads with default values when no env vars set."""
    # Clear all ZNDRAW_* env vars
    for key in list(os.environ.keys()):
        if key.startswith("ZNDRAW_") or key == "FLASK_SECRET_KEY":
            del os.environ[key]

    config = ZnDrawConfig()

    assert config.redis_url is None
    assert config.storage_path == "./zndraw-data"
    assert config.server_host == "localhost"
    assert config.server_port == 5000
    assert config.log_level == "WARNING"
    assert config.flask_secret_key == "dev-secret-key-change-in-production"
    assert config.admin_username is None
    assert config.admin_password is None
    assert config.upload_temp == "/tmp/zndraw_uploads"
    assert config.max_upload_mb == 500
    assert config.simgen_enabled is False
    assert config.file_browser_enabled is False
    assert config.celery_enabled is True


def test_config_loads_from_environment(clean_env):
    """Test that config loads values from environment variables."""
    os.environ["ZNDRAW_REDIS_URL"] = "redis://test:6379"
    os.environ["ZNDRAW_STORAGE_PATH"] = "/custom/storage"
    os.environ["ZNDRAW_SERVER_HOST"] = "example.com"
    os.environ["ZNDRAW_SERVER_PORT"] = "8000"
    os.environ["ZNDRAW_LOG_LEVEL"] = "DEBUG"
    os.environ["FLASK_SECRET_KEY"] = "test-secret-key"
    os.environ["ZNDRAW_SIMGEN_ENABLED"] = "true"
    os.environ["ZNDRAW_FILE_BROWSER_ENABLED"] = "1"
    os.environ["ZNDRAW_CELERY_ENABLED"] = "false"

    config = ZnDrawConfig()

    assert config.redis_url == "redis://test:6379"
    assert config.storage_path == "/custom/storage"
    assert config.server_host == "example.com"
    assert config.server_port == 8000
    assert config.log_level == "DEBUG"
    assert config.flask_secret_key == "test-secret-key"
    assert config.simgen_enabled is True
    assert config.file_browser_enabled is True
    assert config.celery_enabled is False


def test_config_boolean_parsing(clean_env):
    """Test that boolean environment variables are parsed correctly."""
    test_cases = [
        ("true", True),
        ("TRUE", True),
        ("1", True),
        ("yes", True),
        ("false", False),
        ("FALSE", False),
        ("0", False),
        ("no", False),
        ("", False),
        ("invalid", False),
    ]

    for value, expected in test_cases:
        os.environ["ZNDRAW_SIMGEN_ENABLED"] = value
        config = ZnDrawConfig()
        assert config.simgen_enabled == expected


def test_config_integer_parsing(clean_env):
    """Test that integer environment variables are parsed correctly."""
    os.environ["ZNDRAW_SERVER_PORT"] = "3000"
    os.environ["ZNDRAW_MAX_UPLOAD_MB"] = "1000"

    config = ZnDrawConfig()

    assert config.server_port == 3000
    assert config.max_upload_mb == 1000


def test_config_server_url_auto_generated(clean_env):
    """Test that server_url is auto-generated from host and port."""
    os.environ["ZNDRAW_SERVER_HOST"] = "example.com"
    os.environ["ZNDRAW_SERVER_PORT"] = "8080"

    config = ZnDrawConfig()

    assert config.server_url == "http://example.com:8080"


def test_config_server_url_explicit(clean_env):
    """Test that explicit server_url is used when provided."""
    os.environ["ZNDRAW_SERVER_URL"] = "https://custom.example.com"
    os.environ["ZNDRAW_SERVER_HOST"] = "localhost"
    os.environ["ZNDRAW_SERVER_PORT"] = "5000"

    config = ZnDrawConfig()

    assert config.server_url == "https://custom.example.com"


def test_config_server_url_docker_mode(clean_env):
    """Test that server_url uses localhost when host is 0.0.0.0 (Docker)."""
    os.environ["ZNDRAW_SERVER_HOST"] = "0.0.0.0"
    os.environ["ZNDRAW_SERVER_PORT"] = "5000"

    config = ZnDrawConfig()

    assert config.server_url == "http://localhost:5000"


def test_config_admin_validation_both_set(clean_env):
    """Test that admin credentials work when both username and password set."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    config = ZnDrawConfig()

    assert config.admin_username == "admin"
    assert config.admin_password == "secret"


def test_config_admin_validation_only_username(clean_env):
    """Test that only setting admin username raises ValueError."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    # Password not set

    with pytest.raises(
        ValueError,
        match="ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD must both be set or both be unset",
    ):
        ZnDrawConfig()


def test_config_admin_validation_only_password(clean_env):
    """Test that only setting admin password raises ValueError."""
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"
    # Username not set

    with pytest.raises(
        ValueError,
        match="ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD must both be set or both be unset",
    ):
        ZnDrawConfig()


def test_config_port_validation_invalid_low(clean_env):
    """Test that port number below 1 raises ValueError."""
    os.environ["ZNDRAW_SERVER_PORT"] = "0"

    with pytest.raises(ValueError, match="Invalid port number"):
        ZnDrawConfig()


def test_config_port_validation_invalid_high(clean_env):
    """Test that port number above 65535 raises ValueError."""
    os.environ["ZNDRAW_SERVER_PORT"] = "65536"

    with pytest.raises(ValueError, match="Invalid port number"):
        ZnDrawConfig()


def test_config_upload_size_validation(clean_env):
    """Test that upload size less than 1MB raises ValueError."""
    os.environ["ZNDRAW_MAX_UPLOAD_MB"] = "0"

    with pytest.raises(ValueError, match="Invalid max upload size"):
        ZnDrawConfig()


def test_config_log_level_validation_invalid(clean_env):
    """Test that invalid log level falls back to WARNING."""
    os.environ["ZNDRAW_LOG_LEVEL"] = "INVALID"

    config = ZnDrawConfig()

    # Should fall back to WARNING with a warning message
    assert config.log_level == "WARNING"


def test_config_log_level_validation_valid(clean_env):
    """Test that valid log levels are accepted."""
    valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

    for level in valid_levels:
        os.environ["ZNDRAW_LOG_LEVEL"] = level
        config = ZnDrawConfig()
        assert config.log_level == level


def test_get_config_singleton(clean_env):
    """Test that get_config returns singleton instance."""
    config1 = get_config()
    config2 = get_config()

    assert config1 is config2


def test_reload_config(clean_env):
    """Test that reload_config creates a new instance."""
    os.environ["ZNDRAW_STORAGE_PATH"] = "/old/path"

    config1 = get_config()
    assert config1.storage_path == "/old/path"

    # Change environment
    os.environ["ZNDRAW_STORAGE_PATH"] = "/new/path"

    # Reload config
    config2 = reload_config()

    assert config2.storage_path == "/new/path"
    assert config1 is not config2


def test_config_file_browser_root_default(clean_env):
    """Test that file_browser_root defaults to current working directory."""
    config = ZnDrawConfig()

    assert config.file_browser_root == os.getcwd()


def test_config_file_browser_root_custom(clean_env):
    """Test that custom file_browser_root is used when set."""
    os.environ["ZNDRAW_FILE_BROWSER_ROOT"] = "/custom/browser/root"

    config = ZnDrawConfig()

    assert config.file_browser_root == "/custom/browser/root"
