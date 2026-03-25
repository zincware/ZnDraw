# tests/test_dependencies.py
from unittest.mock import MagicMock

from fastapi import FastAPI

from zndraw_joblib.dependencies import (
    get_internal_registry,
    get_joblib_settings,
    get_tsio,
)
from zndraw_joblib.settings import JobLibSettings


def _make_request(app: FastAPI) -> MagicMock:
    """Create a mock Request with the given app."""
    request = MagicMock()
    request.app = app
    return request


def test_get_joblib_settings_returns_settings_from_state():
    """get_joblib_settings reads from app.state.joblib_settings."""
    app = FastAPI()
    app.state.joblib_settings = JobLibSettings()
    request = _make_request(app)
    settings = get_joblib_settings(request)
    assert isinstance(settings, JobLibSettings)


def test_get_joblib_settings_returns_same_instance():
    """get_joblib_settings returns the exact instance from app.state."""
    app = FastAPI()
    app.state.joblib_settings = JobLibSettings()
    request = _make_request(app)
    settings1 = get_joblib_settings(request)
    settings2 = get_joblib_settings(request)
    assert settings1 is settings2


def test_get_internal_registry_import():
    assert callable(get_internal_registry)


def test_get_tsio_returns_none_when_not_set():
    """get_tsio returns None when tsio is not on app.state."""
    app = FastAPI()
    request = _make_request(app)
    result = get_tsio(request)
    assert result is None


def test_get_tsio_returns_tsio_from_state():
    """get_tsio returns tsio from app.state when set."""
    app = FastAPI()
    sentinel = object()
    app.state.tsio = sentinel
    request = _make_request(app)
    result = get_tsio(request)
    assert result is sentinel
