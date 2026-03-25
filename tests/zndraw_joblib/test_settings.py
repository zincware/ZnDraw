# tests/test_settings.py
from zndraw_joblib.settings import JobLibSettings


def test_default_settings():
    settings = JobLibSettings()
    assert settings.allowed_categories == ["modifiers", "selections", "analysis"]
    assert settings.worker_timeout_seconds == 60
    assert settings.sweeper_interval_seconds == 30
    assert settings.long_poll_max_wait_seconds == 60
    assert settings.internal_task_timeout_seconds == 3600


def test_settings_from_env(monkeypatch):
    monkeypatch.setenv("ZNDRAW_JOBLIB_ALLOWED_CATEGORIES", '["custom"]')
    monkeypatch.setenv("ZNDRAW_JOBLIB_WORKER_TIMEOUT_SECONDS", "120")
    settings = JobLibSettings()
    assert settings.allowed_categories == ["custom"]
    assert settings.worker_timeout_seconds == 120


def test_claim_settings_from_env(monkeypatch):
    """Test that claim retry settings can be configured via environment."""
    monkeypatch.setenv("ZNDRAW_JOBLIB_CLAIM_MAX_ATTEMPTS", "20")
    monkeypatch.setenv("ZNDRAW_JOBLIB_CLAIM_BASE_DELAY_SECONDS", "0.02")
    settings = JobLibSettings()
    assert settings.claim_max_attempts == 20
    assert settings.claim_base_delay_seconds == 0.02
