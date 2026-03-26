# src/zndraw_joblib/settings.py
from pydantic_settings import BaseSettings, SettingsConfigDict


class JobLibSettings(BaseSettings):
    model_config = SettingsConfigDict(env_prefix="ZNDRAW_JOBLIB_")

    allowed_categories: list[str] = ["modifiers", "selections", "analysis"]
    worker_timeout_seconds: int = 60
    sweeper_interval_seconds: int = 30
    long_poll_max_wait_seconds: int = 60

    # Task claim retry settings (for handling concurrent claim contention)
    claim_max_attempts: int = 10
    claim_base_delay_seconds: float = 0.01  # 10ms

    # Internal taskiq worker settings
    internal_task_timeout_seconds: int = 3600  # 1 hour

    # Provider settings
    allowed_provider_categories: list[str] | None = None  # None = unrestricted
    provider_result_ttl_seconds: int = 300
    provider_inflight_ttl_seconds: int = 30
    provider_long_poll_default_seconds: int = 5
    provider_long_poll_max_seconds: int = 30
