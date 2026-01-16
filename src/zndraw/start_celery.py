"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess
import typing as t

if t.TYPE_CHECKING:
    from zndraw.config import ZnDrawConfig


def run_celery_worker(config: "ZnDrawConfig") -> subprocess.Popen:
    """Run a celery worker with proper configuration.

    Serializes the config to environment variables using the pydantic-settings
    nested delimiter format (ZNDRAW_STORAGE__TYPE, etc.) so the worker can
    reconstruct the same config.

    Parameters
    ----------
    config : ZnDrawConfig
        Configuration object containing all settings.

    Returns
    -------
    subprocess.Popen
        Running celery worker process.
    """
    from zndraw.config import LMDBStorageConfig, MongoDBStorageConfig

    my_env = os.environ.copy()

    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    # Core configuration
    if config.redis_url is not None:
        my_env["ZNDRAW_REDIS_URL"] = config.redis_url

    # Storage configuration - use nested delimiter format
    match config.storage:
        case LMDBStorageConfig():
            my_env["ZNDRAW_STORAGE__TYPE"] = "lmdb"
            my_env["ZNDRAW_STORAGE__PATH"] = config.storage.path
            my_env["ZNDRAW_STORAGE__MAP_SIZE"] = str(config.storage.map_size)
        case MongoDBStorageConfig():
            my_env["ZNDRAW_STORAGE__TYPE"] = "mongodb"
            my_env["ZNDRAW_STORAGE__URL"] = config.storage.url
            my_env["ZNDRAW_STORAGE__DATABASE"] = config.storage.database

    # Media path
    my_env["ZNDRAW_MEDIA_PATH"] = config.media_path

    # Server settings
    my_env["ZNDRAW_SERVER_HOST"] = config.server_host
    my_env["ZNDRAW_SERVER_PORT"] = str(config.server_port)
    if config.server_url:
        my_env["ZNDRAW_SERVER_URL"] = config.server_url
    my_env["ZNDRAW_LOG_LEVEL"] = config.log_level

    # Upload & storage
    my_env["ZNDRAW_UPLOAD_TEMP"] = config.upload_temp
    my_env["ZNDRAW_MAX_UPLOAD_MB"] = str(config.max_upload_mb)

    # Optional features
    my_env["ZNDRAW_SIMGEN_ENABLED"] = "true" if config.simgen_enabled else "false"
    my_env["ZNDRAW_FILE_BROWSER_ENABLED"] = (
        "true" if config.file_browser_enabled else "false"
    )
    my_env["ZNDRAW_FILE_BROWSER_ROOT"] = config.file_browser_root

    # Flask secret key (needed for JWT validation in tasks)
    my_env["FLASK_SECRET_KEY"] = config.flask_secret_key

    # Admin credentials
    if config.admin_username is not None:
        my_env["ZNDRAW_ADMIN_USERNAME"] = config.admin_username
    if config.admin_password is not None:
        my_env["ZNDRAW_ADMIN_PASSWORD"] = config.admin_password.get_secret_value()

    # Build celery command
    # -q flag suppresses the startup banner (added in Celery 4.0)
    celery_cmd = [
        "celery",
        "-A",
        "zndraw_cli.celery",
        "worker",
        f"--loglevel={config.log_level.lower()}",
        "--pool=eventlet",
    ]

    # Add -q flag to suppress banner when not in debug mode
    if config.log_level.upper() != "DEBUG":
        celery_cmd.insert(1, "-q")  # -q must come before subcommand

    worker = subprocess.Popen(
        celery_cmd,
        env=my_env,
    )
    return worker


if __name__ == "__main__":
    from zndraw.config import get_config

    config = get_config()
    try:
        worker = run_celery_worker(config)
        worker.wait()
    finally:
        worker.terminate()
        worker.wait()
        worker.kill()
        print("Celery worker terminated.")
