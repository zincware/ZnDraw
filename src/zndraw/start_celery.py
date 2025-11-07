"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess


# We use this for running tests for now
def run_celery_worker(config: "ZnDrawConfig | None" = None, redis_url: str | None = None) -> subprocess.Popen:
    """Run a celery worker with proper configuration.

    Parameters
    ----------
    config : ZnDrawConfig | None
        Configuration object containing all settings. If provided, all config
        values are passed as environment variables to the worker.
    redis_url : str | None
        Backwards compatibility parameter. If config is not provided, only redis_url
        is passed to the worker.

    Returns
    -------
    subprocess.Popen
        Running celery worker process.
    """
    my_env = os.environ.copy()

    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    # Handle backwards compatibility: if config is a string, treat it as redis_url
    if isinstance(config, str):
        redis_url = config
        config = None

    # If config object provided, dump all config to environment
    if config is not None:
        # Core configuration
        if config.redis_url is not None:
            my_env["ZNDRAW_REDIS_URL"] = config.redis_url
        my_env["ZNDRAW_STORAGE_PATH"] = config.storage_path
        my_env["ZNDRAW_SERVER_HOST"] = config.server_host
        my_env["ZNDRAW_SERVER_PORT"] = str(config.server_port)
        my_env["ZNDRAW_SERVER_URL"] = config.server_url
        my_env["ZNDRAW_LOG_LEVEL"] = config.log_level

        # Upload & storage
        my_env["ZNDRAW_UPLOAD_TEMP"] = config.upload_temp
        my_env["ZNDRAW_MAX_UPLOAD_MB"] = str(config.max_upload_mb)

        # Optional features
        my_env["ZNDRAW_SIMGEN_ENABLED"] = "true" if config.simgen_enabled else "false"
        my_env["ZNDRAW_FILE_BROWSER_ENABLED"] = "true" if config.file_browser_enabled else "false"
        my_env["ZNDRAW_FILE_BROWSER_ROOT"] = config.file_browser_root
        my_env["ZNDRAW_CELERY_ENABLED"] = "true" if config.celery_enabled else "false"

        # Flask secret key (needed for JWT validation in tasks)
        my_env["FLASK_SECRET_KEY"] = config.flask_secret_key
    elif redis_url is not None:
        # Backwards compatibility: only set redis_url if provided
        my_env["ZNDRAW_REDIS_URL"] = redis_url

    worker = subprocess.Popen(
        # eventlet worker - use zndraw_cli.celery for proper monkey patching
        [
            "celery",
            "-A",
            "zndraw_cli.celery",
            "worker",
            "--loglevel=info",
            "--pool=eventlet",

        ],
        env=my_env,
    )
    return worker


if __name__ == "__main__":
    try:
        worker = run_celery_worker()
        worker.wait()
    finally:
        worker.terminate()
        worker.wait()
        worker.kill()
        print("Celery worker terminated.")
