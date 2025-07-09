"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess


# We use this for running tests for now
def run_celery_worker() -> subprocess.Popen:
    """Run a celery worker."""
    my_env = os.environ.copy()
    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw_app.make_celery",
            "worker",
            "--loglevel=info",
            "--pool=eventlet",
        ],
        env=my_env,
    )
    return worker
