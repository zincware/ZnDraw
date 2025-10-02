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
            "src.zndraw.app.make_celery",
            "worker",
            "--loglevel=info",
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
