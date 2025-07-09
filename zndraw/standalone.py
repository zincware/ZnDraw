"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess
import threading


def run_celery_thread_worker() -> threading.Thread:
    """Run a celery worker."""
    my_env = os.environ.copy()
    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    def run_celery_worker():
        from zndraw_app.make_celery import celery_app

        celery_app.worker_main(
            argv=["worker", "--loglevel=info", "--without-gossip", "--pool=eventlet"]
        )

    worker = threading.Thread(target=run_celery_worker)
    worker.start()
    return worker


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
            "-P",
            "eventlet",
        ],
        env=my_env,
    )
    return worker
