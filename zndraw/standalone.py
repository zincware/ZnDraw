"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess
import threading
import time

import znsocket.exceptions


def run_znsocket(port) -> subprocess.Popen:
    """Run a znsocket server instead of redis."""

    server = subprocess.Popen(["znsocket", "--port", str(port)])

    for trial in range(1000):
        try:
            znsocket.Client.from_url(f"znsocket://localhost:{port}")
            break
        except znsocket.exceptions.ConnectionError:
            time.sleep(0.1)
            if trial % 10 == 0:
                print("Waiting for znsocket to start...")
    else:
        raise RuntimeError("Unable to start ZnSocket server!")

    return server


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
