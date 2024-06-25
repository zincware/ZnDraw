"""Utils for running ZnDraw standalone, without redis or external celery worker."""

import os
import platform
import subprocess
import time

import znsocket.exceptions


def run_znsocket(port) -> subprocess.Popen:
    """Run a znsocket server instead of redis."""

    server = subprocess.Popen(["znsocket", str(port)])

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
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            "-P",
            "eventlet",
        ],
        env=my_env,
    )
    return worker
