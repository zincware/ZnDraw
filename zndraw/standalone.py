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
        import signal

        from zndraw_app.make_celery import celery_app

        # Prevent celery from installing signal handlers by temporarily disabling them
        original_signal = signal.signal

        def dummy_signal(sig, handler):
            if sig in (signal.SIGINT, signal.SIGTERM):
                # Don't let celery install handlers for these signals
                return signal.SIG_DFL
            return original_signal(sig, handler)

        signal.signal = dummy_signal
        try:
            celery_app.worker_main(
                argv=[
                    "worker",
                    "--loglevel=info",
                    "--pool=eventlet",
                ]
            )
        finally:
            signal.signal = original_signal

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
