import random
import shutil
import socket
import subprocess
import time

import eventlet  # noqa - eventlet must be installed for flask-socketio to start a production server
import pytest
import signal


@pytest.fixture
def server():
    port = random.randint(10000, 20000)

    # Start zndraw-server subprocess
    proc = subprocess.Popen(
        ["zndraw-server", "--port", str(port)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Wait for the server to be ready
    for _ in range(100):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                sock.settimeout(0.1)
                sock.connect(("127.0.0.1", port))
                break
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    else:
        proc.kill()
        raise TimeoutError("Server did not start in time")

    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()
            shutil.rmtree("data", ignore_errors=True)
            raise RuntimeError("Server did not shut down in time")
