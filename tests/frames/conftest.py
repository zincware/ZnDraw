import random
import socket
import subprocess
import time

import pytest
from znsocket import Client


@pytest.fixture
def eventlet_memory_server():
    port = random.randint(10000, 20000)

    # Start znsocket server subprocess
    server_proc = subprocess.Popen(
        ["znsocket", "--port", str(port)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # Wait for the server to be ready
    for _ in range(100):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
                sock.settimeout(0.2)
                sock.connect(("127.0.0.1", port))
                break
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    else:
        server_proc.terminate()
        server_proc.wait()
        raise TimeoutError("Server did not start in time")

    yield f"znsocket://127.0.0.1:{port}"

    # Clean up
    server_proc.terminate()
    server_proc.wait()


@pytest.fixture
def znsclient(eventlet_memory_server):
    r = Client.from_url(eventlet_memory_server)
    yield r
    r.flushall()
