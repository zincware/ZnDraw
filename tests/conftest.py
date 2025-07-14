import os
import platform
import random
import signal
import socket
import subprocess
import time

import ase.build
import ase.collections
import numpy as np
import pytest
from ase.calculators.singlepoint import SinglePointCalculator


@pytest.fixture
def server():
    port = random.randint(10000, 20000)

    # Set up environment variables
    my_env = os.environ.copy()
    my_env["FLASK_PORT"] = str(port)
    my_env["FLASK_STORAGE"] = f"znsocket://localhost:{port}"
    my_env["FLASK_SERVER_URL"] = f"http://localhost:{port}"

    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    def kill_child_processes(parent_pid, context=""):
        """Kill all child processes of the given parent PID"""
        try:
            # Get child processes
            result = subprocess.run(
                ["pgrep", "-P", str(parent_pid)],
                capture_output=True,
                text=True,
                timeout=2,
            )
            if result.returncode == 0:
                child_pids = result.stdout.strip().split("\n")
                for child_pid in child_pids:
                    if child_pid:
                        try:
                            # Get process name for logging
                            try:
                                name_result = subprocess.run(
                                    ["ps", "-p", child_pid, "-o", "comm="],
                                    capture_output=True,
                                    text=True,
                                    timeout=1,
                                )
                                process_name = (
                                    name_result.stdout.strip()
                                    if name_result.returncode == 0
                                    else "unknown"
                                )
                            except:
                                process_name = "unknown"

                            print(
                                f"Killing child process {child_pid} ({process_name}){context}"
                            )
                            os.kill(int(child_pid), signal.SIGTERM)
                        except (ProcessLookupError, ValueError):
                            pass
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            pass

    # Start zndraw with standalone=True (this will start celery worker automatically)
    server_proc = subprocess.Popen(
        ["zndraw", "--port", str(port), "--no-browser"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=my_env,
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
        # Clean up if server failed to start
        kill_child_processes(server_proc.pid, " - startup failure cleanup")

        # Then kill the main server process
        server_proc.terminate()
        try:
            server_proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            server_proc.kill()
        raise TimeoutError("Server did not start in time")

    yield f"http://127.0.0.1:{port}"

    # use time.sleep and disbale the no-browser flag to see the results (might interfere with tests)
    # time.sleep(10)

    # Clean up: kill the server process and any child processes
    # Kill child processes first (celery workers, znsocket, etc.)
    kill_child_processes(server_proc.pid)

    # Then kill the main server process
    server_proc.terminate()
    try:
        server_proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        server_proc.kill()

    # Double check - kill any remaining child processes after main process cleanup
    time.sleep(0.5)
    kill_child_processes(server_proc.pid, " - post-cleanup verification")


@pytest.fixture
def s22() -> list[ase.Atoms]:
    """S22 dataset."""
    return list(ase.collections.s22)


@pytest.fixture
def water() -> ase.Atoms:
    """Water molecule."""
    return ase.build.molecule("H2O")


@pytest.fixture
def s22_energy_forces() -> list[ase.Atoms]:
    images = []
    for atoms in ase.collections.s22:
        calc = SinglePointCalculator(
            atoms, energy=np.random.rand(), forces=np.random.rand(len(atoms), 3)
        )
        atoms.calc = calc
        images.append(atoms)
    return images
