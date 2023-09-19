import contextlib
import socket

import ase


def get_port(default: int = 1234) -> int:
    """Get an open port."""
    try:
        sock = socket.socket()
        sock.bind(("", default))
        port = 1234
    except OSError:
        sock = socket.socket()
        sock.bind(("", 0))
        port = sock.getsockname()[1]
    finally:
        sock.close()
    return port


@contextlib.contextmanager
def set_global_atoms(atoms: ase.Atoms):
    """Temporarily create a global 'ATOMS' variable."""
    global ATOMS
    ATOMS = atoms
    yield
    del ATOMS
