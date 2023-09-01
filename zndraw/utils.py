import socket

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