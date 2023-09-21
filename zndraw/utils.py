import logging
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


class ZnDrawLoggingHandler(logging.Handler):
    """Logging handler which emits log messages to the ZnDraw server."""

    def __init__(self, socket):
        super().__init__()
        self.socket = socket

    def emit(self, record):
        try:
            msg = self.format(record)
            self.socket.emit("message:log", msg)
        except RecursionError:  # See StreamHandler
            raise
        except Exception:
            self.handleError(record)
