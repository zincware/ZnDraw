import sys


def run_healthcheck(server: str):
    """Check the health of a server."""
    from zndraw import ZnDraw

    _ = ZnDraw(url=server, token="healthcheck")
    sys.exit(0)
