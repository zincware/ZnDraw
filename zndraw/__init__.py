import importlib.metadata

from zndraw.app import app

__all__ = ["app"]

__version__ = importlib.metadata.version("zndraw")
