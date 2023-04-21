import importlib.metadata

from zndraw.app import app
from zndraw.main import ZnDraw

__all__ = ["app", "ZnDraw"]

__version__ = importlib.metadata.version("zndraw")
