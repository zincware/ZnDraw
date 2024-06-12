import os

if "ZNDRAW_MONKEY_PATCH" in os.environ:
    import eventlet

    eventlet.monkey_patch()


import importlib.metadata

from zndraw.base import Extension
from zndraw.zndraw import ZnDraw

__all__ = ["ZnDraw", "Extension"]

__version__ = importlib.metadata.version("zndraw")
