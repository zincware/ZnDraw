import os

if "ZNDRAW_ENABLE_EVENTLET" in os.environ:
    print("ZnDraw running with eventlet `monkey_patch` enabled.")

    import eventlet

    eventlet.monkey_patch()
    # remove "ZNDRAW_ENABLE_EVENTLET" from the environment
    os.environ.pop("ZNDRAW_ENABLE_EVENTLET")

import importlib.metadata

from zndraw.base import Extension
from zndraw.zndraw import ZnDraw

__all__ = ["ZnDraw", "Extension"]

__version__ = importlib.metadata.version("zndraw")
