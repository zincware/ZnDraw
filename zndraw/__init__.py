import os

import eventlet

eventlet.monkey_patch()

os.environ["FLASK_ENV"] = "development"
os.environ["FLASK_PORT"] = "3141"
# os.environ["FLASK_DEBUG"] = "1"
os.environ["FLASK_STORAGE"] = "redis://localhost:6379/0"
os.environ["FLASK_AUTH_TOKEN"] = "ABC"
os.environ["FLASK_TUTORIAL"] = "https://google.com"
os.environ["FLASK_SIMGEN"] = "TRUE"

from zndraw.zndraw import ZnDraw

__all__ = ["ZnDraw"]
