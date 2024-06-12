# Only for testing, these should be set in the environment.
# import os
# os.environ["FLASK_ENV"] = "development"
# os.environ["FLASK_PORT"] = "3141"
# # os.environ["FLASK_DEBUG"] = "1"
# os.environ["FLASK_STORAGE"] = "redis://localhost:6379/0"
# os.environ["FLASK_AUTH_TOKEN"] = "ABC"
# os.environ["FLASK_TUTORIAL"] = "https://google.com"
# os.environ["FLASK_SIMGEN"] = "TRUE"

from zndraw.base import Extension
from zndraw.zndraw import ZnDraw

__all__ = ["ZnDraw", "Extension"]
