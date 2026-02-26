"""ZnDraw - Interactive visualization for atomistic simulations."""

from zndraw.client import ZnDraw
from zndraw.tqdm import ZnDrawTqdm

try:
    from zndraw._version import __version__, __version_tuple__
except ImportError:
    # Fallback for development without build
    __version__ = "0.0.0-dev"
    __version_tuple__ = (0, 0, 0, "dev")

__all__ = ["ZnDraw", "ZnDrawTqdm", "__version__", "__version_tuple__"]
