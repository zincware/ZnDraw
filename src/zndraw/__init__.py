"""ZnDraw communication package with physical-to-logical frame mapping."""

import logging

from zndraw.client import Client
from zndraw.server import create_app

__all__ = ["Client", "create_app"]

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
log.addHandler(handler)
