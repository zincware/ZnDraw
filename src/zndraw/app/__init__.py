from flask import Blueprint

main = Blueprint("main", __name__)

from . import events, routes  # noqa: E402

__all__ = ["main", "routes", "events"]
