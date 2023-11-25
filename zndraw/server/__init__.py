# example https://github.com/miguelgrinberg/Flask-SocketIO-Chat

from flask import Blueprint

main = Blueprint("main", __name__)

from . import routes, events  # noqa: E402, F401, I001
