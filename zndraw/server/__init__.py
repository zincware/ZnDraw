
# example https://github.com/miguelgrinberg/Flask-SocketIO-Chat

from flask import Blueprint

main = Blueprint('main', __name__)

from . import routes, events
