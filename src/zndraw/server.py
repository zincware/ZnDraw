import logging

import redis
from flask import Flask
from flask_socketio import SocketIO

log = logging.getLogger(__name__)

socketio = SocketIO()


def create_app():
    app = Flask(__name__)

    from .app import main as main_blueprint

    app.register_blueprint(main_blueprint)

    socketio.init_app(app, cors_allowed_origins="*")
    r = redis.Redis(decode_responses=True)

    app.config["SECRET_KEY"] = "your_secret_key"
    app.config["redis"] = r

    return app, socketio
