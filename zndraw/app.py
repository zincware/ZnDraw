import uuid

from flask import Flask
from flask_socketio import SocketIO
from flask_caching import Cache

socketio = SocketIO()
cache = Cache(config={'CACHE_TYPE': 'SimpleCache', 'CACHE_DEFAULT_TIMEOUT': 60 * 60 * 24})


def create_app(
    use_token, upgrade_insecure_requests, compute_bonds, tutorial: str, auth_token: str
) -> Flask:
    """Create the Flask app."""

    app = Flask(__name__)
    app.config["SECRET_KEY"] = str(uuid.uuid4())
    app.config["ROOM_HOSTS"] = {}
    app.config["DEFAULT_PYCLIENT"] = None
    app.config["TUTORIAL"] = tutorial
    # dict of {uuid: sid} for each client
    app.config["pyclients"] = {}
    # dict of {token: dict}
    app.config["PER-TOKEN-DATA"] = {}
    app.config["AUTH_TOKEN"] = auth_token
    app.config["MODIFIER"] = {"default_schema": {}, "active": None, "queue": []}

    if not use_token:  # TODO: handle this differently
        app.config["token"] = "notoken"
    app.config["upgrade_insecure_requests"] = upgrade_insecure_requests
    app.config["compute_bonds"] = compute_bonds

    from .server import main as main_blueprint

    app.register_blueprint(main_blueprint)

    cache.init_app(app)
    socketio.init_app(app, cors_allowed_origins="*")
    return app
