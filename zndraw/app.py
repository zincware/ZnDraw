from flask import Flask
from flask_socketio import SocketIO
import uuid

socketio = SocketIO()

def create_app(token, upgrade_insecure_requests, compute_bonds) -> Flask:
    """Create the Flask app."""

    app = Flask(__name__)
    app.config["SECRET_KEY"] = str(uuid.uuid4())
    app.config["ROOM_HOSTS"] = {}
    app.config["DEFAULT_PYCLIENT"] = None
    app.config["MODIFIER"] = {"default_schema": {}}

    if token is not None:
        app.config["token"] = token
    app.config["upgrade_insecure_requests"] = upgrade_insecure_requests
    app.config["compute_bonds"] = compute_bonds


    from .server import main as main_blueprint
    app.register_blueprint(main_blueprint)

    socketio.init_app(app, cors_allowed_origins="*")
    return app
