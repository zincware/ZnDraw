
# import eventlet

# eventlet.monkey_patch()

import uuid

from flask import Flask, render_template
from flask_socketio import SocketIO, emit

from zndraw.data import get_atomsdict_list

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

socketio = SocketIO(app)



@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/exit")
def exit():
    """Exit the session."""
    socketio.stop()
    return "Server shutting down..."

@socketio.on("atoms:request")
def atoms_request(data):
    """Return the atoms."""
    print(f"atoms:request {data = }")
    for atoms_dict in get_atomsdict_list(app.config["filename"]):
        emit("atoms:upload", atoms_dict)

# @app.route("/download/<int:idx>")
# def download(idx):
#     socketio.emit("atoms:download", [idx])
#     return "OK"

# @socketio.on("atoms:download")
# def atoms_download(data):
#     """Return the atoms."""
#     print(f"atoms:download {data = }")
