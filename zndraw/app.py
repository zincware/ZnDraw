import uuid

import tqdm
import znh5md
from flask import Flask, redirect, render_template, session, url_for
from flask_socketio import SocketIO, emit

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

socketio = SocketIO(app)


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/read/<filename>")
def read(filename):
    """Read a file."""
    filename = filename.replace(".slash.", "/").replace("\\", ".slash.")
    session["filename"] = filename
    # TODO convert file to H5MD, if not already

    # dataset = znh5md.ASEH5MD(filename)
    # for idx in tqdm.tqdm(range(1000)):
    #     _ = dataset.get_atoms_list(idx)

    return redirect(url_for("index"))


@app.route("/session")
def session_view():
    """View the session."""
    return dict(session)


@socketio.on("configuration:id")
def handle_get_id_on_configurations(json):
    print("received json: " + str(json))

    ds = znh5md.ASEH5MD(session["filename"])
    atoms = ds.get_atoms_list(int(json["id"]))

    emit(
        "configuration:id",
        {
            "positions": atoms.get_positions().tolist(),
            "cell": atoms.get_cell().tolist(),
            "atomic_numbers": atoms.get_atomic_numbers().tolist(),
        },
    )
