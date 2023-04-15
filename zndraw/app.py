from flask import Flask
from flask import render_template
from flask import session, request
import uuid
import numpy as np
from zndraw import globals, io
import dataclasses

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/xyz")
def xyz():
    graph = io.read_file(globals.config.file)
    globals.graph = graph

    data = {"nodes": [], "edges": []}
    for node in graph.nodes:
        data["nodes"].append(graph.nodes[node] | {"id": node})

    for edge in graph.edges:
        data["edges"].append(edge)

    return data


@app.route("/config")
def config():
    return dataclasses.asdict(globals.config)


def get_atoms():
    import ase.io

    yield from ase.io.iread(globals.config.file)


atoms_iter = iter(get_atoms())


@app.route("/animation")
def get_position_updates():
    try:
        return [next(atoms_iter).get_positions().tolist() for _ in range(10)]
    except StopIteration:
        return {}


@app.route("/atom/<atom_id>", methods=["GET", "POST"])
def add_message(atom_id):
    content = request.json
    print(content)
    session["selected"] = atom_id
    session["step"] = 0
    return {}


@app.route("/update")
def new_atoms():
    if "selected" not in session:
        return {}
    if session["step"] > 10:
        return {}
    node = globals.graph.nodes[int(session["selected"])]
    positions = np.array(node["position"]) + np.ones(3) * session["step"]
    session["step"] += 1
    return {
        "id": session["selected"],
        "position": positions.tolist(),
        "radius": 0.1,
        "color": "limegreen",
    }
