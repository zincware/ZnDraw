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


@app.route("/config")
def config():
    return dataclasses.asdict(globals.config)


@app.route("/atoms")
def atoms():
    atoms = globals.config.get_atoms(step=0)
    graph = io.get_graph(atoms)
    return [graph.nodes[idx] | {"id": idx} for idx in graph.nodes]


@app.route("/atoms/<step>")
def atoms_step(step):
    atoms = globals.config.get_atoms(step=int(step))
    graph = io.get_graph(atoms)
    return [graph.nodes[idx] | {"id": idx} for idx in graph.nodes]


@app.route("/atoms/<step>/<atom_id>")
def atom_step(step, atom_id):
    return {}


@app.route("/bonds")
def bonds():
    atoms = globals.config.get_atoms(step=0)
    graph = io.get_graph(atoms)
    return list(graph.edges)


@app.route("/bonds/<step>/")
def bonds_step(step):
    atoms = globals.config.get_atoms(step=int(step))
    graph = io.get_graph(atoms)
    return list(graph.edges)


@app.route("/bonds/<step>/<bond_id>")
def bond_step(step, bond_id):
    return {}


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
    # content = request.json
    try:
        if atom_id in session["selected"]:
            session["selected"] = [x for x in session["selected"] if x != atom_id]
        else:
            session["selected"] = session["selected"] + [atom_id]
    except KeyError:
        session["selected"] = [atom_id]

    session["updated"] = True

    print(session["selected"])
    return {}


@app.route("/update", methods=["GET", "POST"])
def update_scene():
    # content = request.json

    if "selected" not in session:
        return []

    if not session["updated"]:
        return []

    function = globals.config.get_update_function()
    atoms = function(
        [int(x) for x in session["selected"]], globals.atoms
    )  # TODO animation

    session["updated"] = False

    return np.array([x.get_positions() for x in atoms]).tolist()

    # graph = io.get_graph(atoms)
    # globals.atoms = atoms
    # globals.graph = graph

    # data = {"nodes": [], "edges": []}
    # for node in graph.nodes:
    #     data["nodes"].append(graph.nodes[node] | {"id": node})

    # for edge in graph.edges:
    #     data["edges"].append(edge)

    # return data
