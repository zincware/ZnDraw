from flask import Flask
from flask import render_template
from flask import session, request
import uuid
from zndraw import globals, io
import dataclasses

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    session["key"] = str(uuid.uuid4())  # TODO use session key e.g. for atoms cache
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
    try:
        atoms = globals.config.get_atoms(step=int(step))
        graph = io.get_graph(atoms)
        return [graph.nodes[idx] | {"id": idx} for idx in graph.nodes]
    except (KeyError, IndexError):
        return []


@app.route("/atoms/<start>&<stop>")
def atoms_steps(start, stop):
    try:
        result = []
        for step in range(int(start), int(stop)):
            atoms = globals.config.get_atoms(step=int(step))
            graph = io.get_graph(atoms)
            result.append([graph.nodes[idx] | {"id": idx} for idx in graph.nodes])
        return result
    except (KeyError, IndexError):
        return []


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


@app.route("/select", methods=["POST"])
def select():
    session["selected"] = request.json
    return {}


@app.route("/update/<step>")
def update_scene(step):
    if "selected" not in session:
        return []

    function = globals.config.get_update_function()
    atoms = function(
        [int(x) for x in session["selected"]], globals.config.get_atoms(step=int(step))
    )  # TODO animation

    globals._atoms_cache |= {
        idx + len(globals._atoms_cache): atom for idx, atom in enumerate(atoms)
    }
    return {}
