import uuid

import networkx as nx
from flask import Flask, render_template, request, session

from zndraw import globals, io

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    session["key"] = str(uuid.uuid4())  # TODO use session key e.g. for atoms cache
    return render_template("index.html", config=globals.config.dict())


@app.route("/config")
def config():
    """Get the zndraw configuration."""
    return {
        **globals.config.dict(),
        "total_frames": len(globals._atoms_cache) - 1,
    }

@app.route("/graph", methods=["POST"])
def get_graph():
    step = request.json
    try:
        atoms = globals.config.get_atoms(step=int(step))
        graph = io.get_graph(atoms)
        return {
            "nodes": [{**graph.nodes[idx], "id": idx} for idx in graph.nodes],
            "edges": list(graph.edges),
        }
    except KeyError:
        return {}


@app.route("/positions", methods=["POST"])
def positions_step():
    params = request.json
    result = []
    try:
        for step in range(params["start"], params["stop"]):
            atoms = globals.config.get_atoms(step=int(step))
            result.append(atoms.get_positions().tolist())
        return result
    except KeyError:
        return result

@app.route("/select", methods=["POST"])
def select() -> list[int]:
    """Update the selected atoms."""
    step = request.json["step"]
    method = request.json["method"]
    selected_ids = request.json["selected_ids"]
    if method == "particles":
        return {"selected_ids": selected_ids, "updated": False}
    elif method == "species":
        atoms = globals.config.get_atoms(step)

        for id in tuple(selected_ids):
            selected_symbol = atoms[id].symbol
            selected_ids += [
                idx for idx, atom in enumerate(atoms) if atom.symbol == selected_symbol
            ]
        return {"selected_ids": list(set(selected_ids)), "updated": True}
    elif method == "connected":
        atoms = globals.config.get_atoms(step)
        graph = io.get_graph(atoms)

        total_ids = []

        for node_id in selected_ids:
            total_ids += list(nx.node_connected_component(graph, node_id))

        return {"selected_ids": list(set(total_ids)), "updated": True}

    else:
        raise ValueError(f"Unknown selection method: {method}")


@app.route("/add_update_function", methods=["POST"])
def add_update_function():
    """Add a function to the config."""
    globals.config.update_function = request.json
    try:
        signature = globals.config.get_update_signature()
    except (ImportError, ValueError) as err:
        return {"error": str(err)}
    return signature


@app.route("/update_function_values", methods=["POST"])
def update_function_values():
    """Update the values of the update function."""
    globals.config.set_update_function_parameters(request.json)
    return {}


@app.route("/select_update_function/<name>")
def select_update_function(name):
    """Select a function from the config."""
    if name == "none":
        name = None
    globals.config.update_function_name = name
    return {}


@app.route("/update", methods=["POST"])
def update_scene():
    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]

    function = globals.config.get_update_function()
    atoms = function(
        [int(x) for x in selected_ids], globals.config.get_atoms(step=int(step))
    )

    offset = len(globals._atoms_cache)

    for idx, atom in enumerate(atoms):
        globals._atoms_cache[idx + offset] = atom

    # this has to return before the scene is automatically updated
    return {}


@app.route("/load")
def load():
    """Function to call asynchronously to load atoms in the background."""
    globals.config.load_atoms()
    return {}
