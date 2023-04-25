import uuid

import networkx as nx
import numpy as np
from flask import Flask, render_template, request, session

from zndraw import globals, io

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    session["key"] = str(uuid.uuid4())  # TODO use session key e.g. for atoms cache
    return render_template("index.html", config=globals.config.dict())


@app.route("/config", methods=["POST", "GET"])
def config():
    """Get the zndraw configuration."""
    if request.method == "POST":
        print(f"Updating config: {request.json}")
        for key, value in request.json.items():
            setattr(globals.config, key, value)
    return {
        **globals.config.dict(),
        "total_frames": len(globals.config._atoms_cache) - 1,
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
            "box": atoms.get_cell().diagonal().tolist(),
        }
    except KeyError:
        return {}


@app.route("/data", methods=["POST"])
def positions_step():
    params = request.json
    result = {"position": [], "force": [], "box": []}
    try:
        for step in range(params["start"], params["stop"]):
            atoms = globals.config.get_atoms(step=int(step))
            result["position"].append(atoms.get_positions().tolist())
            result["box"].append(atoms.get_cell().diagonal().tolist())
            # TODO MAKE THIS OPTIONAL!!, also energy, etc.
            # try:
            #     result["force"].append(atoms.get_forces().tolist())
            # except:
            #     result["force"].append(np.zeros_like(atoms.get_positions()).tolist())
        return result
    except KeyError:
        return result


@app.route("/select", methods=["POST"])
def select() -> list[int]:
    """Update the selected atoms."""
    step = request.json["step"]
    method = request.json["method"]
    try:
        selected_ids = [int(x) for x in request.json["selected_ids"]]
    except TypeError:
        selected_ids = []
    if method in ["particles", "none"]:
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
    try:
        signature = globals.config.add_update_function(request.json)
    except (ImportError, ValueError) as err:
        return {"error": str(err)}
    return signature


@app.route("/set_update_function_parameter", methods=["POST"])
def set_update_function_parameter():
    """Update the values of the update function."""
    globals.config.set_update_function_parameter(request.json)
    return {}


@app.route("/select_update_function/<name>")
def select_update_function(name):
    """Select a function from the config."""
    if name == "none":
        name = None
    globals.config.active_update_function = name
    return {}


@app.route("/update", methods=["POST"])
def update_scene():
    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]
    points = np.array([[x["x"], x["y"], x["z"]] for x in request.json["points"]])
    globals.config.apply_update_function(selected_ids, step, points=points)
    return {}


@app.route("/analyse", methods=["POST"])
def analyse():
    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]
    from zndraw.analyse import get_distance_plot

    fig = get_distance_plot(step, selected_ids)
    return fig.to_json()


@app.route("/load")
def load():
    """Function to call asynchronously to load atoms in the background."""
    globals.config.load_atoms()
    return {}
