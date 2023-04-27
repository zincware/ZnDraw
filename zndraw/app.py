import uuid

import networkx as nx
import numpy as np
from flask import Flask, render_template, request, session

from zndraw import globals, io, tools

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
    return tools.data.serialize_atoms(request.json["start"], request.json["stop"])


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

    atoms = globals.config.get_atoms(step)

    if method == "species":
        return {
            "selected_ids": tools.select.select_identical_species(atoms, selected_ids),
            "updated": True,
        }
    elif method == "connected":
        return {
            "selected_ids": tools.select.select_connected(atoms, selected_ids),
            "updated": True,
        }
    else:
        raise ValueError(f"Unknown selection method: {method}")


@app.route("/add_update_function", methods=["POST"])
def add_update_function():
    """Add a function to the config."""
    try:
        signature = globals.config.get_modifier_schema(request.json)
    except (ImportError, ValueError) as err:
        return {"error": str(err)}
    return signature


@app.route("/update", methods=["POST"])
def update_scene():
    """Update the scene with the selected atoms."""
    modifier = request.json["modifier"]
    modifier_kwargs = request.json["modifier_kwargs"]
    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]
    points = np.array([[x["x"], x["y"], x["z"]] for x in request.json["points"]])
    globals.config.run_modifier(
        modifier, selected_ids, step, modifier_kwargs, points=points
    )
    return {}


@app.route("/add_analysis", methods=["POST"])
def add_analysis():
    """Add a function to the config."""
    import importlib

    try:
        module_name, function_name = request.json.rsplit(".", 1)
        module = importlib.import_module(module_name)
        cls = getattr(module, function_name)()
        schema = cls.schema_from_atoms(globals.config._atoms_cache)
        schema["title"] = request.json
    except (ImportError, ValueError) as err:
        return {"error": str(err)}
    return schema


@app.route("/analyse", methods=["POST"])
def analyse():
    import importlib

    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]
    modifier = request.json["modifier"]
    modifier_kwargs = request.json["modifier_kwargs"]
    module_name, function_name = modifier.rsplit(".", 1)
    module = importlib.import_module(module_name)
    cls = getattr(module, function_name)
    instance = cls(**modifier_kwargs)
    fig = instance.run(selected_ids)
    return fig.to_json()


@app.route("/load")
def load():
    """Function to call asynchronously to load atoms in the background."""
    globals.config.load_atoms()
    return {}
