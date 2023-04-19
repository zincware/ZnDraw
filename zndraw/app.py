import dataclasses
import uuid

from flask import Flask, make_response, render_template, request, session

from zndraw import globals, io

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    session["key"] = str(uuid.uuid4())  # TODO use session key e.g. for atoms cache
    return render_template("index.html", config=dataclasses.asdict(globals.config))


@app.route("/config")
def config():
    return {
        **dataclasses.asdict(globals.config),
        "total_frames": len(globals._atoms_cache) - 1,
    }


@app.route("/atoms", methods=["POST"])
def atoms_step():
    step = request.json
    try:
        atoms = globals.config.get_atoms(step=int(step))
        graph = io.get_graph(atoms)
        return [{**graph.nodes[idx], "id": idx} for idx in graph.nodes]
    except (KeyError, IndexError):
        return []


@app.route("/positions", methods=["POST"])
def positions_step():
    params = request.json
    result = []
    try:
        for step in range(params["start"], params["stop"]):
            atoms = globals.config.get_atoms(step=int(step))
            result.append(atoms.get_positions().tolist())
        return result
    except (KeyError, IndexError):
        return result


@app.route("/bonds", methods=["POST"])
def bonds_step():
    step = request.json
    atoms = globals.config.get_atoms(step=int(step))
    graph = io.get_graph(atoms)
    return list(graph.edges)


@app.route("/select", methods=["POST"])
def select() -> list[int]:
    """Update the selected atoms."""
    return (
        request.json
    )  # + [x + 1 for x in request.json] + [x - 1 for x in request.json]


@app.route("/update", methods=["POST"])
def update_scene():
    selected_ids = request.json["selected_ids"]
    step = request.json["step"]

    function = globals.config.get_update_function()
    atoms = function(
        [int(x) for x in selected_ids], globals.config.get_atoms(step=int(step))
    )

    offset = len(globals._atoms_cache)

    for idx, atom in enumerate(atoms):
        print(f"processing {idx}")
        globals._atoms_cache[idx + offset] = atom

    # this has to return before the scene is automatically updated
    return {}


@app.route("/load")
def load():
    """Function to call asynchronously to load atoms in the background."""
    globals.config.load_atoms()
    return {}
