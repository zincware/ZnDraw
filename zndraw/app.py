import base64
import dataclasses
import uuid
from io import BytesIO

import matplotlib.pyplot as plt
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


@app.route("/atoms/<step>")
def atoms_step(step):
    try:
        atoms = globals.config.get_atoms(step=int(step))
        graph = io.get_graph(atoms)
        return [{**graph.nodes[idx], "id": idx} for idx in graph.nodes]
    except (KeyError, IndexError):
        return []


@app.route("/positions/<start>&<stop>")
def positions_step(start, stop):
    result = []
    try:
        for step in range(int(start), int(stop)):
            atoms = globals.config.get_atoms(step=int(step))
            result.append(atoms.get_positions().tolist())
        return result
    except (KeyError, IndexError):
        return result


@app.route("/bonds/<step>/")
def bonds_step(step):
    atoms = globals.config.get_atoms(step=int(step))
    graph = io.get_graph(atoms)
    return list(graph.edges)


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


@app.route("/distance/<id1>+<id2>&<step>")
def distances(id1, id2, step):
    # Generate the figure **without using pyplot**.

    # TODO this could be a list of python functions mapped to some keys that take
    # selected indices and return some properties e.g. as plots
    fig, ax = plt.subplots(figsize=(4.5, 2.5), dpi=100)

    distances = []
    for atom in globals.config.get_atoms_list():
        try:
            distances.append(atom.get_distance(int(id1), int(id2)))
        except IndexError:
            distances.append(0)
    ax.plot(distances, label=f"avg {sum(distances)/len(distances):.2f}")
    ax.axvline(
        x=int(step), color="black", label=f"step {step} ({distances[int(step)]:.2f}"
    )
    ax.set_xlabel("Frame")
    ax.set_ylabel("Distance")
    ax.set_title(f"Distance between {id1} and {id2}")
    ax.legend()

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    # Embed the result in the html output.
    data = base64.b64encode(buf.getbuffer()).decode("ascii")
    return f"<img src='data:image/png;base64,{data}'/>"
