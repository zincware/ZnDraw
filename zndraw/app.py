from flask import Flask
from flask import render_template
from flask import session
import uuid
import numpy as np
from zndraw import globals, io

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/xyz")
def xyz():
    graph = io.read_file(globals.file)
    globals.graph = graph

    data = {"nodes": [], "edges": []}
    for node in graph.nodes:
        data["nodes"].append(graph.nodes[node] | {"id": node})
    
    for edge in graph.edges:
        data["edges"].append(edge)

    return data


@app.route("/atom/<atom_id>", methods=["GET", "POST"])
def add_message(atom_id):
    # content = request.json
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
