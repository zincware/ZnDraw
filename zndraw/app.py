from flask import Flask
from flask import render_template
from flask import session
import uuid
import numpy as np
from zndraw import globals, io
import networkx as nx

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/xyz")
def xyz():
    graph = io.read_file(globals.file)
    globals.graph = graph

    session["graph"] = nx.node_link_data(graph)
    data = {"nodes": [], "edges": []}
    for node in graph.nodes:
        data["nodes"].append(graph.nodes[node])

    for edge in graph.edges:
        pos_1 = np.array(graph.nodes[edge[0]]["position"])

        pos_2 = np.array(graph.nodes[edge[1]]["position"])

        ref = pos_1 - pos_2

        data["edges"].append(
            {
                "radius": 0.15,
                "height": np.linalg.norm(pos_1 - pos_2).item(),
                "position": ((pos_1 + pos_2) / 2).tolist(),
                "axis": (ref / np.linalg.norm(ref)).tolist(),
            }
        )
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
