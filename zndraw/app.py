from flask import Flask
from flask import render_template
from flask import request
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
    data = {"nodes": [], "edges": []}
    for node in graph.nodes:
        data["nodes"].append(graph.nodes[node])

    for edge in graph.edges:
        pos_1 = np.array(
            [
                graph.nodes[edge[0]]["x"],
                graph.nodes[edge[0]]["y"],
                graph.nodes[edge[0]]["z"],
            ]
        )

        pos_2 = np.array(
            [
                graph.nodes[edge[1]]["x"],
                graph.nodes[edge[1]]["y"],
                graph.nodes[edge[1]]["z"],
            ]
        )
        # TODO why split x/y/z and not have positions

        ref = pos_1 - pos_2

        data["edges"].append(
            {
                "radius": 0.0001,
                "height": np.linalg.norm(pos_1 - pos_2).item(),
                "position": ((pos_1 + pos_2) / 2).tolist(),
                "axis": (ref / np.linalg.norm(ref)).tolist(),
            }
        )
    return data


@app.route("/atom/<uuid>", methods=["GET", "POST"])
def add_message(uuid):
    content = request.json
    _ = content
    return {}
