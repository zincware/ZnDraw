import importlib
import json
import uuid

import networkx as nx
import numpy as np
import tqdm
from flask import (Flask, Response, render_template, request, send_file,
                   session, stream_with_context)

from zndraw import io, shared, tools

app = Flask(__name__)
app.secret_key = str(uuid.uuid4())


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    session["key"] = str(uuid.uuid4())  # TODO use session key e.g. for atoms cache
    session["step"] = 0
    shared.bond_method = tools.data.ASEComputeBonds()
    return render_template("index.html", config=shared.config.dict())


@app.route("/config", methods=["POST", "GET"])
def config():
    """Get the zndraw configuration."""
    if request.method == "POST":
        print(f"Updating config: {request.json}")
        for key, value in request.json.items():
            setattr(shared.config, key, value)
    return {
        **shared.config.dict(),
        "total_frames": len(shared.config._atoms_cache) - 1,
    }


@app.route("/select", methods=["POST"])
def select() -> list[int]:
    """Update the selected atoms."""
    step = request.json["step"]
    method = request.json["method"]
    print(f"Selecting atoms {request.json}")
    try:
        selected_ids = [int(x) for x in request.json["selected_ids"]]
    except TypeError:
        selected_ids = []
    if method in ["particles", "none"]:
        return {"selected_ids": selected_ids, "updated": False}

    atoms = shared.config.get_atoms(step)

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
    elif method == "all":
        return {
            "selected_ids": tools.select.select_all(atoms, selected_ids),
            "updated": True,
        }
    else:
        raise ValueError(f"Unknown selection method: {method}")


@app.route("/add_update_function", methods=["POST"])
def add_update_function():
    """Add a function to the config."""
    try:
        signature = shared.config.get_modifier_schema(request.json)
    except Exception as err:
        return {"error": str(err)}
    return signature


@app.route("/update", methods=["POST"])
def update_scene():
    """Update the scene with the selected atoms."""
    # delete via {'selected_ids': [144], 'step': 38, 'modifier': 'zndraw.examples.Delete', 'modifier_kwargs': {}, 'points': []}
    modifier = request.json["modifier"]
    modifier_kwargs = request.json["modifier_kwargs"]
    selected_ids = list(sorted(request.json["selected_ids"]))
    step = request.json["step"]
    points = np.array([[x["x"], x["y"], x["z"]] for x in request.json["points"]])
    shared.config.run_modifier(
        modifier, selected_ids, step, modifier_kwargs, points=points
    )
    return {}


@app.route("/add_analysis", methods=["POST"])
def add_analysis():
    """Add a function to the config."""
    import importlib

    # we need to load the first atoms to get the schema
    atoms = shared.config.get_atoms(0)

    try:
        module_name, function_name = request.json.rsplit(".", 1)
        module = importlib.import_module(module_name)
        cls = getattr(module, function_name)
        try:
            schema = cls.schema_from_atoms(shared.config.atoms_list)
        except AttributeError as e:
            raise ImportError(
                f"Unable to import {cls}. Can not generate schema."
            ) from e
        schema["title"] = request.json
    except (ImportError, ValueError) as err:
        return {"error": str(err)}
    print(f"Adding analysis {schema}")
    return schema


@app.route("/add_bonds", methods=["POST"])
def add_bonds():
    """Add a function to the config."""
    try:
        signature = shared.config.get_modifier_schema(request.json)
    except Exception as err:
        return {"error": str(err)}
    return signature


@app.route("/set_bonds", methods=["POST"])
def set_bonds():
    """Add a function to the config."""
    print(f"Setting bonds {request.json}")
    module_name, function_name = request.json["method"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    shared.bond_method = getattr(module, function_name)(**request.json["bonds_kwargs"])

    if "order" in request.json:
        shared.bond_method.update_bond_order(
            atoms=shared.config.get_atoms(request.json["step"]),
            particles=request.json["selected_ids"],
            order=request.json["order"],
        )
    return {}


@app.route("/analyse", methods=["POST"])
def analyse():
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
    shared.config.load_atoms()
    return {}


@app.route("/download")
def download():
    """Download the current atoms."""

    b = shared.config.export_atoms()
    b.seek(0)
    return send_file(b, download_name="traj.h5", as_attachment=True)


@app.route("/download-selected/<int:step>/<selected_ids>")
def download_selection(step, selected_ids):
    """Download the current atoms."""
    selected_ids = [int(x) for x in selected_ids.split(",")]
    b = shared.config.export_selection(step, selected_ids)
    b.seek(0)
    return send_file(b, download_name="traj.xyz", as_attachment=True)


@app.route("/frame-set", methods=["POST"])
def frame_set():
    session["step"] = request.json["step"]
    print(f"Setting step to {session['step']}")
    return {}


@app.route("/frame-stream")
def frame_stream():
    def generate(step):
        values = list(range(step, step + shared.config.js_frame_buffer[1])) + list(
            range(step, step - shared.config.js_frame_buffer[0], -1)
        )

        stream_id = uuid.uuid4()
        shared.streaming = stream_id

        pbar = tqdm.tqdm(values, desc=f"Streaming {step}", ncols=80, leave=False)
        for idx in pbar:
            try:
                data = {idx: shared.bond_method.get_frame(idx)}
                pbar.set_description(
                    f"Streaming {step} {'+' if idx-step > 0 else '-'} {str(abs(idx-step)).zfill(3)}"
                )
                yield f"data: {json.dumps(data)}\n\n"
                if shared.streaming != stream_id:
                    break
            except KeyError:
                pbar.set_description(
                    f"Streaming {step} {'+' if idx-step > 0 else '-'} ... "
                )

        pbar.close()

        yield f"data: {json.dumps({})} \nretry: 10\n\n"

    return Response(generate(session["step"]), mimetype="text/event-stream")


@app.route("/reset-scene-modifiers")
def reset_scene_modifiers():
    shared.config.reset_scene_modifiers()
    return {}
