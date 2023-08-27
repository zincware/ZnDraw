import importlib
import uuid
from io import StringIO

import numpy as np
import tqdm
from flask import Flask, render_template
from flask_socketio import SocketIO, emit

from zndraw.data import atoms_from_json, atoms_to_json, get_atomsdict_list
from zndraw.settings import GlobalConfig

app = Flask(__name__)
app.config["SECRET_KEY"] = str(uuid.uuid4())

io = SocketIO(app, max_http_buffer_size=int(1e10))
# 10 GB Upload limit


@app.route("/")
def index():
    """Render the main ZnDraw page."""
    return render_template("index.html")


@app.route("/exit")
def exit_route():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()
    return "Server shutting down..."


@io.on("exit")
def exit_io():
    """Exit the session."""
    print("Server shutting down...")
    io.stop()


@io.on("atoms:request")
def atoms_request(data):
    """Return the atoms."""
    print(f"atoms:request {data = }")
    if "filename" in app.config:
        for idx, atoms_dict in enumerate(get_atomsdict_list(app.config["filename"])):
            emit("atoms:upload", atoms_dict)
            # At some point we should just emit all messages and call this function again for the remainder ?
            # if idx > 1000:
            #     break
    else:
        emit("atoms:upload", {})


@io.on("modifier:schema")
def modifier_schema():
    config = GlobalConfig.load()

    for modifier in config.modify_functions:
        module_name, function_name = modifier.rsplit(".", 1)
        module = importlib.import_module(module_name)
        modifier_cls = getattr(module, function_name)
        schema = modifier_cls.model_json_schema()

        if modifier in config.function_schema:
            kwargs = config.function_schema[modifier]
            for key, value in kwargs.items():
                schema["properties"][key]["default"] = value

        io.emit(
            "modifier:schema",
            {"name": modifier, "schema": schema},
        )


@io.on("modifier:run")
def modifier_run(data):
    import ase

    points = np.array([[val["x"], val["y"], val["z"]] for val in data["points"]])
    segments = np.array(data["segments"])

    if "atoms" in data:
        atoms = atoms_from_json(data["atoms"])
    else:
        atoms = ase.Atoms()

    module_name, function_name = data["name"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    modifier_cls = getattr(module, function_name)
    modifier = modifier_cls(**data["params"])
    # available_methods = {x.__name__: x for x in [Explode, Duplicate]}

    # modifier = available_methods[data["name"]](**data["params"])
    print(f"modifier:run {modifier = } from {data['params'] = }")
    atoms_list = modifier.run(
        atom_ids=data["selection"], atoms=atoms, points=points, segments=segments
    )
    io.emit("atoms:clear", int(data["step"]) + 1)
    for idx, atoms in tqdm.tqdm(enumerate(atoms_list)):
        atoms_dict = atoms_to_json(atoms)
        io.emit("atoms:upload", {idx + 1 + int(data["step"]): atoms_dict})

    io.emit("view:set", int(data["step"]) + 1)
    io.emit("view:play")


@io.on("analysis:schema")
def analysis_schema(data):
    config = GlobalConfig.load()
    if "atoms" not in data:
        return
    atoms = atoms_from_json(data["atoms"])

    for modifier in config.analysis_functions:
        module_name, function_name = modifier.rsplit(".", 1)
        module = importlib.import_module(module_name)
        modifier_cls = getattr(module, function_name)
        schema = modifier_cls.schema_from_atoms(atoms)

        io.emit(
            "analysis:schema",
            {"name": modifier, "schema": schema},
        )


@io.on("analysis:run")
def analysis_run(data):
    # print(f"analysis:run {data = }")
    atoms_list = [atoms_from_json(x) for x in data["atoms_list"].values()]

    print(f"Analysing {len(atoms_list)} frames")

    module_name, function_name = data["name"].rsplit(".", 1)
    module = importlib.import_module(module_name)
    cls = getattr(module, function_name)
    instance = cls(**data["params"])

    fig = instance.run(atoms_list, data["selection"])
    return fig.to_json()


@io.on("config")
def config(data):
    pass


# @app.route("/download/<int:idx>")
# def download(idx):
#     socketio.emit("atoms:download", [idx])
#     return "OK"

# @app.route("/download")
# def download():
#     from zndraw.view import ZnDraw
#     import ase.io

#     url = request.base_url
#     url = url.replace("/download", "")
#     print(f"Downloading {url = }")

#     zndraw = ZnDraw(url=url)
#     print(f"Downloading {zndraw = }")

#     file = BytesIO()
#     ase.io.write(file, zndraw[0], format="xyz")
#     file.seek(0)
#     return send_file(file, mimetype="text/plain", as_attachment=True)


@io.on("download")
def download(data):
    atoms = [atoms_from_json(x) for x in data["atoms_list"].values()]
    import ase.io

    file = StringIO()
    ase.io.write(file, atoms, format="xyz")
    file.seek(0)
    return file.read()


@io.on("atoms:download")
def atoms_download(data):
    """Return the atoms."""
    emit("atoms:download", data, broadcast=True, include_self=False)


@io.on("atoms:upload")
def atoms_upload(data):
    """Return the atoms."""
    emit("atoms:upload", data, broadcast=True, include_self=False)


@io.on("view:set")
def display(data):
    """Display the atoms at the given index"""
    emit("view:set", data, broadcast=True, include_self=False)


@io.on("upload")
def upload(data):
    from io import StringIO

    import ase.io
    import tqdm

    # tested with small files only

    format = data["filename"].split(".")[-1]
    if format == "h5":
        print("H5MD format not supported for uploading yet")
        # import znh5md
        # stream = BytesIO(data["content"].encode("utf-8"))
        # atoms = znh5md.ASEH5MD(stream).get_atoms_list()
        # for idx, atoms in tqdm.tqdm(enumerate(atoms)):
        #     atoms_dict = atoms_to_json(atoms)
        #     io.emit("atoms:upload", {idx: atoms_dict})
    else:
        stream = StringIO(data["content"])
        io.emit("atoms:clear", 0)
        for idx, atoms in tqdm.tqdm(enumerate(ase.io.iread(stream, format=format))):
            atoms_dict = atoms_to_json(atoms)
            io.emit("atoms:upload", {idx: atoms_dict})
    emit("view:set", 0)
