import numpy as np
import znsocket

from zndraw import ZnDraw
from zndraw.analyse import Properties1D


def run_queue(vis, key, msg: dict):
    modifier_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key=f"queue:{vis.token}:{key}",
    )
    modifier_queue.update(msg)
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(7)


def test_run_analysis_distance(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.selection = [0, 1]

    run_queue(vis, "analysis", {"Distance": {}})

    fig = vis.figures["Distance"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"


def test_run_analysis_Properties1D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}

    run_queue(vis, "analysis", {"Properties1D": {"value": "energy"}})

    fig = vis.figures["Properties1D"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"
    assert fig.layout.yaxis.title.text == "energy"


def test_run_analysis_Properties2D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}

    run_queue(
        vis,
        "analysis",
        {"Properties2D": {"x_data": "energy", "y_data": "step", "color": "energy"}},
    )

    fig = vis.figures["Properties2D"]
    # assert that the x-axis label is "step"
    assert fig.layout.yaxis.title.text == "step"
    assert fig.layout.xaxis.title.text == "energy"


def test_run_analysis_DihedralAngle(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.selection = [0, 1, 3, 4]

    run_queue(vis, "analysis", {"DihedralAngle": {}})

    fig = vis.figures["DihedralAngle"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"


def test_analysis_Properties1D_json_schema(s22_energy_forces):
    # add custom info keys
    for atoms in s22_energy_forces:
        atoms.info["custom"] = 42
        atoms.info["custom2"] = np.random.rand(10)
        atoms.info["custom3"] = np.random.rand(10, 5)
        atoms.calc.results["custom4"] = np.random.rand(10)
        atoms.arrays["arr"] = np.zeros_like(atoms.get_positions())

    schema = Properties1D.model_json_schema_from_atoms(s22_energy_forces[0])
    assert set(schema["properties"]["value"]["enum"]) == {
        "energy",
        "forces",
        "custom",
        "numbers",
        "positions",
        "arr",
        "custom2",
        "custom3",
        "custom4",
    }
