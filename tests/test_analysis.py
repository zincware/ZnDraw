import znsocket

from zndraw import ZnDraw


def run_queue(vis, key, msg: dict):
    modifier_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key=f"queue:{vis.token}:{key}",
    )
    modifier_queue.update(msg)
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(5)


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
