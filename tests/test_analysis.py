from zndraw import ZnDraw


def test_run_analysis_distance(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.selection = [0, 1]

    vis.socket.emit("analysis:run", {"method": {"discriminator": "Distance"}})
    vis.socket.sleep(7)
    fig = vis.figures["Distance"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"


def test_run_analysis_Properties1D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.socket.emit(
        "analysis:run", {"method": {"discriminator": "Properties1D", "value": "energy"}}
    )
    vis.socket.sleep(7)
    fig = vis.figures["Properties1D"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"
    assert fig.layout.yaxis.title.text == "energy"


def test_run_analysis_Properties2D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.socket.emit(
        "analysis:run",
        {
            "method": {
                "discriminator": "Properties2D",
                "x_data": "energy",
                "y_data": "step",
                "color": "energy",
            }
        },
    )
    vis.socket.sleep(7)
    fig = vis.figures["Properties2D"]
    # assert that the x-axis label is "step"
    assert fig.layout.yaxis.title.text == "step"
    assert fig.layout.xaxis.title.text == "energy"


def test_run_analysis_DihedralAngle(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figures == {}
    vis.selection = [0, 1, 3, 4]

    vis.socket.emit("analysis:run", {"method": {"discriminator": "DihedralAngle"}})
    vis.socket.sleep(7)
    fig = vis.figures["DihedralAngle"]
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"
