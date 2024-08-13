from zndraw import ZnDraw
import plotly.io as pio


def test_run_analysis_distance(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figure is None
    vis.selection = [0, 1]

    vis.socket.emit("analysis:run", {"method": {"discriminator": "Distance"}})
    vis.socket.sleep(7)
    fig = pio.from_json(vis.figure)
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"
    

def test_run_analysis_Properties1D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figure is None
    vis.socket.emit("analysis:run", {"method": {"discriminator": "Properties1D", "value": "energy"}})
    vis.socket.sleep(7)
    fig = pio.from_json(vis.figure)
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"
    assert fig.layout.yaxis.title.text == "energy"


def test_run_analysis_Properties2D(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figure is None
    vis.socket.emit("analysis:run", {"method": {"discriminator": "Properties2D", "x_data": "energy", "y_data": "step", "color": "energy"}})
    vis.socket.sleep(7)
    fig = pio.from_json(vis.figure)
    # assert that the x-axis label is "step"
    assert fig.layout.yaxis.title.text == "step"
    assert fig.layout.xaxis.title.text == "energy"


def test_run_analysis_DihedralAngle(server, s22_energy_forces):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22_energy_forces)
    assert vis.figure is None
    vis.selection = [0, 1, 3, 4]

    vis.socket.emit("analysis:run", {"method": {"discriminator": "DihedralAngle"}})
    vis.socket.sleep(7)
    fig = pio.from_json(vis.figure)
    # assert that the x-axis label is "step"
    assert fig.layout.xaxis.title.text == "step"