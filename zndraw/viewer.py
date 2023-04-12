import networkx as nx
import plotly.graph_objs as go
import dash
from dash import dcc
from dash import html
import numpy as np
from dash.dependencies import Input, Output
import dataclasses


def get_max_vals(graph: nx.Graph):
    x_max = 0
    y_max = 0
    z_max = 0
    for node in graph.nodes:
        x_max = max(max(x_max, graph.nodes[node]["x"]), x_max)
        y_max = max(max(y_max, graph.nodes[node]["y"]), y_max)
        z_max = max(max(z_max, graph.nodes[node]["z"]), z_max)

    return x_max, y_max, z_max


@dataclasses.dataclass
class DashApp:
    graph: nx.Graph
    fig: go.Figure = dataclasses.field(default_factory=go.Figure)
    app: dash.Dash = dataclasses.field(default_factory=dash.Dash)

    select_canvas: bool = True
    select_atoms: bool = True

    _atoms_scatter: go.Scatter3d = None
    _canvas_surface: go.Surface = None  # there seems to be a copy of this in fig.data

    def update_layout(self):
        self.fig.update_layout(showlegend=False)
        self.fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
            )
        )

    def plot_atoms(self):
        x = []
        y = []
        z = []
        number = []

        for node in self.graph.nodes:
            x.append(self.graph.nodes[node]["x"])
            y.append(self.graph.nodes[node]["y"])
            z.append(self.graph.nodes[node]["z"])
            number.append(self.graph.nodes[node]["number"])

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        number = np.array(number)

        self._atoms_scatter = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="markers",
            marker={"color": number, "size": 10},
            hoverinfo=None if self.select_atoms else "skip",
        )

        self.fig.add_trace(self._atoms_scatter)

    def plot_bonds(self):
        traces = []
        for edge in self.graph.edges:
            traces.append(
                go.Scatter3d(
                    x=[self.graph.nodes[edge[0]]["x"], self.graph.nodes[edge[1]]["x"]],
                    y=[self.graph.nodes[edge[0]]["y"], self.graph.nodes[edge[1]]["y"]],
                    z=[self.graph.nodes[edge[0]]["z"], self.graph.nodes[edge[1]]["z"]],
                    mode="lines",
                    line={"width": 10, "color": "black"},
                    hoverinfo="skip",
                )
            )

        self.fig.add_traces(traces)

    def add_canvas_slider(self):
        if not self.select_canvas:
            return
        x_max, y_max, z_max = get_max_vals(self.graph)

        x = np.linspace(0, x_max, 10)
        y = np.linspace(0, y_max, 10)
        z = np.ones((10, 10))

        self._canvas_surface = go.Surface(
            x=x,
            y=y,
            z=z,
            opacity=0.3,
            surfacecolor=z * [0],
            showscale=False,
            contours={
                "x": {
                    "show": True,
                    "color": "white",
                    "size": x_max / 10,
                    "start": 0,
                    "end": x_max,
                },
                "y": {
                    "show": True,
                    "color": "white",
                    "size": y_max / 10,
                    "start": 0,
                    "end": y_max,
                },
            },
        )

        self.fig.add_trace(self._canvas_surface)

    def callback_graph_clickData(self):
        @self.app.callback(
            Output("graph", "figure", allow_duplicate=True),
            [Input("graph", "clickData")],
            prevent_initial_call=True,
        )
        def add_scatter(hoverData):
            if hoverData is not None:
                self.fig.add_trace(
                    go.Scatter3d(
                        x=[hoverData["points"][0]["x"]],
                        y=[hoverData["points"][0]["y"]],
                        z=[hoverData["points"][0]["z"]],
                        mode="markers",
                        marker={"color": "black", "size": 20},
                        hoverinfo="skip",
                    )
                )
            self.fig.update_layout(uirevision="constant")
            return self.fig

    def callback_slider_canvas(self):
        @self.app.callback(
            Output("graph", "figure", allow_duplicate=True),
            [Input("slider-canvas", "value")],
            prevent_initial_call=True,
        )
        def move_canvas(value):
            # self._canvas_surface["z"] = np.ones_like(self._canvas_surface["z"]) * value
            self.fig.data[0]["z"] = np.ones_like(self.fig.data[0]["z"]) * value
            self.fig.update_layout(uirevision="constant")
            return self.fig

    def run_dash_server_click(self):
        x_max, y_max, z_max = get_max_vals(self.graph)

        div = [
            dcc.Graph(
                figure=self.fig, style={"width": "90vw", "height": "90vh"}, id="graph"
            ),
        ]

        if self.select_canvas:
            div.append(
                dcc.Slider(0, z_max, id="slider-canvas", updatemode="drag"),
            )
            self.callback_slider_canvas()
        self.app.layout = html.Div(div)
        self.callback_graph_clickData()
        self.app.run_server(debug=True, use_reloader=False)

    def run_dash_server(self):
        self.app.layout = html.Div(
            [dcc.Graph(figure=self.fig, style={"width": "90vw", "height": "90vh"})]
        )

        self.app.run_server(debug=True, use_reloader=True)
