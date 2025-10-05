import typing as t
from collections.abc import MutableMapping

import plotly.graph_objs as go
import plotly.io as pio

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Figures(MutableMapping):
    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def _dict_to_figure(self, figure: dict) -> go.Figure:
        if not isinstance(figure, dict):
            raise ValueError("Figure data must be a dictionary")

        if figure.get("type") == "plotly":
            fig = pio.from_json(figure["data"])
            return fig

        raise ValueError("Unsupported figure type or invalid figure data")

    def _figure_to_dict(self, figure: go.Figure) -> dict:
        if not isinstance(figure, go.Figure):
            raise ValueError("Only plotly.graph_objs.Figure instances are supported")
        response = figure.to_json()
        if response is None:
            raise ValueError("Failed to convert figure to JSON")
        return {"data": response, "type": "plotly"}

    def __getitem__(self, key: str) -> go.Figure:
        if key not in self.vis._figures:
            response = self.vis.api.get_figure(key=key)
            if response is None:
                raise KeyError(f"Figure with key '{key}' does not exist")
            self.vis._figures[key] = response
        if key not in self.vis._figures:
            raise KeyError(f"Figure with key '{key}' does not exist")
        return self._dict_to_figure(figure=self.vis._figures[key])

    def __setitem__(self, key: str, value: t.Any) -> None:
        self.vis._figures[key] = self._figure_to_dict(figure=value)
        self.vis.api.add_figure(
            key=key,
            figure=self.vis._figures[key],
        )

    def __delitem__(self, key: str) -> None:
        self.vis._figures.pop(key, None)
        self.vis.api.delete_figure(key=key)

    def __iter__(self):
        return iter(self.vis.api.list_figures())

    def __len__(self) -> int:
        return len(self.vis.api.list_figures())
