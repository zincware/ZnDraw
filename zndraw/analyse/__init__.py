import itertools
import logging
import typing as t

import ase
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pydantic import Field

from zndraw.base import Extension

try:
    from zndraw.analyse import mda  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass


log = logging.getLogger(__name__)


def _get_data_from_frames(key, frames: list[ase.Atoms]):
    if frames[0].calc is not None and key in frames[0].calc.results:
        data = np.array([x.calc.results[key] for x in frames])
    elif key in frames[0].arrays:
        data = np.array([x.arrays[key] for x in frames])
    elif key in frames[0].info:
        data = np.array([x.info[key] for x in frames])
    else:
        raise ValueError(f"Property '{key}' not found in atoms")

    return data


class AnaylsisMethod(Extension):
    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        return cls.model_json_schema()


class DihedralAngle(AnaylsisMethod):
    def run(self, vis):
        atoms_lst = list(vis)
        dihedral_angles = []

        if len(vis.selection) != 4:
            raise ValueError("Please select exactly 4 atoms")
        for atoms in atoms_lst:
            dihedral_angles.append(
                atoms.get_dihedrals(indices=[vis.selection], mic=True)[0]
            )
        df = pd.DataFrame(
            {"step": list(range(len(atoms_lst))), "dihedral": dihedral_angles}
        )
        fig = px.line(df, x="step", y="dihedral", render_mode="svg")

        meta_step = np.arange(len(atoms_lst))
        # meta_idx = np.full_like(meta_step, np.nan)

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
        )

        vis.figures.update({"DihedralAngle": fig})


class Distance(AnaylsisMethod):
    smooth: bool = False
    mic: bool = True

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()

        schema["properties"]["smooth"]["format"] = "checkbox"
        schema["properties"]["mic"]["format"] = "checkbox"

        return schema

    def run(self, vis):
        atoms_lst, ids = list(vis), vis.selection
        distances = {}
        for x in itertools.combinations(ids, 2):
            distances[f"{tuple(x)}"] = []
        for atoms in atoms_lst:
            for x in itertools.combinations(ids, 2):
                distances[f"{tuple(x)}"].append(
                    atoms.get_distance(x[0], x[1], mic=self.mic)
                )

        df = pd.DataFrame({"step": list(range(len(atoms_lst)))} | distances)

        fig = px.line(
            df,
            x="step",
            y=df.columns,
            title="Distance between selected particles",
            render_mode="svg",  # This is important, otherwise openGL will be used
            # and there can/will be issues with three.js
        )
        if self.smooth:
            smooth_df = df.rolling(window=100).mean().dropna()
            for col in smooth_df.columns:
                if col != "step":
                    fig.add_scatter(
                        x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}"
                    )
        meta_step = np.arange(len(atoms_lst))
        # meta_idx = np.full_like(meta_step, np.nan)

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
        )

        vis.figures.update({"Distance": fig})


class Properties2D(AnaylsisMethod):
    x_data: str
    y_data: str
    color: str
    fix_aspect_ratio: bool = True

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()

        available_properties = list(atoms.arrays.keys())
        available_properties += list(atoms.info.keys())
        if atoms.calc is not None:
            available_properties += list(atoms.calc.results.keys())  # global ATOMS object

        available_properties += ["step"]
        schema["properties"]["x_data"]["enum"] = available_properties
        schema["properties"]["y_data"]["enum"] = available_properties
        schema["properties"]["color"]["enum"] = available_properties
        schema["properties"]["fix_aspect_ratio"]["format"] = "checkbox"

        return schema

    def run(self, vis):
        atoms_lst = list(vis)
        log.info(f"run {self}")

        if self.x_data == "step":
            x_data = list(range(len(atoms_lst)))
        else:
            x_data = _get_data_from_frames(self.x_data, atoms_lst)

        if self.y_data == "step":
            y_data = list(range(len(atoms_lst)))
        else:
            y_data = _get_data_from_frames(self.y_data, atoms_lst)

        if self.color == "step":
            color = list(range(len(atoms_lst)))
        else:
            color = _get_data_from_frames(self.color, atoms_lst)

        y_data = np.array(y_data).reshape(-1)
        x_data = np.array(x_data).reshape(-1)
        color = np.array(color).reshape(-1)

        df = pd.DataFrame({self.x_data: x_data, self.y_data: y_data, self.color: color})
        fig = px.scatter(df, x=self.x_data, y=self.y_data, color=self.color)
        if self.fix_aspect_ratio:
            fig.update_yaxes(
                scaleanchor="x",
                scaleratio=1,
            )

        meta_step = np.arange(len(atoms_lst))
        # meta_idx = np.full_like(meta_step, np.nan)

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
        )

        vis.figures.update({"Properties2D": fig})


class ForceCorrelation(AnaylsisMethod):
    """Compute the correlation between two properties for the current frame."""

    x_data: str
    y_data: str

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()

        available_properties = list(atoms.arrays.keys())
        available_properties += list(atoms.info.keys())
        if atoms.calc is not None:
            available_properties += list(atoms.calc.results.keys())
        schema["properties"]["x_data"]["enum"] = available_properties
        schema["properties"]["y_data"]["enum"] = available_properties

        return schema

    def run(self, vis):
        atoms = vis.atoms
        x_data = _get_data_from_frames(self.x_data, [atoms])
        y_data = _get_data_from_frames(self.y_data, [atoms])

        x_data = np.linalg.norm(x_data, axis=-1)
        y_data = np.linalg.norm(y_data, axis=-1)

        vis.log(f"x_data: {x_data.shape}, y_data: {y_data.shape}")

        x_data = x_data.reshape(-1)
        y_data = y_data.reshape(-1)

        meta_step = [None for _ in range(len(x_data))]
        meta_idx = list(range(len(x_data)))

        df = pd.DataFrame(
            {
                self.x_data: x_data,
                self.y_data: y_data,
            }
        )

        fig = px.scatter(df, x=self.x_data, y=self.y_data, render_mode="svg")
        fig.update_traces(customdata=np.stack([meta_step, meta_idx], axis=-1))

        vis.figures.update({"ForceCorrelation": fig})


class Properties1D(AnaylsisMethod):
    value: str
    smooth: bool = False

    aggregation: t.Literal["mean", "median", "max", ""] = Field(
        "",
        description="For multidimensional data, aggregate over all dimensions, except the first one.",
    )

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()

        available_properties = list(atoms.arrays.keys())
        available_properties += list(atoms.info.keys())
        if atoms.calc is not None:
            available_properties += list(atoms.calc.results.keys())
        log.critical(f"AVAILABLE PROPERTIES: {available_properties=}")
        schema["properties"]["value"]["enum"] = available_properties

        return schema

    def run(self, vis):
        vis.log("Downloading data...")
        atoms_lst = list(vis)
        data = _get_data_from_frames(self.value, atoms_lst)

        if data.ndim > 1:
            axis = tuple(range(1, data.ndim))
            if self.aggregation == "mean":
                data = np.mean(data, axis=axis)
            elif self.aggregation == "median":
                data = np.median(data, axis=axis)
            elif self.aggregation == "max":
                data = np.max(data, axis=axis)

        df = pd.DataFrame({"step": list(range(len(atoms_lst))), self.value: data})

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=df["step"],
                y=df[self.value],
                mode="lines+markers",
            )
        )
        # set xlabel to be step
        fig.update_layout(
            xaxis_title="step",
            yaxis_title=self.value,
        )

        if self.smooth:
            smooth_df = df.rolling(window=100).mean().dropna()
            for col in smooth_df.columns:
                if col != "step":
                    fig.add_scatter(
                        x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}"
                    )

        meta_step = np.arange(len(atoms_lst))
        # meta_idx = np.full_like(meta_step, np.nan)

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
        )

        vis.figures.update({"Properties1D": fig})


analyses: dict[str, t.Type[AnaylsisMethod]] = {
    Properties1D.__name__: Properties1D,
    DihedralAngle.__name__: DihedralAngle,
    Distance.__name__: Distance,
    Properties2D.__name__: Properties2D,
    ForceCorrelation.__name__: ForceCorrelation,
}
