import itertools
import logging
import typing as t

import numpy as np
import pandas as pd
import plotly.express as px
from pydantic import ConfigDict, Field

from zndraw.base import Extension, MethodsCollection

try:
    from zndraw.analyse import mda  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass


log = logging.getLogger(__name__)


def _schema_from_atoms(schema, cls):
    return cls.model_json_schema_from_atoms(schema)


class DihedralAngle(Extension):
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
        vis.figure = fig.to_json()


class Distance(Extension):
    smooth: bool = False

    def run(self, vis):
        atoms_lst, ids = list(vis), vis.selection
        distances = {}
        for x in itertools.combinations(ids, 2):
            distances[f"{tuple(x)}"] = []
        for atoms in atoms_lst:
            positions = atoms.get_positions()
            for x in itertools.combinations(ids, 2):
                distances[f"{tuple(x)}"].append(
                    np.linalg.norm(positions[x[0]] - positions[x[1]])
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
        vis.figure = fig.to_json()


class Properties2D(Extension):
    x_data: str
    y_data: str
    color: str
    fix_aspect_ratio: bool = True

    model_config = ConfigDict(json_schema_extra=_schema_from_atoms)

    @classmethod
    def model_json_schema_from_atoms(cls, schema: dict) -> dict:
        ATOMS = cls.get_atoms()
        log.debug(f"GATHERING PROPERTIES FROM {ATOMS=}")
        try:
            available_properties = list(ATOMS.calc.results)
            available_properties += list(ATOMS.arrays)
            available_properties += ["step"]
            schema["properties"]["x_data"]["enum"] = available_properties
            schema["properties"]["y_data"]["enum"] = available_properties
            schema["properties"]["color"]["enum"] = available_properties
        except AttributeError:
            pass
        return schema

    def run(self, vis):
        atoms_lst = list(vis)
        log.info(f"run {self}")

        if self.x_data == "step":
            x_data = list(range(len(atoms_lst)))
        else:
            try:
                x_data = [x.calc.results[self.x_data] for x in atoms_lst]
            except KeyError:
                x_data = [x.arrays[self.x_data] for x in atoms_lst]

        if self.y_data == "step":
            y_data = list(range(len(atoms_lst)))
        else:
            try:
                y_data = [x.calc.results[self.y_data] for x in atoms_lst]
            except KeyError:
                y_data = [x.arrays[self.y_data] for x in atoms_lst]

        if self.color == "step":
            color = list(range(len(atoms_lst)))
        else:
            try:
                color = [x.calc.results[self.color] for x in atoms_lst]
            except KeyError:
                color = [x.arrays[self.color] for x in atoms_lst]

        y_data = np.array(y_data).reshape(-1)
        x_data = np.array(x_data).reshape(-1)
        color = np.array(color).reshape(-1)

        df = pd.DataFrame({self.x_data: x_data, self.y_data: y_data, self.color: color})
        fig = px.scatter(
            df, x=self.x_data, y=self.y_data, color=self.color, render_mode="svg"
        )
        if self.fix_aspect_ratio:
            fig.update_yaxes(
                scaleanchor="x",
                scaleratio=1,
            )
        vis.figure = fig.to_json()


class Properties1D(Extension):
    value: str
    smooth: bool = False

    model_config = ConfigDict(json_schema_extra=_schema_from_atoms)
    aggregation: t.Literal["mean", "median", "max", ""] = Field(
        "",
        description="For multidimensional data, aggregate over all dimensions, except the first one.",
    )

    @classmethod
    def model_json_schema_from_atoms(cls, schema: dict) -> dict:
        ATOMS = cls.get_atoms()
        try:
            available_properties = list(ATOMS.calc.results.keys())  # global ATOMS object
            log.debug(f"AVAILABLE PROPERTIES: {available_properties=}")
            schema["properties"]["value"]["enum"] = available_properties
        except AttributeError:
            print(f"{ATOMS=}")
        return schema

    def run(self, vis):
        vis.log("Downloading data...")
        atoms_lst = list(vis)
        data = np.array([x.calc.results[self.value] for x in atoms_lst])

        if data.ndim > 1:
            axis = tuple(range(1, data.ndim))
            if self.aggregation == "mean":
                data = np.mean(data, axis=axis)
            elif self.aggregation == "median":
                data = np.median(data, axis=axis)
            elif self.aggregation == "max":
                data = np.max(data, axis=axis)

        df = pd.DataFrame({"step": list(range(len(atoms_lst))), self.value: data})

        fig = px.line(df, x="step", y=self.value, render_mode="svg")

        if self.smooth:
            smooth_df = df.rolling(window=100).mean().dropna()
            for col in smooth_df.columns:
                if col != "step":
                    fig.add_scatter(
                        x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}"
                    )

        vis.figure = fig.to_json()


methods = t.Union[
    Properties1D,
    DihedralAngle,
    Distance,
    Properties2D,
]


class Analysis(MethodsCollection):
    method: methods = Field(description="Analysis method", discriminator="discriminator")
