import itertools
import typing as t
from typing import Any

import ase
import numpy as np
import pandas as pd
import plotly.express as px
from pydantic import BaseModel, Field

from zndraw.utils import set_global_atoms


class Distance(BaseModel):
    method: t.Literal["Distance"] = "Distance"

    smooth: bool = False

    def run(self, atoms_lst, ids):
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
            render_mode="svg"  # This is important, otherwise openGL will be used
            # and there can/will be issues with three.js
        )
        if self.smooth:
            smooth_df = df.rolling(window=100).mean().dropna()
            for col in smooth_df.columns:
                if col != "step":
                    fig.add_scatter(
                        x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}"
                    )
        return fig


class Properties2D(BaseModel):
    method: t.Literal["Properties2D"] = "Properties2D"

    x_data: str = "step"
    y_data: str = "energy"
    color: str = "energy"
    fix_aspect_ratio: bool = True

    @classmethod
    def model_json_schema(cls, *args, **kwargs) -> dict[str, Any]:
        schema = super().model_json_schema(*args, **kwargs)
        print(f"GATHERING PROPERTIES FROM {ATOMS=}")  # noqa: F821
        try:
            available_properties = list(ATOMS.calc.results)  # noqa: F821
            available_properties += list(ATOMS.arrays)  # noqa: F821
            available_properties += ["step"]
            schema["properties"]["x_data"]["enum"] = available_properties
            schema["properties"]["y_data"]["enum"] = available_properties
            schema["properties"]["color"]["enum"] = available_properties
        except AttributeError:
            pass
        return schema

    def run(self, atoms_lst, ids):
        print(f"run {self}")

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
        return fig


class Properties1D(BaseModel):
    method: t.Literal["Properties1D"] = "Properties1D"

    value: str = "energy"
    smooth: bool = False

    @classmethod
    def model_json_schema(cls, *args, **kwargs) -> dict[str, Any]:
        schema = super().model_json_schema(*args, **kwargs)
        try:
            available_properties = list(
                ATOMS.calc.results.keys()  # noqa: F821
            )  # global ATOMS object
            print(f"AVAILABLE PROPERTIES: {available_properties=}")
            schema["properties"]["value"]["enum"] = available_properties
        except AttributeError:
            pass
        return schema

    def run(self, atoms_lst, ids):
        data = np.array([x.calc.results[self.value] for x in atoms_lst])

        df = pd.DataFrame({"step": list(range(len(atoms_lst))), self.value: data})

        fig = px.line(df, x="step", y=self.value, render_mode="svg")

        if self.smooth:
            smooth_df = df.rolling(window=100).mean().dropna()
            for col in smooth_df.columns:
                if col != "step":
                    fig.add_scatter(
                        x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}"
                    )

        return fig


def get_analysis_class(methods):
    class Analysis(BaseModel):
        method: methods = Field(
            ..., description="Analysis method", discriminator="method"
        )

        def run(self, *args, **kwargs) -> list[ase.Atoms]:
            return self.method.run(*args, **kwargs)

        @classmethod
        def model_json_schema_from_atoms(
            cls, atoms, *args, **kwargs
        ) -> dict[str, t.Any]:
            with set_global_atoms(atoms):
                result = cls.model_json_schema(*args, **kwargs)
            return result

        @classmethod
        def model_json_schema(cls, *args, **kwargs) -> dict[str, t.Any]:
            schema = super().model_json_schema(*args, **kwargs)
            for prop in [x.__name__ for x in t.get_args(methods)]:
                schema["$defs"][prop]["properties"]["method"]["options"] = {
                    "hidden": True
                }
                schema["$defs"][prop]["properties"]["method"]["type"] = "string"

            return schema

    return Analysis
