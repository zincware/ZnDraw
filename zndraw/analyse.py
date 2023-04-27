import itertools

import ase
import numpy as np
import pandas as pd
import plotly.express as px
from pydantic import BaseModel

from zndraw import globals


class Distance(BaseModel):
    smooth: bool = False

    @classmethod
    def schema_from_atoms(cls, atoms):
        return cls.schema()

    def run(self, ids):
        atoms_lst = list(globals.config._atoms_cache.values())
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
    x_data: str = "step"
    y_data: str = "energy"
    color: str = "energy"
    fix_aspect_ratio: bool = True

    @classmethod
    def schema_from_atoms(cls, atoms):
        schema = cls.schema()
        available_properties = list(atoms[0].calc.results.keys())
        available_properties += ["step"]
        schema["properties"]["x_data"]["enum"] = available_properties
        schema["properties"]["y_data"]["enum"] = available_properties
        schema["properties"]["color"]["enum"] = available_properties
        return schema

    def run(self, ids):
        print(f"run {self}")
        atoms_lst = list(globals.config._atoms_cache.values())

        if self.x_data == "step":
            x_data = list(range(len(atoms_lst)))
        else:
            x_data = [x.calc.results[self.x_data] for x in atoms_lst]

        if self.y_data == "step":
            y_data = list(range(len(atoms_lst)))
        else:
            y_data = [x.calc.results[self.y_data] for x in atoms_lst]

        if self.color == "step":
            color = list(range(len(atoms_lst)))
        else:
            color = [x.calc.results[self.color] for x in atoms_lst]

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
    value: str = "energy"
    smooth: bool = False

    @classmethod
    def schema_from_atoms(cls, atoms):
        schema = cls.schema()
        available_properties = list(atoms[0].calc.results.keys())
        schema["properties"]["value"]["enum"] = available_properties
        return schema

    def run(self, ids):
        atoms_lst = list(globals.config._atoms_cache.values())

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
