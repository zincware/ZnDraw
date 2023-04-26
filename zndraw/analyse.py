import itertools

import ase
import numpy as np
import pandas as pd
import plotly.express as px
from pydantic import BaseModel

from zndraw import globals


class Distance(BaseModel):
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
        # smooth_df = df.rolling(window=100).mean().dropna()
        # for col in smooth_df.columns:
        #     if col != "step":
        #         fig.add_scatter(x=smooth_df["step"], y=smooth_df[col], name=f"smooth_{col}")
        return fig


class Properties2D(BaseModel):
    horizontal: str = "CV1"
    vertical: str = "CV2"
    color: str = "ω1"

    def run(self, ids):
        atoms_lst = list(globals.config._atoms_cache.values())

        cv1 = [x.calc.results[self.horizontal] for x in atoms_lst]
        cv2 = [x.calc.results[self.vertical] for x in atoms_lst]
        omega1 = [x.calc.results[self.color] for x in atoms_lst]

        df = pd.DataFrame({"cv1": cv1, "cv2": cv2, "omega1": omega1})
        fig = px.scatter(df, x="cv1", y="cv2", color="omega1", render_mode="svg")
        fig.update_yaxes(
            scaleanchor="x",
            scaleratio=1,
        )

        return fig
