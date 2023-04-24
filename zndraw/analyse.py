import itertools

import ase
import numpy as np
import pandas as pd
import plotly.express as px

from zndraw import globals


def get_distance_plot(step: int, ids: list):
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
        df.dropna().head(step),
        x="step",
        y=df.columns,
        title="Distance between selected particles",
    )
    return fig.to_json()
