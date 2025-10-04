import random
import typing as t

import networkx as nx
import numpy as np
from pydantic import Field
import pandas as pd
import plotly.express as px

from zndraw.extensions.abc import Extension, ExtensionType

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Analysis(Extension):
    """The base class for all analysis extensions."""

    category: t.ClassVar[ExtensionType] = ExtensionType.ANALYSIS

class DihedralAngle(Analysis):
    def run(self, vis: "ZnDraw") -> None:
        atoms_lst = vis[:]
        dihedral_angles = []

        selection = list(vis.selection)

        if len(selection) != 4:
            raise ValueError("Please select exactly 4 atoms")
        for atoms in atoms_lst:
            dihedral_angles.append(
                atoms.get_dihedrals(indices=[selection], mic=True)[0]
            )
        df = pd.DataFrame(
            {"step": list(range(len(atoms_lst))), "dihedral": dihedral_angles}
        )
        fig = px.line(df, x="step", y="dihedral", render_mode="svg")

        meta_step = np.arange(len(atoms_lst))

        fig.update_traces(
            customdata=np.stack([meta_step], axis=-1),
        )
        # update_figure_layout(fig)

        vis.figures["DihedralAngle"] = fig


analysis: dict[str, t.Type[Analysis]] = {
    "DihedralAngle": DihedralAngle,
}