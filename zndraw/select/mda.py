import io
import typing as t

import ase.io
import MDAnalysis as mda
from pydantic import BaseModel, Field


class MDAnalysis(BaseModel):
    """Select Particles using MDAnalysis selection syntax."""

    discriminator: t.Literal["MDAnalysis"] = Field("MDAnalysis")
    selection: str = Field(..., description="MDAnalysis selection string")
    append: bool = Field(False, description="Append to current selection")

    def run(self, vis) -> None:
        with io.StringIO() as f:
            ase.io.write(f, vis[vis.step], format="xyz")
            u = mda.Universe(f, format="XYZ", in_memory=True)

        selection = u.select_atoms(self.selection).ids.tolist()
        selection = [i - 1 for i in selection]
        if self.append:
            vis.selection = list(set(vis.selection + selection))
        else:
            vis.selection = selection


from zndraw.settings import _SELECTION_FUNCTIONS

_SELECTION_FUNCTIONS.append("zndraw.select.mda.MDAnalysis")
