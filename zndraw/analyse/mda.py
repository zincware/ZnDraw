from pydantic import BaseModel, Field
import ase.io
import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import io
import typing as t
import numpy as np

class MDAInterRDF(BaseModel):
    """Select Particles using MDAnalysis selection syntax."""

    discriminator: t.Literal["MDAInterRDF"] = Field("MDAInterRDF")
    selection_a: str = Field(..., description="MDAnalysis selection string")
    selection_b: str = Field(..., description="MDAnalysis selection string")
    nbins: int = Field(100, description="Number of bins")
    range: t.Tuple[float, float] = Field((0, 10), description="Range of histogram")
    norm: str = Field("rdf", description="Normalization of histogram")
    

    def run(self, vis) -> None:
        with io.StringIO() as f:
            ase.io.write(f, list(vis), format='xyz')
            u = mda.Universe(f, format='XYZ', in_memory=True)
        
        grp1 = u.select_atoms(self.selection_a)
        grp2 = u.select_atoms(self.selection_b)

        box = np.diag(vis[vis.step].cell)
        u.dimensions = [box[0], box[1], box[2], 90, 90, 90]

        result = rdf.InterRDF(grp1, grp2, nbins=self.nbins, range=self.range, norm=self.norm).run()

        print(result.rdf)


from zndraw.settings import _ANALYSIS_FUNCTIONS
_ANALYSIS_FUNCTIONS.append("zndraw.analyse.mda.MDAInterRDF")