import enum
import io
import typing as t

import ase.io
import MDAnalysis as mda
import numpy as np
import plotly.express as px
from MDAnalysis.analysis import msd, rdf
from pydantic import BaseModel, Field

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class MDAInterRDFNorm(str, enum.Enum):
    """Normalization of RDF."""

    rdf = "rdf"
    density = "density"
    none = "none"


class MDAInterRDF(BaseModel):
    """Select Particles using MDAnalysis selection syntax."""

    discriminator: t.Literal["MDAInterRDF"] = Field("MDAInterRDF")
    selection_a: str = Field(..., description="MDAnalysis selection string")
    selection_b: str = Field(..., description="MDAnalysis selection string")
    nbins: int = Field(100, description="Number of bins")
    range: str = Field("0, 10", description="Range of histogram")
    norm: MDAInterRDFNorm = Field(MDAInterRDFNorm.rdf, description="Normalization of RDF")

    def run(self, vis) -> None:
        self.range = tuple(map(float, self.range.split(",")))
        with io.StringIO() as f:
            ase.io.write(f, list(vis), format="xyz")
            u = mda.Universe(f, format="XYZ", in_memory=True)

        grp1 = u.select_atoms(self.selection_a)
        grp2 = u.select_atoms(self.selection_b)

        box = np.diag(vis[vis.step].cell)
        u.dimensions = [box[0], box[1], box[2], 90, 90, 90]

        result = (
            rdf.InterRDF(
                grp1, grp2, nbins=self.nbins, range=self.range, norm=self.norm.value
            )
            .run()
            .results
        )

        fig = px.line(
            x=result.bins,
            y=result.rdf,
            # title="Distance between selected particles",
            render_mode="svg",  # This is important, otherwise openGL will be used
            # and there can/will be issues with three.js
        )

        vis.figure = fig.to_json()


class MDAEinsteinMSD(BaseModel):
    discriminator: t.Literal["MDAEinsteinMSD"] = Field("MDAEinsteinMSD")
    selection: str = Field("all", description="MDAnalysis selection string")
    delta_t: float = Field(0.1, description="Time step")

    def run(self, vis: "ZnDraw") -> None:
        vis.log("Running EinsteinMSD - Loading trajectory...")
        with io.StringIO() as f:
            ase.io.write(f, list(vis), format="xyz")
            u = mda.Universe(f, format="XYZ", in_memory=True)
        box = np.diag(vis[vis.step].cell)
        u.dimensions = [box[0], box[1], box[2], 90, 90, 90]
        # u.trajectory.dt = self.delta_t

        vis.log("Running EinsteinMSD - Running analysis...")
        MSD = msd.EinsteinMSD(u, select=self.selection, msd_type="xyz", fft=True)
        MSD.run()

        nframes = MSD.n_frames
        lagtimes = np.arange(nframes) * self.delta_t

        fig = px.line(
            x=lagtimes,
            y=MSD.results.timeseries,
            render_mode="svg",  # This is important, otherwise openGL will be used
            # and there can/will be issues with three.js
        )
        vis.figure = fig.to_json()


from zndraw.settings import _ANALYSIS_FUNCTIONS

_ANALYSIS_FUNCTIONS.append("zndraw.analyse.mda.MDAInterRDF")
_ANALYSIS_FUNCTIONS.append("zndraw.analyse.mda.MDAEinsteinMSD")
