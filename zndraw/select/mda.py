from pydantic import BaseModel, Field
import ase.io
import MDAnalysis as mda
import io

class SelectionBase(BaseModel):
    selection: str = Field(..., description="MDAnalysis selection string")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        with io.StringIO() as f:
            ase.io.write(f, atoms, format='xyz')
            u = mda.Universe(f, format='XYZ', in_memory=True)