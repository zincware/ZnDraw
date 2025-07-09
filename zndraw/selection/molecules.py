import rdkit2ase
from zndraw.base import Extension

class SelectSmarts(Extension):
    smarts: str
    
    def run(self, vis: "ZnDraw", **kwargs) -> None:
        box = vis.atoms
        selection = rdkit2ase.match_substructure(
            atoms=box,
            pattern=self.smarts,
        )
        # flatten the selection
        selection = [item for sublist in selection for item in sublist]
        vis.selection = selection
        
