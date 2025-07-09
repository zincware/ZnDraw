import rdkit2ase
from zndraw.base import Extension

class CreateBox(Extension):
    smiles: list[str]
    count: list[int]
    density: float = 1000
    packmol: str = "packmol.jl"
    

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        if len(self.smiles) != len(self.count):
            raise ValueError("The length of smiles and count must be the same.")
        frames = []
        for smiles in self.smiles:
            frames.append(
                rdkit2ase.smiles2conformers(
                    smiles=smiles,
                    numConfs=1
                )
            )
        
        box = rdkit2ase.pack(
            data=frames,
            counts=self.count,
            density=self.density,
            packmol=self.packmol,
        )
        vis.append(box)
