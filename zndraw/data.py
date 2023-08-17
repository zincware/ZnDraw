from dask.distributed import Lock, Variable, Client
import dataclasses
import ase
import znh5md
import pandas as pd
import dask.dataframe as dd
import numpy as np
from ase.data.colors import jmol_colors

def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)

def _get_radius(value):
    return 0.25 * (2 - np.exp(-0.2 * value)),

@dataclasses.dataclass
class DataHandler:
    client: Client

    def create_dataset(self, filename):
        # TODO this should happen on a worker
        atoms = znh5md.ASEH5MD(filename).get_atoms_list()
        df = dd.DataFrame.from_dict({"atoms": atoms}, npartitions=10)
        self.client.persist(df)
        self.client.publish_dataset(atoms=df)

    def __len__(self):
        return len(self.get_dataset())

    def get_dataset(self):
        return self.client.get_dataset("atoms") 

    def get_atoms(self, _slice) -> dict[str, ase.Atoms]:
        df = self.get_dataset()
        return df.loc[_slice]["atoms"].compute().to_dict()
    
    def get_atoms_json(self, _slice) -> dict[str, ase.Atoms]:
        data = {}
        for idx, atoms in self.get_atoms(_slice).items():
            atoms_dict = atoms.todict()
            for key in list(atoms_dict):
                # includes ['numbers', 'positions', 'cell', 'pbc']
                if isinstance(atoms_dict[key], np.ndarray):
                    atoms_dict[key] = atoms_dict[key].tolist()
            
            atoms_dict["colors"] = [_rgb2hex(jmol_colors[number]) for number in atoms_dict["numbers"]]
            atoms_dict["radii"] = [_get_radius(number) for number in atoms_dict["numbers"]]
            data[idx] = atoms_dict
        
        return data

