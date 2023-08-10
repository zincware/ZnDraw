from dask.distributed import Lock, Variable, Client
import dataclasses
import ase
import znh5md
import pandas as pd
import dask.dataframe as dd
import numpy as np


@dataclasses.dataclass
class DataHandler:
    client: Client

    def create_dataset(self, filename):
        atoms = znh5md.ASEH5MD(filename).get_atoms_list()
        df = dd.from_pandas(pd.DataFrame({"atoms": atoms}), npartitions=10)
        self.client.persist(df)
        self.client.publish_dataset(atoms=df)

    def get_dataset(self):
        return self.client.get_dataset("atoms") 

    def get_atoms(self, _slice) -> dict[str, ase.Atoms]:
        df = self.get_dataset()
        return df.loc[_slice]["atoms"].compute().to_dict()
    
    def get_atoms_json(self, _slice) -> dict[str, ase.Atoms]:
        _dict = {key: val.todict() for key, val in self.get_atoms(_slice).items()}
        for frame in _dict:
            for key in _dict[frame]:
                if isinstance(_dict[frame][key], np.ndarray):
                    _dict[frame][key] = _dict[frame][key].tolist()
        return _dict
