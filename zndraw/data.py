from dask.distributed import Lock, Variable, Client
import dataclasses
import ase
import znh5md
import pandas as pd
import dask.dataframe as dd


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

    def get_atoms(self, idx) -> ase.Atoms:
        df = self.get_dataset()
        return df.loc[idx:idx+100]["atoms"].values.compute()
