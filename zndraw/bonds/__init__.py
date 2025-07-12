import ase
import networkx as nx
import numpy as np
from networkx.exception import NetworkXError
from pydantic import BaseModel, Field
from rdkit2ase import ase2networkx


class ASEComputeBonds(BaseModel):
    def get_bonds(self, atoms: ase.Atoms, graph: nx.Graph = None):
        """Generate bonds and store connectivity in atoms.info"""
        graph = ase2networkx(atoms, pbc=False)

        # Extract connectivity from networkx graph
        connectivity = []
        for u, v, data in graph.edges(data=True):
            bond_order = data.get("bond_order", None)
            connectivity.append((u, v, bond_order))

        # Store connectivity in atoms.info instead of atoms.connectivity
        atoms.info["connectivity"] = connectivity

        return connectivity
