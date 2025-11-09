import typing as t

import ase
import numpy as np
import vesin
from ase.neighborlist import natural_cutoffs

T_CONNECTIVITY = t.List[t.Tuple[int, int, int | float | None]]


def add_connectivity(atoms: ase.Atoms, scale: float = 1.2) -> None:
    """Add connectivity information to an ASE Atoms object.

    This function determines atomic connectivity based on scaled natural covalent radii.
    A bond is considered to exist between two atoms if their distance is less than
    the sum of their scaled covalent radii. The bond order is set to 1 by default.

    The connectivity information is stored in the `atoms.info['connectivity']`
    attribute as a list of tuples, where each tuple is (atom_i, atom_j, bond_order).
    Bond order can be 1 (single), 2 (double), 3 (triple), 1.5 (aromatic), or None.

    Parameters
    ----------
    atoms
        The ASE Atoms object to analyze.
    scale
        A scaling factor for the covalent radii to determine the cutoff distance.
        A larger value results in more bonds being detected.

    Returns
    -------
    None
        The function modifies the Atoms object in-place.
    """
    if len(atoms) == 0:
        return
    # Calculate scaled covalent radii for each atom
    atom_radii = np.array(natural_cutoffs(atoms, mult=scale))

    # Create a matrix of pairwise cutoff distances.
    # The cutoff for a pair (i, j) is the sum of their scaled radii.
    pairwise_cutoffs = atom_radii[:, None] + atom_radii[None, :]
    max_cutoff = np.max(pairwise_cutoffs)

    # Store original PBC settings and temporarily disable them to avoid
    # bonds being drawn across periodic boundaries
    original_pbc = atoms.pbc.copy()
    atoms.set_pbc(False)

    # Get the neighbor list using the maximum possible cutoff.
    # 'i' and 'j' are indices, 'd' is the distance.
    i_list, j_list, d_list, _ = vesin.ase_neighbor_list(
        "ijdS", atoms, cutoff=max_cutoff, self_interaction=False
    )

    # Restore original PBC settings
    atoms.set_pbc(original_pbc)

    connectivity: T_CONNECTIVITY = []

    # Iterate through all neighbor pairs found by vesin
    for i, j, d in zip(i_list, j_list, d_list):
        # To avoid double counting (e.g., adding both (0,1) and (1,0)), we only
        # consider pairs where the first index is smaller than the second.
        if i < j:
            # Check if the distance is within the specific cutoff for this atom pair
            if d < pairwise_cutoffs[i, j]:
                # Bond is found, add it to the list with bond order 1
                connectivity.append((i, j, 1))

    # there is a bug, that connectivity must be an array, otherwise it fails!
    atoms.info["connectivity"] = np.array(connectivity, dtype=np.int32)
