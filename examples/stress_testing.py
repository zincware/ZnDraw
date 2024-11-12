from rdkit2ase import smiles2conformers
import tqdm
#### Import ZnDraw ####
from zndraw import ZnDraw

vis = ZnDraw(url="https://zndraw.icp.uni-stuttgart.de/", token="6b1d9f77")
#### ------------- ####

conformers = smiles2conformers("CCCCCCCCCO", numConfs=100)

# append
for atoms in tqdm.tqdm(conformers, desc="append", ncols=80):
    vis.append(atoms)

# read
for i in tqdm.trange(len(vis), desc="getitem", ncols=80):
    _ = vis[i]

# delete
for i in tqdm.tqdm(range(len(vis) - 1, -1, -1), desc="delete", ncols=80):
    del vis[i]

# extend
vis.extend(conformers)

# read_all
print(f"len(vis[:]): {len(vis[:])}")
