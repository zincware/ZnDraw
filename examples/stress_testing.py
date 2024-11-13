# Run this test with and without eventlet.monkey_patch()
import eventlet
eventlet.monkey_patch()

import tqdm
from rdkit2ase import smiles2conformers
import uuid
import datetime # noqa

#### Import ZnDraw ####
from zndraw import ZnDraw

vis = ZnDraw(url="https://zndraw.icp.uni-stuttgart.de/", token=uuid.uuid4().hex)
#### ------------- ####

# vis._refresh_client.delay_between_calls = datetime.timedelta(milliseconds=10)

conformers = smiles2conformers("CCCCCCCCCO", numConfs=1000)

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
