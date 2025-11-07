from zndraw import ZnDraw
import rdkit2ase
from datetime import datetime
import uuid

vis = ZnDraw(
  url="http://localhost:5000/",
  room=str(uuid.uuid4()),
  user="user-fb9f0d37"
)

etoh = rdkit2ase.smiles2conformers("CCO", numConfs=1000)

start = datetime.now()
vis.extend(etoh)
end = datetime.now()
print(f"Uploaded {len(etoh)} ethanol molecules in {end - start}")

print(vis[0].__dict__)
print(etoh[0].__dict__)