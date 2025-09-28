from celery import shared_task
import znh5md
from tqdm import tqdm
from zndraw.utils import update_colors_and_radii, atoms_to_dict


@shared_task
def read_file() -> None:
    from zndraw import Client

    client = Client(room="testroom", url="http://localhost:5000")
    client.connect()
    io = znh5md.IO("/Users/fzills/tools/zndraw-communication-testing/structures.h5")
    frames = io[:]
    for atoms in tqdm(frames, desc="Uploading frames"):
        update_colors_and_radii(atoms)
        # atoms = atoms.repeat((4, 4, 4))
        client.append(atoms_to_dict(atoms))
