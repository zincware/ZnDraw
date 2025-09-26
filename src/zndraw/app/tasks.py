from celery import shared_task
import znh5md
import ase
from ase.calculators.singlepoint import SinglePointCalculator
from tqdm import tqdm

def atoms_to_dict(atoms: ase.Atoms) -> dict:
    if not atoms.calc:
        return atoms.todict()
    results = atoms.todict()
    results["<SinglePointCalculator>"] = atoms.calc.results
    return results

def atoms_from_dict(d: dict) -> ase.Atoms:
    if "<SinglePointCalculator>" not in d:
        return ase.Atoms.fromdict(d)
    calc_results = d.pop("<SinglePointCalculator>")
    atoms = ase.Atoms.fromdict(d)
    atoms.calc = SinglePointCalculator(atoms)
    atoms.calc.results = calc_results
    return atoms

@shared_task(ignore_result=True)
def read_file() -> None:
    from zndraw import Client

    client = Client(room="testroom", url="http://localhost:5000")
    client.connect()
    io = znh5md.IO("/Users/fzills/tools/zndraw-communication-testing/structures.h5")
    frames = io[:]
    for atoms in tqdm(frames, desc="Uploading frames"):
        client.append(atoms_to_dict(atoms))
