import dataclasses
import importlib
import pathlib

import ase.io
import tqdm
import znh5md


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None
    sphere_size: float = None
    bond_size: float = None
    max_fps: int = None
    update_function: str = None
    frames_per_post: int = 100
    restart_animation: bool = False
    resolution: int = 5
    repeat: tuple = (1, 1, 1)

    def get_update_function(self):
        if self.update_function is None:
            return None
        module_name, function_name = self.update_function.rsplit(".", 1)
        module = importlib.import_module(module_name)
        return getattr(module, function_name)

    def load_atoms(self, item=None):
        if item == 0:
            if pathlib.Path(self.file).suffix == ".h5":
                _atoms_cache[0] = znh5md.ASEH5MD(self.file)[0]
            else:
                _atoms_cache[0] = ase.io.read(self.file)
        elif self.update_function is not None:
            return
        else:
            if pathlib.Path(self.file).suffix == ".h5":
                _atoms_cache.update(
                    dict(enumerate(znh5md.ASEH5MD(self.file).get_atoms_list()))
                )
            else:
                for idx, atom in enumerate(
                    tqdm.tqdm(ase.io.iread(self.file), desc="File Reading")
                ):
                    _atoms_cache[idx] = atom

    def get_atoms(self, step=0) -> ase.Atoms:
        try:
            return _atoms_cache[step].repeat(self.repeat)
        except KeyError:
            if step != 0:
                raise
            self.load_atoms(0)
            return _atoms_cache[0]


# TODO set defaults here and load in typer?

_atoms_cache: dict = {}
config = Config()

graph = None
atoms = None
