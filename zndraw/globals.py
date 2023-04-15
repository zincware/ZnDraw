import dataclasses
import importlib
import functools
import ase.io


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None
    sphere_size: float = None
    bond_size: float = None
    max_fps: int = None
    update_function: str = None

    def get_update_function(self):
        if self.update_function is None:
            return None
        module_name, function_name = self.update_function.rsplit(".", 1)
        module = importlib.import_module(module_name)
        return getattr(module, function_name)

    def get_atoms(self, step=0) -> ase.Atoms:
        if step in _atoms_cache:
            return _atoms_cache[step]
        for idx, atoms in enumerate(ase.io.iread(self.file)):
            _atoms_cache[idx] = atoms.copy()
            if idx == step:
                return atoms


# TODO set defaults here and load in typer?

_atoms_cache: dict = {}
config = Config()

graph = None
atoms = None
