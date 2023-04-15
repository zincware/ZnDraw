import dataclasses
import importlib
import ase.io


@dataclasses.dataclass
class Config:
    file: str = None
    animate: bool = None
    sphere_size: float = None
    bond_size: float = None
    max_fps: int = None
    update_function: str = None
    frames_per_post: int = 10
    restart_animation: bool = False
    resolution: int = 5
    repeat: tuple = (1, 1, 1)
    frame_buffer: int = 100

    def get_update_function(self):
        if self.update_function is None:
            return None
        module_name, function_name = self.update_function.rsplit(".", 1)
        module = importlib.import_module(module_name)
        return getattr(module, function_name)

    def get_atoms(self, step=0) -> ase.Atoms:
        try:
            return _atoms_cache[step]
        except KeyError:
            if not self.animate and step != 0:
                raise
            for idx, atoms in enumerate(ase.io.iread(self.file)):
                _atoms_cache[idx] = atoms.copy().repeat(self.repeat)
                if step == 0:
                    break
        return _atoms_cache[step]


# TODO set defaults here and load in typer?

_atoms_cache: dict = {}
config = Config()

graph = None
atoms = None
