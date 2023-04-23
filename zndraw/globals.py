import dataclasses
import importlib
import pathlib

import ase.io
import pydantic
import tqdm
import znh5md
from pydantic import BaseModel, Field
import tqdm


class Config(BaseModel):
    file: str = Field(..., description="Trajectory File")
    bond_size: float = Field(1.0, description="Bond Size")
    sphere_size: float = Field(1.0, description="Sphere Size")
    resolution: int = Field(5, description="Resolution")

    def get_atoms(self, step: int = 0):
        try:
            return _atoms_cache[step]
        except KeyError:
            # TODO ZnH5MD
            for idx, atoms in enumerate(tqdm.tqdm(ase.io.iread(self.file))):
                _atoms_cache[idx] = atoms
                if idx == step:
                    return atoms


# @dataclasses.dataclass
# class Config:
#     file: str = None
#     animate: bool = None
#     sphere_size: float = 1.0
#     bond_size: float = 1.0
#     max_fps: int = 100
#     update_function: str = None
#     frames_per_post: int = 100
#     restart_animation: bool = False
#     resolution: int = 5
#     repeat: tuple = (1, 1, 1)

#     _update_function_name: str = None

#     @property
#     def update_function_name(self):
#         if self._update_function_name is not None:
#             return self._update_function_name
#         return self.update_function.rsplit(".", 1)[1]

#     @update_function_name.setter
#     def update_function_name(self, value):
#         self._update_function_name = value

#     def get_update_signature(self):
#         if self.update_function is None:
#             return {}
#         module_name, function_name = self.update_function.rsplit(".", 1)
#         if function_name in _update_functions:
#             return _update_functions[self.update_function_name].schema()
#         module = importlib.import_module(module_name)
#         instance: pydantic.BaseModel = getattr(module, function_name)()
#         _update_functions[function_name] = instance
#         return _update_functions[self.update_function_name].schema()

#     def get_update_function(self):
#         module_name, function_name = self.update_function.rsplit(".", 1)
#         if self.update_function is None:
#             return None
#         if function_name in _update_functions:
#             return _update_functions[self.update_function_name].run

#         module = importlib.import_module(module_name)
#         _update_functions[self.update_function_name] = getattr(module, function_name)()
#         return _update_functions[self.update_function_name].run

#     def set_update_function_parameters(self, value):
#         instance = _update_functions[value["function_id"]]
#         attribute = value["property"].lower()
#         value = value["value"]
#         if instance.__annotations__[attribute] == float:
#             value = float(value)
#         elif instance.__annotations__[attribute] == int:
#             value = int(value)
#         elif instance.__annotations__[attribute] == bool:
#             value = bool(value)
#         else:
#             value = value
#         setattr(instance, attribute, value)

#     def load_atoms(self, item=None):
#         if item == 0:
#             if pathlib.Path(self.file).suffix == ".h5":
#                 _atoms_cache[0] = znh5md.ASEH5MD(self.file)[0]
#             else:
#                 _atoms_cache[0] = ase.io.read(self.file)
#         elif self.update_function is not None:
#             return
#         else:
#             if pathlib.Path(self.file).suffix == ".h5":
#                 _atoms_cache.update(
#                     dict(enumerate(znh5md.ASEH5MD(self.file).get_atoms_list()))
#                 )
#             else:
#                 for idx, atom in enumerate(
#                     tqdm.tqdm(ase.io.iread(self.file), desc="File Reading")
#                 ):
#                     _atoms_cache[idx] = atom

#     def get_atoms(self, step=0) -> ase.Atoms:
#         try:
#             return _atoms_cache[step].repeat(self.repeat)
#         except KeyError:
#             if step != 0:
#                 raise
#             self.load_atoms(0)
#             return _atoms_cache[0].repeat(self.repeat)


# TODO set defaults here and load in typer?

_update_functions = {}

_atoms_cache: dict = {}
config = None

graph = None
atoms = None
