import dataclasses
import importlib
import pathlib
import typing

import ase.io
import pydantic
import tqdm
import znh5md
from pydantic import BaseModel, Field, PrivateAttr


class Config(BaseModel):
    file: str = Field(..., description="Trajectory File")
    bond_size: float = Field(1.0, description="Bond Size")
    sphere_size: float = Field(1.0, description="Sphere Size")
    resolution: int = Field(10, description="Resolution")
    max_fps: int = Field(100, description="Maximum Frames Per Second")
    frames_per_post: int = Field(100, description="Frames per JS POST request")
    active_update_function: str = Field(None, description="Active Update Function")
    material: str = Field("MeshPhongMaterial", description="Material")
    antialias: bool = Field(True, description="Antialias")

    _update_functions = PrivateAttr(default_factory=dict)
    _atoms_cache = PrivateAttr(default_factory=dict)

    def add_update_function(self, update_function) -> dict:
        """Add an update function to the config.

        Arguments
        ---------
        update_function: str
            The name of the update function 'module.cls'

        Returns
        -------
        dict: The signature of the function.
        """
        module_name, function_name = update_function.rsplit(".", 1)
        if function_name in self._update_functions:
            raise ValueError(f"Function {function_name} already exists.")

        self.active_update_function = function_name

        module = importlib.import_module(module_name)
        instance: pydantic.BaseModel = getattr(module, function_name)()
        self._update_functions[function_name] = instance

        return self._update_functions[function_name].schema()

    def set_update_function_parameter(self, value):
        """Set a parameter of the update function."""
        instance = self._update_functions[value["function_id"]]
        attribute = value["property"].lower()
        value = value["value"]
        if instance.__annotations__[attribute] == float:
            value = float(value)
        elif instance.__annotations__[attribute] == int:
            value = int(value)
        elif instance.__annotations__[attribute] == bool:
            value = bool(value)
        else:
            value = value
        setattr(instance, attribute, value)

    def apply_update_function(self, selected_ids, step, **kwargs):
        """Apply the update function to the selected atoms.

        Arguments
        ---------
        selected_ids: list
            The selected ids.
        step: int
            The current step.

        Returns
        -------
        Save the updated atoms in the cache.
        """
        if self.active_update_function is None:
            return

        step = int(step)
        function = self._update_functions[self.active_update_function].run

        atoms = function(selected_ids, self.get_atoms(step=step), **kwargs)

        for key in list(self._atoms_cache.keys()):
            # we remove all steps after the current one
            if key > step:
                del self._atoms_cache[key]

        for idx, atom in enumerate(atoms):
            self._atoms_cache[idx + step + 1] = atom

    def get_atoms(self, step: int = 0) -> ase.Atoms:
        """Get the atoms for a given step.

        Raises:
            KeyError: If the step could not be loaded.
        """
        try:
            return self._atoms_cache[step]
        except KeyError:
            self.load_atoms(step)
        return self._atoms_cache[step]

    def load_atoms(self, step: int = 999999999):
        """Load the atoms up to a given step."""
        # TODO ZnH5MD
        if self.active_update_function is not None:
            return  # We don't want to load any further from file at this point
        if pathlib.Path(self.file).suffix == ".h5":
            # We load all at once here
            self._atoms_cache.update(
                dict(enumerate(znh5md.ASEH5MD(self.file).get_atoms_list()[: step + 1]))
            )
        else:
            for idx, atoms in enumerate(tqdm.tqdm(ase.io.iread(self.file))):
                self._atoms_cache[idx] = atoms
                if idx == step:
                    break


config: Config = None

graph = None
atoms = None
