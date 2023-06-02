import dataclasses
import enum
import importlib
import io
import pathlib
import typing

import ase.io
import pydantic
import tqdm
import znh5md
from pydantic import BaseModel, Field, PrivateAttr

_ANALYSIS_FUNCTIONS = [
    "zndraw.analyse.Properties1D",
    "zndraw.analyse.Properties2D",
    "zndraw.analyse.Distance",
]
_MODIFY_FUNCTIONS = [
    "zndraw.examples.Explode",
    "zndraw.examples.Delete",
    "zndraw.examples.Move",
    "zndraw.examples.Duplicate",
    "zndraw.examples.AddLineParticles",
]

_BONDS_FUNCTIONS = [
    "zndraw.tools.data.ASEComputeBonds",
]
_SELECTION_FUNCTIONS = []


class CameraChoices(str, enum.Enum):
    OrthographicCamera = "OrthographicCamera"
    PerspectiveCamera = "PerspectiveCamera"


class Config(BaseModel):
    file: str = Field(..., description="Trajectory File")
    camera: CameraChoices = Field(
        CameraChoices.PerspectiveCamera, description="Camera style"
    )
    bond_size: float = Field(1.0, description="Bond Size")
    sphere_size: float = Field(1.0, description="Sphere Size")
    resolution: int = Field(10, description="Resolution")
    max_fps: int = Field(100, description="Maximum Frames Per Second")
    frames_per_post: int = Field(100, description="Frames per JS POST request")
    active_update_function: str = Field(None, description="Active Update Function")
    material: str = Field("MeshPhongMaterial", description="Material")
    material_wireframe: bool = Field(False, description="Material Wireframe")
    antialias: bool = Field(True, description="Antialias")
    continuous_loading: bool = Field(
        False, description="Continuous Loading of the trajectory"
    )
    auto_restart: bool = Field(False, description="Auto Restart")
    analysis_functions: typing.List[str] = _ANALYSIS_FUNCTIONS
    modify_functions: typing.List[str] = _MODIFY_FUNCTIONS
    selection_functions: typing.List[str] = _SELECTION_FUNCTIONS
    bonds_functions: typing.List[str] = _BONDS_FUNCTIONS
    js_frame_buffer: tuple = Field(
        (50, 50), description="Javascript frame buffer in negative/positive direction"
    )
    _atoms_cache = PrivateAttr(default_factory=dict)
    _loaded_modifiers: typing.Dict[str, typing.Any] = PrivateAttr(default_factory=dict)
    _modifier_applied: bool = PrivateAttr(False)

    def get_modifier_schema(self, update_function) -> dict:
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
        module = importlib.import_module(module_name)
        instance: pydantic.BaseModel = getattr(module, function_name)()
        schema = instance.schema()
        schema["title"] = update_function
        return schema

    def run_modifier(self, modifier, selected_ids, step, modifier_kwargs, **kwargs):
        """Run the update function.

        Arguments
        ---------
        update_function: str
            The name of the update function 'module.cls'
        selected_ids: list
            The selected ids.
        step: int
            The current step.

        Returns
        -------
        Save the updated atoms in the cache.
        """
        self._modifier_applied = True
        if modifier not in self._loaded_modifiers:
            module_name, function_name = modifier.rsplit(".", 1)
            module = importlib.import_module(module_name)
            cls: pydantic.BaseModel = getattr(module, function_name)
            instance = cls(**modifier_kwargs)
            self._loaded_modifiers[modifier] = instance
        instance = self._loaded_modifiers[modifier]
        kwargs = kwargs | modifier_kwargs
        atoms = instance.run(selected_ids, self.get_atoms(step=step).copy(), **kwargs)
        for key in list(self._atoms_cache.keys()):
            # we remove all steps after the current one
            if key > step:
                del self._atoms_cache[key]

        for idx, atom in enumerate(atoms):
            self._atoms_cache[idx + step + 1] = atom

    def reset_scene_modifiers(self) -> None:
        """Reset the scene modifiers."""
        for key in list(self._loaded_modifiers):
            del self._loaded_modifiers[key]
        self._loaded_modifiers = {}

    def export_atoms(self):
        file = io.BytesIO()
        db = znh5md.io.DataWriter(file)
        db.initialize_database_groups()
        db.add(znh5md.io.AtomsReader(self.atoms_list))
        return file

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

    @property
    def atoms_list(self) -> typing.List[ase.Atoms]:
        """Get a list of the atoms in the atoms cache."""
        return list(self._atoms_cache.values())

    def load_atoms(self, step: int = 999999999):
        """Load the atoms up to a given step."""
        # TODO ZnH5MD
        if self._modifier_applied:
            return  # We don't want to load any further from file at this point
        # print("Loading atoms")
        if pathlib.Path(self.file).suffix == ".h5":
            # We load all at once here
            self._atoms_cache.update(
                dict(enumerate(znh5md.ASEH5MD(self.file).get_atoms_list()))
            )
        else:
            for idx, atoms in enumerate(tqdm.tqdm(ase.io.iread(self.file))):
                self._atoms_cache[idx] = atoms


config: Config = None

bond_method = None
