import dataclasses
import enum
import functools
import importlib
import io
import logging
import pathlib
import time
import typing

import ase.io
import pydantic
import tqdm
import znh5md
from pydantic import BaseModel, Field, PrivateAttr

log = logging.getLogger(__name__)

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
_SELECTION_FUNCTIONS = []


class CameraChoices(str, enum.Enum):
    OrthographicCamera = "OrthographicCamera"
    PerspectiveCamera = "PerspectiveCamera"


class LoadingState(str, enum.Enum):
    NotLoaded = "NotLoaded"
    Loading = "Loading"
    Loaded = "Loaded"


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
    antialias: bool = Field(True, description="Antialias")
    continuous_loading: bool = Field(
        False, description="Continuous Loading of the trajectory"
    )
    analysis_functions: typing.List[str] = _ANALYSIS_FUNCTIONS
    modify_functions: typing.List[str] = _MODIFY_FUNCTIONS
    selection_functions: typing.List[str] = _SELECTION_FUNCTIONS

    _atoms_cache = PrivateAttr(default_factory=list)
    _modifier_applied: bool = PrivateAttr(False)
    _loaded: LoadingState = PrivateAttr(LoadingState.NotLoaded)

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
        module_name, function_name = modifier.rsplit(".", 1)
        module = importlib.import_module(module_name)
        cls: pydantic.BaseModel = getattr(module, function_name)
        instance = cls(**modifier_kwargs)
        atoms = instance.run(selected_ids, self.atoms_list[step].copy(), **kwargs)

        self._atoms_cache = self._atoms_cache[: step + 1]
        self._atoms_cache.extend(atoms)

    def export_atoms(self):
        file = io.BytesIO()
        db = znh5md.io.DataWriter(file)
        db.initialize_database_groups()
        db.add(znh5md.io.AtomsReader(self.atoms_list))
        return file

    @property
    def atoms_list(self) -> typing.List[ase.Atoms]:
        """Get a list of the atoms in the atoms cache."""
        return self._atoms_cache

    def load_atoms(self):
        """Load the atoms up to a given step."""
        # TODO ZnH5MD
        while self._loaded == LoadingState.Loading:
            time.sleep(0.1)
            log.debug("Waiting for loading to finish")
        self._loaded = LoadingState.Loading
        log.debug("Loading atoms")
        if self._modifier_applied:
            return  # We don't want to load any further from file at this point
        if pathlib.Path(self.file).suffix == ".h5":
            # We load all at once here
            self._atoms_cache = znh5md.ASEH5MD(self.file).get_atoms_list()
        else:
            atoms = list(ase.io.iread(self.file))
            self._atoms_cache = atoms
        self._loaded = LoadingState.Loaded


config: Config = None
