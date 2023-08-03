import dataclasses
import enum
import importlib
import io
import pathlib
import typing
import uuid

import ase.io
import pydantic
import tqdm
import znh5md
from pydantic import BaseModel, Field, PrivateAttr

from zndraw.settings import GlobalConfig

if pathlib.Path("~/.zincware/zndraw/config.json").expanduser().exists():
    SETTINGS = GlobalConfig.from_file()
else:
    SETTINGS = GlobalConfig()


class CameraChoices(str, enum.Enum):
    OrthographicCamera = "OrthographicCamera"
    PerspectiveCamera = "PerspectiveCamera"


class Config(BaseModel):
    file: str = Field(..., description="Trajectory File")
    remote: str = Field(None, description="Remote to use for loading data via ZnTrack.")
    rev: str = Field(None, description="Revision to use for loading data via ZnTrack.")
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
    analysis_functions: typing.List[str] = SETTINGS.analysis_functions
    modify_functions: typing.List[str] = SETTINGS.modify_functions
    selection_functions: typing.List[str] = SETTINGS.selection_functions
    bonds_functions: typing.List[str] = SETTINGS.bonds_functions
    js_frame_buffer: tuple = Field(
        (250, 250), description="Javascript frame buffer in negative/positive direction"
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
        if update_function in SETTINGS.function_schema:
            kwargs = SETTINGS.function_schema[update_function]
            for key, value in kwargs.items():
                schema["properties"][key]["default"] = value
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

    def export_selection(self, step, selected_ids):
        string_io = io.StringIO()
        atoms = self.get_atoms(step=step)
        atoms = atoms[selected_ids]
        ase.io.write(string_io, atoms, format="xyz")
        byte_string = string_io.getvalue().encode("utf-8")
        bytes_io = io.BytesIO(byte_string)
        return bytes_io

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
        if self.remote is not None and self.rev is not None:
            from zndraw.tools.zntrack_data import get_atoms_via_dvc

            self._atoms_cache.update(
                dict(enumerate(get_atoms_via_dvc(self.file, self.remote, self.rev)))
            )
        else:
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
streaming: uuid.UUID = None
