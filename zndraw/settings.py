import importlib
import json
import logging
import pathlib
import typing as t

import pydantic

log = logging.getLogger(__name__)

_ANALYSIS_FUNCTIONS = [
    "zndraw.analyse.Properties1D",
    "zndraw.analyse.Properties2D",
    "zndraw.analyse.Distance",
]
_MODIFY_FUNCTIONS = [
    "zndraw.modify.Explode",
    "zndraw.modify.Delete",
    "zndraw.modify.Move",
    "zndraw.modify.Duplicate",
    "zndraw.modify.AddLineParticles",
    "zndraw.modify.Rotate",
    "zndraw.modify.ChangeType",
    "zndraw.modify.Wrap",
    "zndraw.modify.Center",
    "zndraw.modify.Replicate",
    # "zndraw.modify.CustomModifier",
]

_BONDS_FUNCTIONS = [
    "zndraw.tools.data.ASEComputeBonds",
]
_SELECTION_FUNCTIONS = [
    "zndraw.select.NoneSelection",
    "zndraw.select.All",
    "zndraw.select.Invert",
    "zndraw.select.Range",
    "zndraw.select.Random",
    "zndraw.select.IdenticalSpecies",
    "zndraw.select.ConnectedParticles",
    "zndraw.select.Neighbour",
]


class GlobalConfig(pydantic.BaseModel):
    analysis_functions: t.List[str] = _ANALYSIS_FUNCTIONS
    modify_functions: t.List[str] = _MODIFY_FUNCTIONS
    bonds_functions: t.List[str] = _BONDS_FUNCTIONS
    selection_functions: t.List[str] = _SELECTION_FUNCTIONS
    function_schema: t.Dict[str, dict] = {}

    def save(self, path="~/.zincware/zndraw/config.json"):
        save_path = pathlib.Path(path).expanduser()
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, "w") as f:
            f.write(self.model_dump_json(indent=4))

    @classmethod
    def from_file(cls, path="~/.zincware/zndraw/config.json"):
        load_path = pathlib.Path(path).expanduser()
        with open(load_path, "r") as f:
            return cls(**json.load(f))

    @classmethod
    def load(cls):
        if pathlib.Path("~/.zincware/zndraw/config.json").expanduser().exists():
            return cls.from_file()
        else:
            return cls()

    def get_selection_methods(self):
        classes = []
        for method in self.selection_functions:
            module_name, cls_name = method.rsplit(".", 1)
            module = importlib.import_module(module_name)
            cls = getattr(module, cls_name)
            classes.append(cls)

        return t.Union[tuple(classes)]

    def get_analysis_methods(self):
        classes = []
        for method in self.analysis_functions:
            module_name, cls_name = method.rsplit(".", 1)
            module = importlib.import_module(module_name)
            cls = getattr(module, cls_name)
            classes.append(cls)

        return t.Union[tuple(classes)]

    def get_modify_methods(self, include: list = None):
        if include is None:
            classes = []
        else:
            classes = include
        for method in self.modify_functions:
            module_name, cls_name = method.rsplit(".", 1)
            try:
                module = importlib.import_module(module_name)
                cls = getattr(module, cls_name)
                classes.append(cls)
            except ModuleNotFoundError:
                log.critical(f"Module {module_name} not found - skipping")

        return t.Union[tuple(classes)]
