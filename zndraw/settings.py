import pathlib
import typing

import pydantic

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


class GlobalConfig(pydantic.BaseModel):
    analysis_functions: typing.List[str] = _ANALYSIS_FUNCTIONS
    modify_functions: typing.List[str] = _MODIFY_FUNCTIONS
    bonds_functions: typing.List[str] = _BONDS_FUNCTIONS
    selection_functions: typing.List[str] = _SELECTION_FUNCTIONS
    function_schema: typing.Dict[str, dict] = {}

    def save(self, path="~/.zincware/zndraw/config.json"):
        save_path = pathlib.Path(path).expanduser()
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, "w") as f:
            f.write(self.json())

    @classmethod
    def from_file(cls, path="~/.zincware/zndraw/config.json"):
        load_path = pathlib.Path(path).expanduser()
        with open(load_path, "r") as f:
            return cls.parse_raw(f.read())
