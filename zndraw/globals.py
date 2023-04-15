import dataclasses
import importlib


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


# TODO set defaults here and load in typer?


config = Config()

graph = None
atoms = None
