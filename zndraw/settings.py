import importlib
import json
import logging
import pathlib
import typing as t

import pydantic

from zndraw.utils import ensure_path

log = logging.getLogger(__name__)

_ANALYSIS_FUNCTIONS = [
    "zndraw.analyse.Properties1D",
    "zndraw.analyse.Properties2D",
    "zndraw.analyse.Distance",
    "zndraw.analyse.DihedralAngle",
]
_MODIFY_FUNCTIONS = [
    "zndraw.modify.Delete",
    "zndraw.modify.Move",
    "zndraw.modify.Duplicate",
    "zndraw.modify.AddLineParticles",
    "zndraw.modify.Rotate",
    "zndraw.modify.ChangeType",
    "zndraw.modify.Wrap",
    "zndraw.modify.Center",
    "zndraw.modify.Replicate",
    "zndraw.modify.Connect",
    "zndraw.modify.ClearScene",
    # "zndraw.modify.CustomModifier",
]

_BONDS_FUNCTIONS = [
    "zndraw.tools.data.ASEComputeBonds",
]
_SELECTION_FUNCTIONS = [
    "zndraw.select.ConnectedParticles",
    "zndraw.select.NoneSelection",
    "zndraw.select.All",
    "zndraw.select.Invert",
    "zndraw.select.Range",
    "zndraw.select.Random",
    "zndraw.select.IdenticalSpecies",
    "zndraw.select.Neighbour",
    "zndraw.select.UpdateSelection",
]


class CacheSettings(pydantic.BaseModel):
    backend: str = "FileSystemCache"
    backend_options: dict = {}
    timeout: int = 60 * 60 * 24
    dir: str = "~/.zincware/zndraw/cache"
    threshold: int = 50000

    def to_dict(self):
        return dict(
            CACHE_TYPE=self.backend,
            CACHE_DEFAULT_TIMEOUT=self.timeout,
            CACHE_DIR=ensure_path(self.dir),
            CACHE_THRESHOLD=self.threshold,
        )


class CeleryBaseConfig(pydantic.BaseModel):
    name: t.Literal["CeleryBaseConfig"] = "CeleryBaseConfig"
    broker: str = "filesystem://"
    result_backend: str = "cache"
    cache_backend: str = "memory"
    task_ignore_result: bool = True

    @property
    def task_routes(self):
        return {
            "*._run_global_modifier": {"queue": "slow"},
            "*.update_atoms": {"queue": "io"},
            "*.get_selection_schema": {"queue": "io"},
            "*.scene_schema": {"queue": "io"},
            "*.geometries_schema": {"queue": "io"},
            "*.analysis_schema": {"queue": "io"},
            "*.handle_room_get": {"queue": "io"},
            "*.handle_room_set": {"queue": "io"},
            "*.activate_modifier": {"queue": "io"},
            "*on_disconnect": {"queue": "io"},
        }

    def to_dict(self):
        return dict(
            broker_url=self.broker,
            result_backend=self.result_backend,
            cache_backend=self.cache_backend,
            task_ignore_result=self.task_ignore_result,
            task_routes=self.task_routes,
        )


class CeleryFileSystemConfig(CeleryBaseConfig):
    name: t.Literal["CeleryFileSystemConfig"] = "CeleryFileSystemConfig"
    data_folder: str = "~/.zincware/zndraw/celery/out"
    data_folder_processed: str = "~/.zincware/zndraw/celery/processed"
    task_ignore_result: bool = True

    def to_dict(self):
        return dict(
            broker_url=self.broker,
            broker_transport_options=dict(
                data_folder_in=ensure_path(self.data_folder),
                data_folder_out=ensure_path(self.data_folder),
                data_folder_processed=ensure_path(self.data_folder_processed),
            ),
            result_backend=self.result_backend,
            cache_backend=self.cache_backend,
            task_ignore_result=self.task_ignore_result,
            task_routes=self.task_routes,
        )


class DatabaseSQLiteConfig(pydantic.BaseModel):
    name: t.Literal["DatabaseSQLiteConfig"] = "DatabaseSQLiteConfig"

    path: str = "~/.zincware/zndraw/database.sqlite"

    def get_path(self) -> str:
        path = pathlib.Path(self.path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        return f"sqlite:///{path}"


class DatabasePostgresConfig(pydantic.BaseModel):
    name: t.Literal["DatabasePostgresConfig"] = "DatabasePostgresConfig"

    host: str = "localhost"
    port: int = 5432
    database: str = "zndraw"
    user: str = "zndraw"
    password: str = "zndraw"

    def get_path(self) -> str:
        return f"postgresql://{self.user}:{self.password}@{self.host}:{self.port}/{self.database}"


class GlobalConfig(pydantic.BaseModel):
    cache: CacheSettings = pydantic.Field(default_factory=CacheSettings)
    celery: t.Union[CeleryBaseConfig, CeleryFileSystemConfig] = pydantic.Field(
        default_factory=CeleryFileSystemConfig, discriminator="name"
    )
    database: t.Union[DatabasePostgresConfig, DatabaseSQLiteConfig] = pydantic.Field(
        default_factory=DatabaseSQLiteConfig, discriminator="name"
    )

    # Socket settings
    read_batch_size: int = 1
    max_socket_data_size: int | float = 1024 * 8

    # Webclient Interface
    analysis_functions: t.List[str] = _ANALYSIS_FUNCTIONS
    modify_functions: t.List[str] = _MODIFY_FUNCTIONS
    bonds_functions: t.List[str] = _BONDS_FUNCTIONS
    selection_functions: t.List[str] = _SELECTION_FUNCTIONS
    function_schema: t.Dict[str, dict] = {}
    read_batch_size: int = 8

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
            try:
                module = importlib.import_module(module_name)
                cls = getattr(module, cls_name)
                classes.append(cls)
            except (ModuleNotFoundError, AttributeError):
                log.critical(f"Module {module_name} not found - skipping")

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
