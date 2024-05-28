import logging
import typing as t

import pydantic

from zndraw.utils import ensure_path

log = logging.getLogger(__name__)


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
            "*.update_atoms": {"queue": "fast"},
            "*.get_selection_schema": {"queue": "fast"},
            "*.scene_schema": {"queue": "fast"},
            "*.geometries_schema": {"queue": "fast"},
            "*.analysis_schema": {"queue": "fast"},
            "*.handle_room_get": {"queue": "fast"},
            "*.activate_modifier": {"queue": "fast"},
            "*on_disconnect": {"queue": "fast"},
            "*upload_file": {"queue": "fast"},
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
