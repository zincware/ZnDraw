import dataclasses
import logging
import typing as t
from abc import abstractmethod

from pydantic import BaseModel

if t.TYPE_CHECKING:
    from zndraw import ZnDraw

log = logging.getLogger(__name__)


class Extension(BaseModel):
    @abstractmethod
    def run(self, vis: "ZnDraw", **kwargs) -> None:
        raise NotImplementedError("run method must be implemented in subclass")


@dataclasses.dataclass  # TODO: move to a separate file, so it can be imported in other files
class FileIO:
    name: str | None = None
    start: int = 0
    stop: int | None = None
    step: int = 1
    remote: str | None = None
    rev: str | None = None
    convert_nan: bool = False

    def to_dict(self):
        return dataclasses.asdict(self)
