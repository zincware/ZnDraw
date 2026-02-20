"""Read-only filesystem provider via fsspec."""

import typing as t

from pydantic import BaseModel
from zndraw_joblib import Provider


class FileItem(BaseModel):
    """Normalized filesystem entry."""

    name: str
    path: str
    size: int = 0
    type: t.Literal["file", "directory"] = "file"


class FilesystemRead(Provider):
    """Read-only filesystem operations via fsspec."""

    category: t.ClassVar[str] = "filesystem"

    path: str = "/"
    glob: str | None = None

    def read(self, handler: t.Any) -> list[dict[str, t.Any]]:
        """List or glob files from the filesystem handler."""
        if self.glob:
            base = self.path.rstrip("/")
            pattern = f"{base}/{self.glob}" if base else self.glob
            items = [_from_info(handler.info(m)) for m in handler.glob(pattern)]
        else:
            items = [_from_info(info) for info in handler.ls(self.path, detail=True)]
        return [item.model_dump() for item in items]


def _from_info(info: dict[str, t.Any]) -> FileItem:
    """Normalize an fsspec info dict to a FileItem."""
    name: str = info.get("name", "")
    return FileItem(
        name=name.rsplit("/", 1)[-1],
        path=name,
        size=info.get("size", 0) or 0,
        type="directory" if info.get("type") == "directory" else "file",
    )
