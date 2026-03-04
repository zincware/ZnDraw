"""Accessor classes for ZnDraw client properties.

Thin wrappers implementing ``collections.abc`` interfaces over the REST API.
Each accessor takes an ``APIManager`` and delegates all I/O to it.
"""

from __future__ import annotations

import base64
import json
import time
from collections.abc import (
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, overload

from zndraw.geometries import geometries as geometry_models
from zndraw.geometries.base import BaseGeometry
from zndraw.geometries.camera import Camera
from zndraw.schemas import Preset, PresetApplyResult

if TYPE_CHECKING:
    import plotly.graph_objects as go

    from zndraw.client import APIManager
    from zndraw.schemas import MessageResponse


# =============================================================================
# ScreenshotImage
# =============================================================================


@dataclass
class ScreenshotImage:
    """Screenshot image with Jupyter display support.

    Parameters
    ----------
    data
        Raw PNG bytes of the captured screenshot.
    """

    data: bytes = field(repr=False)

    def _repr_png_(self) -> bytes:
        """Jupyter PNG display hook."""
        return self.data

    def _repr_html_(self) -> str:
        """Jupyter HTML display hook (inline base64 image)."""
        b64 = base64.b64encode(self.data).decode("ascii")
        return f'<img src="data:image/png;base64,{b64}" />'

    def save(self, path: str | Path) -> None:
        """Save the screenshot to a file.

        Parameters
        ----------
        path
            File path to save to.
        """
        Path(path).write_bytes(self.data)


# =============================================================================
# Selections
# =============================================================================


class Selections(MutableMapping[str, tuple[int, ...]]):
    """Accessor for per-geometry selections."""

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, geometry: str) -> tuple[int, ...]:
        response = self._api.get_selection(geometry)
        return tuple(response.get("selection", []))

    def __setitem__(self, geometry: str, indices: Iterable[int]) -> None:
        self._api.update_selection(geometry, list(indices))

    def __delitem__(self, geometry: str) -> None:
        self._api.update_selection(geometry, [])

    def __iter__(self) -> Iterator[str]:
        return iter(self._api.list_geometries().keys())

    def __len__(self) -> int:
        return len(self._api.list_geometries())


# =============================================================================
# SelectionGroups
# =============================================================================


class SelectionGroups(MutableMapping[str, dict[str, list[int]]]):
    """Accessor for named selection groups."""

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, group_name: str) -> dict[str, list[int]]:
        response = self._api.get_selection_group(group_name)
        return response["group"]

    def __setitem__(self, group_name: str, selections: dict[str, list[int]]) -> None:
        self._api.set_selection_group(group_name, selections)

    def __delitem__(self, group_name: str) -> None:
        self._api.delete_selection_group(group_name)

    def __iter__(self) -> Iterator[str]:
        return iter(self._api.list_selection_groups().keys())

    def __len__(self) -> int:
        return len(self._api.list_selection_groups())


# =============================================================================
# Bookmarks
# =============================================================================


class Bookmarks(MutableMapping[int, str]):
    """Accessor for frame bookmarks."""

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, index: int) -> str:
        response = self._api.get_bookmark(index)
        return response["label"]

    def __setitem__(self, index: int, label: str) -> None:
        self._api.set_bookmark(index, label)

    def __delitem__(self, index: int) -> None:
        self._api.delete_bookmark(index)

    def __iter__(self) -> Iterator[int]:
        return iter(int(k) for k in self._api.get_all_bookmarks())

    def __len__(self) -> int:
        return len(self._api.get_all_bookmarks())


# =============================================================================
# Geometries
# =============================================================================


class Geometries(MutableMapping[str, BaseGeometry]):
    """Accessor for custom geometries via Pydantic models."""

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, key: str) -> BaseGeometry:
        geom = self._api.get_geometry(key)
        if geom is None:
            raise KeyError(key)
        model_cls = geometry_models.get(geom["type"])
        if model_cls is None:
            raise ValueError(f"Unknown geometry type: {geom['type']}")
        return model_cls(**geom["data"])

    def __setitem__(self, key: str, value: BaseGeometry) -> None:
        geometry_type = type(value).__name__
        if geometry_type not in geometry_models:
            raise ValueError(f"Unknown geometry type: {geometry_type}")
        self._api.set_geometry(key, geometry_type, value.model_dump())

    def __delitem__(self, key: str) -> None:
        self._api.delete_geometry(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self._api.list_geometries().keys())

    def __len__(self) -> int:
        return len(self._api.list_geometries())


# =============================================================================
# Figures
# =============================================================================


class Figures(MutableMapping[str, "go.Figure"]):
    """Accessor for Plotly figures stored on the server.

    Uses ``FigureData`` Pydantic model for serialization/deserialization.
    """

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, key: str) -> go.Figure:
        from zndraw.schemas import FigureData

        figure = self._api.get_figure(key)
        if figure is None:
            raise KeyError(key)
        return FigureData(**figure).to_figure()

    def __setitem__(self, key: str, value: go.Figure) -> None:
        from zndraw.schemas import FigureData

        self._api.set_figure(key, FigureData.from_figure(value).model_dump())

    def __delitem__(self, key: str) -> None:
        self._api.delete_figure(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self._api.list_figures())

    def __len__(self) -> int:
        return len(self._api.list_figures())


# =============================================================================
# Presets
# =============================================================================


class Presets(MutableMapping[str, Preset]):
    """Accessor for visual presets stored per-room.

    Provides dict-like access to presets plus ``apply``, ``load``, and ``export``.

    Examples
    --------
    >>> vis.presets["matt"]
    Preset(name="matt", ...)

    >>> vis.presets["custom"] = Preset(
    ...     name="custom",
    ...     rules=[PresetRule(pattern="fog", config={"active": True})],
    ... )

    >>> result = vis.presets.apply("matt")
    >>> result.geometries_updated
    ["particles", "bonds", "fog", ...]

    >>> del vis.presets["custom"]
    """

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, name: str) -> Preset:
        data = self._api.get_preset(name)
        if data is None:
            raise KeyError(name)
        return Preset(**data)

    def __setitem__(self, name: str, value: Preset) -> None:
        self._api.update_preset(name, value.model_dump())

    def __delitem__(self, name: str) -> None:
        self._api.delete_preset(name)

    def __iter__(self) -> Iterator[str]:
        return iter(p["name"] for p in self._api.list_presets())

    def __len__(self) -> int:
        return len(self._api.list_presets())

    def apply(self, name: str) -> PresetApplyResult:
        """Apply a preset to all matching geometries in the room."""
        result = self._api.apply_preset(name)
        return PresetApplyResult(**result)

    def load(self, path: Path) -> Preset:
        """Load a preset from a JSON file into the room."""
        data = json.loads(path.read_text())
        result = self._api.create_preset(data)
        return Preset(**result)

    def export(self, name: str, path: Path) -> None:
        """Export a preset to a JSON file."""
        data = self._api.get_preset(name)
        if data is None:
            raise KeyError(name)
        data.pop("created_at", None)
        data.pop("updated_at", None)
        path.write_text(json.dumps(data, indent=2))


# =============================================================================
# RoomMetadata
# =============================================================================


class RoomMetadata(MutableMapping[str, str]):
    """Accessor for room metadata."""

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, key: str) -> str:
        info = self._api.get_room_info()
        metadata = info.get("metadata", {})
        if key not in metadata:
            raise KeyError(key)
        return metadata[key]

    def __setitem__(self, key: str, value: str) -> None:
        raise NotImplementedError("Room metadata updates not yet implemented")

    def __delitem__(self, key: str) -> None:
        raise NotImplementedError("Room metadata deletion not yet implemented")

    def __iter__(self) -> Iterator[str]:
        info = self._api.get_room_info()
        return iter(info.get("metadata", {}).keys())

    def __len__(self) -> int:
        info = self._api.get_room_info()
        return len(info.get("metadata", {}))


# =============================================================================
# Session / Sessions
# =============================================================================


@dataclass
class Session:
    """Proxy to a single frontend browser session.

    Validates on creation that the session exists (raises ``KeyError``
    otherwise). Each subsequent property access is a live HTTP call —
    no caching. Camera access resolves through ``active_camera``.
    """

    _api: APIManager = field(repr=False)
    sid: str

    def __post_init__(self) -> None:
        """Validate the session exists in active-cameras."""
        self._api.get_active_camera(self.sid)

    @property
    def camera(self) -> Camera:
        """Camera the session is currently viewing through."""
        key = self.active_camera
        geom = self._api.get_geometry(key)
        if geom is None:
            raise KeyError(f"Camera geometry '{key}' not found")
        return Camera(**geom["data"])

    @camera.setter
    def camera(self, value: Camera) -> None:
        key = self.active_camera
        self._api.set_geometry(key, "Camera", value.model_dump())

    @property
    def active_camera(self) -> str:
        """Geometry key of the camera this session views through."""
        return self._api.get_active_camera(self.sid)

    @active_camera.setter
    def active_camera(self, value: str) -> None:
        self._api.set_active_camera(self.sid, value)

    def screenshot(self, timeout: float = 30.0) -> ScreenshotImage:
        """Capture a screenshot from this frontend session.

        Parameters
        ----------
        timeout
            Maximum seconds to wait for capture completion.

        Returns
        -------
        ScreenshotImage
            The captured image with Jupyter display support.

        Raises
        ------
        TimeoutError
            If the screenshot is not completed within the timeout.
        """
        result = self._api.create_screenshot_capture(self.sid)
        screenshot_id = result["id"]
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            data = self._api.get_screenshot(screenshot_id)
            if data["status"] == "completed" and data["data"] is not None:
                return ScreenshotImage(data=base64.b64decode(data["data"]))
            time.sleep(0.2)
        raise TimeoutError(f"Screenshot not completed within {timeout}s")

    def __str__(self) -> str:
        return f"Session({self.sid!r})"


@dataclass(frozen=True)
class Sessions(Mapping[str, Session]):
    """Room-scoped mapping of active frontend sessions.

    Read-only: sessions are created by browsers, not Python.
    """

    _api: APIManager = field(repr=False)

    def __getitem__(self, sid: str) -> Session:
        return Session(_api=self._api, sid=sid)

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        return any(item.sid == key for item in self._api.list_sessions())

    def __iter__(self) -> Iterator[str]:
        return iter(item.sid for item in self._api.list_sessions())

    def __len__(self) -> int:
        return len(self._api.list_sessions())

    def __repr__(self) -> str:
        sids = list(self)
        return f"Sessions({sids!r})"

    def __str__(self) -> str:
        return f"Sessions(n={len(self)})"


# =============================================================================
# ChatMessages (new)
# =============================================================================


class ChatMessages(Sequence["MessageResponse"]):
    """Read-only sequence of chat messages with send capability.

    Fetches messages from the server on each access (no caching).
    Supports indexing, slicing, len, and iteration.
    """

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def _fetch_all(self) -> list[MessageResponse]:
        from zndraw.schemas import MessageResponse, MessagesResponse

        data = self._api.list_chat_messages()
        resp = MessagesResponse.model_validate(data)
        return list(reversed(resp.items))

    @overload
    def __getitem__(self, index: int) -> MessageResponse: ...
    @overload
    def __getitem__(self, index: slice) -> list[MessageResponse]: ...
    def __getitem__(
        self, index: int | slice
    ) -> MessageResponse | list[MessageResponse]:
        items = self._fetch_all()
        return items[index]

    def __len__(self) -> int:
        from zndraw.schemas import MessagesResponse

        data = self._api.list_chat_messages(limit=1)
        resp = MessagesResponse.model_validate(data)
        return resp.metadata.total_count

    def __iter__(self) -> Iterator[MessageResponse]:
        return iter(self._fetch_all())

    def send(self, message: str) -> None:
        """Send a chat message to the room.

        Parameters
        ----------
        message
            The message content (supports Markdown).
        """
        self._api.create_chat_message(message)

    def __repr__(self) -> str:
        return f"ChatMessages(n={len(self)})"


# =============================================================================
# Screenshots (new)
# =============================================================================


class Screenshots(Mapping[int, ScreenshotImage]):
    """Read-only mapping of screenshots by ID.

    ``list(vis.screenshots)`` returns screenshot IDs.
    ``vis.screenshots[id]`` returns a ``ScreenshotImage``.
    """

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, screenshot_id: int) -> ScreenshotImage:
        data = self._api.get_screenshot(screenshot_id)
        if data.get("data") is None:
            raise KeyError(f"Screenshot {screenshot_id} has no image data")
        return ScreenshotImage(data=base64.b64decode(data["data"]))

    def __iter__(self) -> Iterator[int]:
        return iter(item["id"] for item in self._api.list_screenshots())

    def __len__(self) -> int:
        return len(self._api.list_screenshots())

    def __repr__(self) -> str:
        return f"Screenshots(n={len(self)})"


# =============================================================================
# Extensions
# =============================================================================


class Extensions(Mapping[str, dict[str, Any]]):
    """Read-only mapping of available extensions by full_name.

    ``vis.extensions[name]`` returns the job metadata dict
    (name, category, schema_, ...).
    ``list(vis.extensions)`` returns extension full_names.
    """

    def __init__(self, api: APIManager) -> None:
        self._api = api

    def __getitem__(self, full_name: str) -> dict[str, Any]:
        try:
            return self._api.get_extension(full_name)
        except Exception:
            raise KeyError(full_name)

    def __iter__(self) -> Iterator[str]:
        data = self._api.list_extensions()
        return iter(item["full_name"] for item in data.get("items", []))

    def __len__(self) -> int:
        data = self._api.list_extensions()
        return data.get("total", len(data.get("items", [])))

    def __repr__(self) -> str:
        return f"Extensions(n={len(self)})"


# =============================================================================
# Tasks
# =============================================================================


@dataclass
class TaskHandle:
    """Handle for a submitted task with polling support."""

    id: str
    _api: APIManager = field(repr=False)

    @property
    def status(self) -> str:
        """Current task status."""
        return self._fetch().status

    def wait(self, *, timeout: float = 300, poll: float = 0.5) -> TaskHandle:
        """Poll until completed/failed. Returns self. Raises TimeoutError."""
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            data = self._fetch()
            if data.status in ("completed", "failed"):
                return self
            time.sleep(poll)
        raise TimeoutError(f"Task {self.id} did not complete within {timeout}s")

    def _fetch(self) -> Any:
        from zndraw_joblib.schemas import TaskResponse

        return TaskResponse.model_validate(self._api.get_task(self.id))

    def __str__(self) -> str:
        return self.id


class Tasks(Mapping[str, Any]):
    """Read-only mapping of tasks by ID with callable filter.

    ``list(vis.tasks)`` returns task IDs.
    ``vis.tasks[task_id]`` returns task details.
    ``vis.tasks(status='running')`` returns a filtered view.
    """

    def __init__(self, api: APIManager, status: str | None = None) -> None:
        self._api = api
        self._status = status

    def __call__(self, status: str) -> Tasks:
        """Return a filtered view of tasks."""
        return Tasks(self._api, status=status)

    def __getitem__(self, task_id: str) -> TaskHandle:
        try:
            self._api.get_task(task_id)
        except Exception:
            raise KeyError(task_id)
        return TaskHandle(id=task_id, _api=self._api)

    def __iter__(self) -> Iterator[str]:
        data = self._api.list_tasks(status=self._status)
        return iter(str(item["id"]) for item in data.get("items", []))

    def __len__(self) -> int:
        data = self._api.list_tasks(status=self._status)
        return data.get("total", len(data.get("items", [])))

    def __repr__(self) -> str:
        suffix = f", status={self._status!r}" if self._status else ""
        return f"Tasks(n={len(self)}{suffix})"
