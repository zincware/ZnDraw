"""ZnDraw Python client for interacting with the server.

This module provides a synchronous client that mirrors the reference
ZnDraw implementation. It includes:
- ZnDraw: Main client class implementing MutableSequence for frame operations
- APIManager: HTTP REST API wrapper
- SocketManager: WebSocket real-time synchronization
- ZnDrawLock: Distributed locking context manager
- Helper classes: Selections, Bookmarks, Geometries, Figures, RoomMetadata
"""

from __future__ import annotations

import base64
import contextlib
import json
import logging
import threading
import uuid
from collections.abc import (
    Generator,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    MutableSequence,
)
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast, overload

if TYPE_CHECKING:
    from zndraw.providers.frame_source import FrameSource
    from zndraw.tqdm import ZnDrawTqdm

import ase
import httpx
import msgpack
import socketio
import typing_extensions
from asebytes import decode, encode
from pydantic import SecretStr
from zndraw_joblib.client import ClaimedTask, Extension as JoblibExtension, JobManager
from zndraw_joblib.exceptions import ProviderTimeoutError
from zndraw_socketio import SyncClientWrapper, wrap

from zndraw.enrichment import add_colors, add_radii
from zndraw.exceptions import PROBLEM_TYPES, ProblemDetail
from zndraw.geometries import geometries as geometry_models
from zndraw.geometries.base import BaseGeometry
from zndraw.geometries.camera import Camera
from zndraw.settings import RoomConfig

if TYPE_CHECKING:
    import plotly.graph_objects as go

    from zndraw.extensions.abc import Extension


log = logging.getLogger(__name__)


# =============================================================================
# Frame Serialization Helpers
# =============================================================================


def atoms_to_json_dict(
    atoms: ase.Atoms, connectivity_threshold: int = 1000
) -> dict[str, Any]:
    """Convert ase.Atoms to a JSON-compatible dictionary.

    Ensures colors, radii, and connectivity are present, then uses asebytes.encode()
    to get dict[bytes, bytes], then converts:
    - Byte keys to base64-encoded strings with "b64:" prefix
    - Byte values to base64-encoded strings

    Parameters
    ----------
    atoms
        The ASE Atoms object to convert.
    connectivity_threshold
        Maximum number of atoms for automatic connectivity calculation.
        Connectivity is only computed if the number of atoms is below this
        threshold and connectivity is not already present. Default: 1000.
    """
    add_colors(atoms)
    add_radii(atoms)

    if len(atoms) < connectivity_threshold and "connectivity" not in atoms.info:
        from zndraw.connectivity import add_connectivity

        add_connectivity(atoms)

    encoded = encode(atoms)
    result: dict[str, Any] = {}
    for key, value in encoded.items():
        # Encode key as base64 string
        key_str = "b64:" + base64.b64encode(key).decode("ascii")
        # Encode value as base64 string
        value_str = base64.b64encode(value).decode("ascii")
        result[key_str] = value_str
    return result


def json_dict_to_atoms(data: dict[str, Any]) -> ase.Atoms:
    """Convert a JSON dictionary back to ase.Atoms.

    Reverses atoms_to_json_dict(): decodes base64 keys/values back to bytes,
    then uses asebytes.decode().
    """
    encoded: dict[bytes, bytes] = {}
    for key_str, value_str in data.items():
        # Decode key from base64 (strip "b64:" prefix if present)
        if key_str.startswith("b64:"):
            key = base64.b64decode(key_str[4:])
        else:
            # Legacy format or string key - convert to bytes
            key = key_str.encode("utf-8") if isinstance(key_str, str) else key_str
        # Decode value from base64
        value = base64.b64decode(value_str) if isinstance(value_str, str) else value_str
        encoded[key] = value
    return decode(encoded)


def raw_frame_to_atoms(frame: dict[bytes, bytes]) -> ase.Atoms:
    """Convert a raw msgpack frame (dict[bytes, bytes]) to ase.Atoms.

    This is used when frames come directly from msgpack (server response).
    """
    return decode(frame)


# Target byte size per upload chunk. Frames are accumulated until adding
# the next frame would exceed this limit, then the chunk is flushed.
_TARGET_CHUNK_BYTES = 2_000_000  # 2 MB

# Hard cap on frames per chunk — must not exceed the server's
# FrameCreateRequest.max_length (1000).
_MAX_CHUNK_FRAMES = 1000


def _estimate_frame_size(frame: dict[str, Any]) -> int:
    """Estimate serialized size of a JSON-encoded frame dict.

    Uses the sum of base64 value lengths as a cheap proxy — no
    extra serialization needed since we already have the strings.
    Intentionally ignores key lengths (~20-30% undercount) since this
    is only used for chunking heuristics, not exact measurement.
    """
    return sum(len(v) for v in frame.values() if isinstance(v, str))


# =============================================================================
# Exceptions
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


class ZnDrawError(Exception):
    """Base exception for ZnDraw client errors."""


class NotConnectedError(ZnDrawError):
    """Raised when operation requires connection but client is disconnected."""


class RoomLockedError(ZnDrawError):
    """Raised when the room is locked (admin lock or edit lock by another user)."""


# =============================================================================
# ZnDrawLock - Edit Lock Context Manager
# =============================================================================


@dataclass
class ZnDrawLock:
    """Context manager for room-level edit lock via REST API.

    Acquires the edit lock on enter, releases on exit.
    Automatically refreshes the lock every EDIT_LOCK_REFRESH seconds.

    Parameters
    ----------
    api : APIManager
        The API manager instance.
    msg : str | None
        Optional message describing the lock purpose.
    """

    api: APIManager
    msg: str | None = None
    _stop: threading.Event = field(default_factory=threading.Event, init=False)
    _refresh_thread: threading.Thread | None = field(default=None, init=False)
    _lock_token: str | None = field(default=None, init=False)

    def __enter__(self) -> ZnDrawLock:
        """Acquire the edit lock."""
        result = self.api.edit_lock_acquire(self.msg)
        self._lock_token = result["lock_token"]
        self._stop.clear()
        self._refresh_thread = threading.Thread(target=self._refresh, daemon=True)
        self._refresh_thread.start()
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> bool:
        """Release the edit lock."""
        self._stop.set()
        if self._refresh_thread is not None:
            self._refresh_thread.join(timeout=2.0)
        try:
            self.api.edit_lock_release(self._lock_token)
        except Exception as e:
            log.warning(f"Failed to release edit lock: {e}")
        self._lock_token = None
        return False

    def _refresh(self) -> None:
        """Periodically refresh the edit lock."""
        from zndraw.schemas import EDIT_LOCK_REFRESH

        while not self._stop.wait(EDIT_LOCK_REFRESH):
            if self._lock_token is None:
                break
            try:
                self.api.edit_lock_refresh(self._lock_token, self.msg)
            except Exception as e:
                log.warning(f"Failed to refresh edit lock: {e}")
                break

    @property
    def lock_token(self) -> str | None:
        """The current lock token, or None if not held."""
        return self._lock_token


# =============================================================================
# APIManager - REST API Wrapper
# =============================================================================


@dataclass
class APIManager:
    """Synchronous HTTP client for ZnDraw REST API.

    Handles authentication, request building, and error handling.
    """

    url: str
    room_id: str
    token: str | None = None
    session_id: str | None = None
    user_id: int | None = None
    username: str | None = None

    http: httpx.Client = field(init=False)

    def __post_init__(self) -> None:
        """Initialize HTTP client."""
        self.http = httpx.Client(base_url=self.url, timeout=30.0)

    @property
    def base_url(self) -> str:
        """Base URL (satisfies ApiManager protocol)."""
        return self.url

    def get_headers(self) -> dict[str, str]:
        """Build request headers (satisfies ApiManager protocol)."""
        return self._headers()

    def close(self) -> None:
        """Close the HTTP client."""
        self.http.close()

    def _headers(self) -> dict[str, str]:
        """Build request headers."""
        headers: dict[str, str] = {}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        if self.session_id:
            headers["X-Session-ID"] = self.session_id
        return headers

    def raise_for_status(self, response: httpx.Response) -> None:
        """Raise appropriate exception for error responses.

        Uses the RFC 9457 ``type`` URI to look up the ``ProblemType`` in
        the registry and delegates to its ``raise_for_client`` method.
        Falls back to ``ZnDrawError`` for non-RFC 9457 responses.
        """
        if response.status_code < 400:
            return
        content_type = response.headers.get("content-type", "")
        if "application/problem+json" in content_type:
            error_data = response.json()
            type_uri = error_data.get("type", "")
            if type_uri:
                problem_id = type_uri.rsplit("/", 1)[-1]
                problem_cls = PROBLEM_TYPES.get(problem_id)
                if problem_cls:
                    problem = ProblemDetail(**error_data)
                    problem_cls.raise_for_client(problem)
        try:
            error_data = response.json()
        except (ValueError, json.JSONDecodeError):
            response.raise_for_status()
            return

        # Fallback: use status code for unregistered or non-RFC 9457 responses
        detail = str(error_data.get("detail", response.text))
        status_code = response.status_code
        if status_code == 404:
            raise KeyError(detail)
        if status_code in {401, 403}:
            raise PermissionError(detail)
        if status_code in {409, 422}:
            raise ValueError(detail)
        if status_code == 423:
            raise RoomLockedError(detail)
        raise ZnDrawError(f"{status_code}: {detail}")

    # -------------------------------------------------------------------------
    # Authentication
    # -------------------------------------------------------------------------

    def create_guest_session(self) -> dict[str, Any]:
        """Create a guest session and get JWT token."""
        response = self.http.post("/v1/auth/guest")
        self.raise_for_status(response)
        data = response.json()
        self.token = data["access_token"]
        return data

    def login(self, email: str, password: str) -> dict[str, Any]:
        """Login with email and password."""
        response = self.http.post(
            "/v1/auth/jwt/login",
            data={"username": email, "password": password},
        )
        self.raise_for_status(response)
        data = response.json()
        self.token = data["access_token"]
        return data

    # -------------------------------------------------------------------------
    # Room Operations
    # -------------------------------------------------------------------------

    def create_room(
        self,
        description: str | None = None,
        copy_from: str | None = None,
    ) -> dict[str, Any]:
        """Create a new room.

        Parameters
        ----------
        description : str, optional
            Room description.
        copy_from : str, optional
            Room ID to copy frames from, or an @-prefixed preset
            (``@empty`` for one empty frame, ``@none`` for zero frames).
        """
        payload: dict[str, Any] = {"room_id": self.room_id}
        if description is not None:
            payload["description"] = description
        if copy_from is not None:
            payload["copy_from"] = copy_from

        response = self.http.post(
            "/v1/rooms",
            json=payload,
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_room_info(self) -> dict[str, Any]:
        """Get room information."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_room(self, updates: dict[str, Any]) -> dict[str, Any]:
        """Update room metadata via PATCH.

        Parameters
        ----------
        updates
            Fields to update (e.g. ``{"frame_count": 100}``).
        """
        response = self.http.patch(
            f"/v1/rooms/{self.room_id}",
            json=updates,
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Frame Operations
    # -------------------------------------------------------------------------

    def get_frame(self, index: int) -> dict[bytes, bytes]:
        """Get a single frame by index via the individual frame endpoint.

        Uses ``GET /v1/rooms/{room_id}/frames/{index}`` which supports
        provider dispatch. Retries on 504 Retry-After until the provider
        delivers the frame.

        Parameters
        ----------
        index
            Frame index (0-based).
        """
        import time

        while True:
            response = self.http.get(
                f"/v1/rooms/{self.room_id}/frames/{index}",
                headers=self._headers(),
            )
            try:
                self.raise_for_status(response)
            except ProviderTimeoutError:
                retry_after = int(response.headers.get("Retry-After", "2"))
                time.sleep(retry_after)
                continue
            frames: list[dict[bytes, bytes]] = msgpack.unpackb(
                response.content, strict_map_key=False
            )
            return frames[0]

    def get_frames(
        self,
        indices: list[int] | None = None,
        start: int | None = None,
        stop: int | None = None,
        keys: list[str] | None = None,
    ) -> list[dict[bytes, bytes]]:
        """Get frames by indices or range.

        Retries on 504 Retry-After for provider-backed rooms.
        Returns raw msgpack-decoded frames: list[dict[bytes, bytes]].
        Each frame has bytes keys (e.g., b"arrays.positions") and bytes values.
        """
        import time

        params: dict[str, Any] = {}
        if indices is not None:
            params["indices"] = ",".join(str(i) for i in indices)
        else:
            if start is not None:
                params["start"] = start
            if stop is not None:
                params["stop"] = stop
        if keys is not None:
            params["keys"] = ",".join(keys)

        while True:
            response = self.http.get(
                f"/v1/rooms/{self.room_id}/frames",
                params=params,
                headers=self._headers(),
            )
            try:
                self.raise_for_status(response)
            except ProviderTimeoutError:
                retry_after = int(response.headers.get("Retry-After", "2"))
                time.sleep(retry_after)
                continue
            return msgpack.unpackb(response.content, strict_map_key=False)

    def append_frames(self, frames: list[dict[str, Any]]) -> dict[str, Any]:
        """Append frames to the room."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/frames",
            json={"frames": frames},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_frame(self, index: int, data: dict[str, Any]) -> dict[str, Any]:
        """Update a single frame."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/frames/{index}",
            json={"data": data},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_frame(self, index: int) -> None:
        """Delete a single frame."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/frames/{index}",
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Step Operations
    # -------------------------------------------------------------------------

    def get_step(self) -> dict[str, Any]:
        """Get current step."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/step",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_step(self, step: int) -> dict[str, Any]:
        """Update current step."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/step",
            json={"step": step},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Frame Selection Operations
    # -------------------------------------------------------------------------

    def get_frame_selection(self) -> list[int] | None:
        """Get selected frame indices."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/frame-selection",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["frame_selection"]

    def update_frame_selection(self, indices: list[int]) -> dict[str, Any]:
        """Update selected frame indices."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/frame-selection",
            json={"indices": indices},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Geometry Operations
    # -------------------------------------------------------------------------

    def list_geometries(self) -> dict[str, Any]:
        """List all geometries."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/geometries",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_geometry(self, key: str) -> dict[str, Any] | None:
        """Get a specific geometry."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/geometries/{key}",
            headers=self._headers(),
        )
        if response.status_code == 404:
            return None
        self.raise_for_status(response)
        return response.json()["geometry"]

    def set_geometry(
        self, key: str, geometry_type: str, data: dict[str, Any]
    ) -> dict[str, Any]:
        """Create or update a geometry."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/geometries/{key}",
            json={"type": geometry_type, "data": data},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_geometry(self, key: str) -> None:
        """Delete a geometry."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/geometries/{key}",
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Selection Operations
    # -------------------------------------------------------------------------

    def get_selection(self, geometry: str) -> dict[str, Any]:
        """Get selection for a specific geometry."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/geometries/{geometry}/selection",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_selection(self, geometry: str, indices: list[int]) -> dict[str, Any]:
        """Update selection for a geometry."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/geometries/{geometry}/selection",
            json={"indices": indices},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def list_selection_groups(self) -> dict[str, dict[str, list[int]]]:
        """List all selection groups."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/selection-groups",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_selection_group(self, group_name: str) -> dict[str, Any]:
        """Get a selection group."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/selection-groups/{group_name}",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def set_selection_group(
        self, group_name: str, selections: dict[str, list[int]]
    ) -> dict[str, Any]:
        """Create or update a selection group."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/selection-groups/{group_name}",
            json={"selections": selections},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_selection_group(self, group_name: str) -> None:
        """Delete a selection group."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/selection-groups/{group_name}",
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Bookmark Operations
    # -------------------------------------------------------------------------

    def get_all_bookmarks(self) -> dict[str, str]:
        """Get all bookmarks."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/bookmarks",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_bookmark(self, index: int) -> dict[str, Any]:
        """Get a specific bookmark."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def set_bookmark(self, index: int, label: str) -> dict[str, Any]:
        """Set a bookmark."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            json={"label": label},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_bookmark(self, index: int) -> None:
        """Delete a bookmark."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Figure Operations
    # -------------------------------------------------------------------------

    def list_figures(self) -> list[str]:
        """List all figure keys."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/figures",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_figure(self, key: str) -> dict[str, Any] | None:
        """Get a specific figure."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/figures/{key}",
            headers=self._headers(),
        )
        if response.status_code == 404:
            return None
        self.raise_for_status(response)
        return response.json()["figure"]

    def set_figure(self, key: str, figure: dict[str, Any]) -> dict[str, Any]:
        """Create or update a figure."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/figures/{key}",
            json={"figure": figure},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_figure(self, key: str) -> None:
        """Delete a figure."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/figures/{key}",
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Edit Lock Operations
    # -------------------------------------------------------------------------

    def edit_lock_acquire(self, msg: str | None = None) -> dict[str, Any]:
        """Acquire the room edit lock."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/edit-lock",
            json={"msg": msg},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def edit_lock_refresh(
        self, lock_token: str, msg: str | None = None
    ) -> dict[str, Any]:
        """Refresh an existing edit lock using its token."""
        headers = {**self._headers(), "Lock-Token": lock_token}
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/edit-lock",
            json={"msg": msg},
            headers=headers,
        )
        self.raise_for_status(response)
        return response.json()

    def edit_lock_release(self, lock_token: str | None = None) -> None:
        """Release the room edit lock."""
        headers = self._headers()
        if lock_token:
            headers["Lock-Token"] = lock_token
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/edit-lock",
            headers=headers,
        )
        self.raise_for_status(response)

    def edit_lock_status(self) -> dict[str, Any]:
        """Get the room edit lock status."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/edit-lock",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Session Operations
    # -------------------------------------------------------------------------

    def list_sessions(self) -> list[str]:
        """List the current user's active frontend sessions."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/sessions",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_active_camera(self, sid: str) -> str:
        """Get active camera key for a session."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/sessions/{sid}/active-camera",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        active = response.json()["active_camera"]
        if active is None:
            raise KeyError(f"No active camera for session {sid}")
        return active

    def set_active_camera(self, sid: str, camera_key: str) -> None:
        """Set active camera key for a session."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/sessions/{sid}/active-camera",
            json={"active_camera": camera_key},
            headers=self._headers(),
        )
        self.raise_for_status(response)

    def get_session_settings(self, sid: str) -> dict[str, Any]:
        """Get session settings data."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/sessions/{sid}/settings",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()["data"]

    def set_session_settings(self, sid: str, settings: dict[str, Any]) -> None:
        """Update session settings."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/sessions/{sid}/settings",
            json=settings,
            headers=self._headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Chat Operations
    # -------------------------------------------------------------------------

    def create_chat_message(self, content: str) -> dict[str, Any]:
        """Send a chat message to the room."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/chat/messages",
            json={"content": content},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Screenshot Operations
    # -------------------------------------------------------------------------

    def create_screenshot_capture(self, session_id: str) -> dict[str, Any]:
        """Request a programmatic screenshot capture from a frontend session."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/screenshots",
            json={"session_id": session_id},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_screenshot(self, screenshot_id: int) -> dict[str, Any]:
        """Get a screenshot by ID."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/screenshots/{screenshot_id}",
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Progress Operations
    # -------------------------------------------------------------------------

    def progress_start(
        self, progress_id: str, description: str, unit: str = "it"
    ) -> dict[str, Any]:
        """Create a new progress tracker in the room."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/progress",
            json={
                "progress_id": progress_id,
                "description": description,
                "unit": unit,
            },
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def progress_update(
        self,
        progress_id: str,
        description: str | None = None,
        n: int | None = None,
        total: int | None = None,
        elapsed: float | None = None,
        unit: str | None = None,
    ) -> dict[str, Any]:
        """Update an existing progress tracker."""
        payload: dict[str, Any] = {}
        if description is not None:
            payload["description"] = description
        if n is not None:
            payload["n"] = n
        if total is not None:
            payload["total"] = total
        if elapsed is not None:
            payload["elapsed"] = elapsed
        if unit is not None:
            payload["unit"] = unit
        response = self.http.patch(
            f"/v1/rooms/{self.room_id}/progress/{progress_id}",
            json=payload,
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def progress_complete(self, progress_id: str) -> None:
        """Complete and remove a progress tracker."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/progress/{progress_id}",
            headers=self._headers(),
        )
        self.raise_for_status(response)


# =============================================================================
# SocketManager - WebSocket Real-time Sync
# =============================================================================


@dataclass
class SocketManager:
    """Manages Socket.IO connection and real-time synchronization."""

    zndraw: ZnDraw
    _tsio: SyncClientWrapper = field(init=False)
    _connected: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        """Initialize wrapped Socket.IO client and register handlers."""
        self._tsio = wrap(socketio.Client())
        self._register_handlers()

    def _register_handlers(self) -> None:
        """Register Socket.IO event handlers using typed models."""
        from zndraw.socket_events import FramesInvalidate

        self._tsio.on("connect", self._on_connect)
        self._tsio.on("disconnect", self._on_disconnect)
        self._tsio.on(FramesInvalidate, self._on_frames_invalidate)

    @property
    def connected(self) -> bool:
        """Check if connected."""
        return self._connected and self._tsio.connected

    def connect(self) -> None:
        """Connect to the server and join the room."""
        from zndraw.socket_events import RoomJoin, RoomJoinResponse

        if self._tsio.connected:
            log.debug("Already connected")
            return

        # Connect with JWT token
        self._tsio.connect(
            self.zndraw.url,
            auth={"token": self.zndraw.api.token},
            wait=True,
        )

        # Join room and get session info
        if self.zndraw.room is None:
            raise NotConnectedError("Cannot join: no room set")
        join_request = RoomJoin(room_id=self.zndraw.room, client_type="pyclient")
        raw_response = self._tsio.call(join_request)
        response: dict[str, Any] = raw_response if raw_response else {}

        if "type" in response:
            # RFC 9457 error response - room might not exist
            if response.get("status") == 404:
                self.zndraw.api.create_room(copy_from="@none")
                # Retry join
                raw_response = self._tsio.call(join_request)
                response = raw_response if raw_response else {}

        if "type" in response:
            error_msg = (
                response.get("detail") or response.get("title") or "Unknown error"
            )
            raise ZnDrawError(f"Failed to join room: {error_msg}")

        # Validate response
        join_response = RoomJoinResponse.model_validate(response)
        self.zndraw.api.session_id = join_response.session_id

        # Seed the frame count cache from join response
        self.zndraw._cached_length = join_response.frame_count

        self._connected = True
        log.debug(f"Connected to room {self.zndraw.room}")

    def disconnect(self) -> None:
        """Disconnect from the server."""
        if self._tsio.connected:
            self._tsio.disconnect()
        self._connected = False
        log.debug("Disconnected")

    def wait(self) -> None:
        """Block until disconnected."""
        self._tsio.wait()

    def _on_connect(self) -> None:
        """Handle connection event. Re-register providers if mounted."""
        log.debug("Socket connected")
        if self.zndraw._mount is not None:
            from zndraw.providers.frame_source import FrameSourceLength, FrameSourceRead

            source = self.zndraw._mount
            room = self.zndraw.room
            if room is None:
                return
            mount_name = self.zndraw._mount_name
            if mount_name is None:
                return
            self.zndraw.register_provider(
                FrameSourceRead, name=mount_name, handler=source, room=room
            )
            self.zndraw.register_provider(
                FrameSourceLength, name=mount_name, handler=source, room=room
            )
            self.zndraw.api.update_room({"frame_count": len(source)})

    def _on_disconnect(self, reason: str = "") -> None:
        """Handle disconnection event."""
        log.debug("Socket disconnected: %s", reason)
        self._connected = False

    def _on_frames_invalidate(self, data: Any) -> None:
        """Handle FramesInvalidate broadcast — update cached length."""
        from zndraw.socket_events import FramesInvalidate

        event = FramesInvalidate.model_validate(data)
        if event.room_id != self.zndraw.room:
            return
        if event.count is not None:
            self.zndraw._cached_length = event.count
        else:
            self.zndraw._cached_length = None


# =============================================================================
# Helper Classes - Dict-like Accessors
# =============================================================================


class Selections(MutableMapping[str, tuple[int, ...]]):
    """Accessor for per-geometry selections."""

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, geometry: str) -> tuple[int, ...]:
        response = self._zndraw.api.get_selection(geometry)
        return tuple(response.get("selection", []))

    def __setitem__(self, geometry: str, indices: Iterable[int]) -> None:
        indices_list = list(indices)
        self._zndraw.api.update_selection(geometry, indices_list)

    def __delitem__(self, geometry: str) -> None:
        self._zndraw.api.update_selection(geometry, [])

    def __iter__(self):
        geometries = self._zndraw.api.list_geometries()
        return iter(geometries.keys())

    def __len__(self) -> int:
        geometries = self._zndraw.api.list_geometries()
        return len(geometries)


class SelectionGroups(MutableMapping[str, dict[str, list[int]]]):
    """Accessor for named selection groups."""

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, group_name: str) -> dict[str, list[int]]:
        response = self._zndraw.api.get_selection_group(group_name)
        return response["group"]

    def __setitem__(self, group_name: str, selections: dict[str, list[int]]) -> None:
        self._zndraw.api.set_selection_group(group_name, selections)

    def __delitem__(self, group_name: str) -> None:
        self._zndraw.api.delete_selection_group(group_name)

    def __iter__(self):
        groups = self._zndraw.api.list_selection_groups()
        return iter(groups.keys())

    def __len__(self) -> int:
        groups = self._zndraw.api.list_selection_groups()
        return len(groups)


class Bookmarks(MutableMapping[int, str]):
    """Accessor for frame bookmarks."""

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, index: int) -> str:
        response = self._zndraw.api.get_bookmark(index)
        return response["label"]

    def __setitem__(self, index: int, label: str) -> None:
        self._zndraw.api.set_bookmark(index, label)

    def __delitem__(self, index: int) -> None:
        self._zndraw.api.delete_bookmark(index)

    def __iter__(self):
        return iter(int(k) for k in self._zndraw.api.get_all_bookmarks())

    def __len__(self) -> int:
        return len(self._zndraw.api.get_all_bookmarks())


class Geometries(MutableMapping[str, BaseGeometry]):
    """Accessor for custom geometries via Pydantic models."""

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, key: str) -> BaseGeometry:
        geom = self._zndraw.api.get_geometry(key)
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
        self._zndraw.api.set_geometry(key, geometry_type, value.model_dump())

    def __delitem__(self, key: str) -> None:
        self._zndraw.api.delete_geometry(key)

    def __iter__(self):
        return iter(self._zndraw.api.list_geometries().keys())

    def __len__(self) -> int:
        return len(self._zndraw.api.list_geometries())


class Figures(MutableMapping[str, "go.Figure"]):
    """Accessor for Plotly figures stored on the server.

    Uses ``FigureData`` Pydantic model for serialization/deserialization.
    """

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, key: str) -> go.Figure:
        from zndraw.schemas import FigureData

        figure = self._zndraw.api.get_figure(key)
        if figure is None:
            raise KeyError(key)
        return FigureData(**figure).to_figure()

    def __setitem__(self, key: str, value: go.Figure) -> None:
        from zndraw.schemas import FigureData

        self._zndraw.api.set_figure(key, FigureData.from_figure(value).model_dump())

    def __delitem__(self, key: str) -> None:
        self._zndraw.api.delete_figure(key)

    def __iter__(self):
        return iter(self._zndraw.api.list_figures())

    def __len__(self) -> int:
        return len(self._zndraw.api.list_figures())


class RoomMetadata(MutableMapping[str, str]):
    """Accessor for room metadata."""

    def __init__(self, zndraw: ZnDraw) -> None:
        self._zndraw = zndraw

    def __getitem__(self, key: str) -> str:
        info = self._zndraw.api.get_room_info()
        metadata = info.get("metadata", {})
        if key not in metadata:
            raise KeyError(key)
        return metadata[key]

    def __setitem__(self, key: str, value: str) -> None:
        # Room metadata updates via PATCH
        raise NotImplementedError("Room metadata updates not yet implemented")

    def __delitem__(self, key: str) -> None:
        raise NotImplementedError("Room metadata deletion not yet implemented")

    def __iter__(self):
        info = self._zndraw.api.get_room_info()
        return iter(info.get("metadata", {}).keys())

    def __len__(self) -> int:
        info = self._zndraw.api.get_room_info()
        return len(info.get("metadata", {}))


@dataclass
class Session:
    """Proxy to a single frontend browser session.

    Validates on creation that the session exists and belongs to the
    current user (raises ``KeyError`` otherwise). Each subsequent
    property access is a live HTTP call — no caching.
    Camera access resolves through ``active_camera``.
    """

    _api: APIManager = field(repr=False)
    sid: str

    def __post_init__(self) -> None:
        """Validate the session exists and belongs to the current user.

        Uses the settings endpoint as a lightweight existence check.
        ``VerifiedSessionDep`` returns 404 for non-existent or non-owned
        sessions, which ``raise_for_status`` maps to ``KeyError``.
        """
        self._api.get_session_settings(self.sid)

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

    @property
    def settings(self) -> RoomConfig:
        """Session rendering settings (frozen -- use model_copy to modify)."""
        data = self._api.get_session_settings(self.sid)
        return RoomConfig(**data)

    @settings.setter
    def settings(self, value: RoomConfig) -> None:
        self._api.set_session_settings(self.sid, value.model_dump())

    def screenshot(self, timeout: float = 30.0) -> ScreenshotImage:
        """Capture a screenshot from this frontend session.

        Requests the frontend to render and upload a screenshot,
        then polls until the image is ready.

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
        import time

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
    """User-scoped mapping of active frontend sessions.

    Read-only: sessions are created by browsers, not Python.
    """

    _api: APIManager = field(repr=False)

    def __getitem__(self, sid: str) -> Session:
        return Session(_api=self._api, sid=sid)

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        return key in self._api.list_sessions()

    def __iter__(self) -> Iterator[str]:
        return iter(self._api.list_sessions())

    def __len__(self) -> int:
        return len(self._api.list_sessions())

    def __repr__(self) -> str:
        sids = list(self)
        return f"Sessions({sids!r})"

    def __str__(self) -> str:
        return f"Sessions(n={len(self)})"


# =============================================================================
# ZnDraw - Main Client Class
# =============================================================================


@dataclass
class ZnDraw(MutableSequence[ase.Atoms]):
    """A synchronous client for interacting with the ZnDraw server.

    Implements MutableSequence for frame operations, providing list-like access
    to trajectory frames as ase.Atoms objects.

    Parameters
    ----------
    url : str
        URL of the ZnDraw server.
    room : str | None
        Room ID to connect to. If None, generates a random UUID.
    user : str | None
        User email for authentication. If None, creates a guest session.
    password : SecretStr | str | None
        Password for login. Accepts ``str`` (auto-wrapped to ``SecretStr``)
        or ``SecretStr``. If None, inferred from ``Settings.guest_password``.
    auto_connect : bool
        If True, connect immediately on creation.

    Examples
    --------
    >>> import ase
    >>> vis = ZnDraw(url="http://localhost:8000")
    >>> len(vis)  # Number of frames
    >>> vis.append(ase.Atoms("H2O", positions=[[0,0,0], [1,0,0], [0,1,0]]))
    >>> atoms = vis[0]  # Get first frame as ase.Atoms
    >>> vis[0] = ase.Atoms("CO2")  # Update frame
    >>> del vis[0]  # Delete frame
    """

    url: str
    room: str | None = None
    user: str | None = None
    password: SecretStr | str | None = None
    auto_connect: bool = True
    auto_pickup: bool = True
    polling_interval: float = 5.0
    heartbeat_interval: float = 30.0

    # Internal state
    api: APIManager = field(init=False)
    socket: SocketManager = field(init=False)
    _jobs: JobManager = field(init=False)
    _cached_length: int | None = field(default=None, init=False, repr=False)
    _mount: FrameSource | None = field(default=None, init=False, repr=False)
    _mount_name: str | None = field(default=None, init=False, repr=False)

    # Accessors (lazy initialized)
    _selections: Selections | None = field(default=None, init=False)
    _selection_groups: SelectionGroups | None = field(default=None, init=False)
    _bookmarks: Bookmarks | None = field(default=None, init=False)
    _geometries: Geometries | None = field(default=None, init=False)
    _figures: Figures | None = field(default=None, init=False)
    _metadata: RoomMetadata | None = field(default=None, init=False)

    def __post_init__(self) -> None:
        """Initialize the client."""
        # Normalize password to SecretStr
        if isinstance(self.password, str):
            self.password = SecretStr(self.password)

        # Generate room ID if not provided
        if self.room is None:
            self.room = str(uuid.uuid4())

        # Ensure URL doesn't have trailing slash
        self.url = self.url.rstrip("/")

        # Create API manager
        self.api = APIManager(url=self.url, room_id=self.room)

        # Authenticate
        if self.user is not None:
            pw = self.password  # SecretStr | None after normalization above
            if pw is None:
                from zndraw.config import Settings

                pw = Settings().guest_password
            self.api.login(self.user, pw.get_secret_value())
        else:
            data = self.api.create_guest_session()
            self.user = data["email"]

        # Create socket manager
        self.socket = SocketManager(zndraw=self)

        # Create job manager (zero-cost until first register())
        self._jobs = JobManager(
            api=self.api,
            tsio=self.socket._tsio,
            execute=self._execute_task if self.auto_pickup else None,
            heartbeat_interval=self.heartbeat_interval,
            polling_interval=self.polling_interval,
        )

        # Auto-connect if requested
        if self.auto_connect:
            self.connect()

    # -------------------------------------------------------------------------
    # Connection Management
    # -------------------------------------------------------------------------

    def connect(self) -> None:
        """Connect to the server."""
        self.socket.connect()

    def disconnect(self) -> None:
        """Disconnect from the server."""
        self.socket.disconnect()
        self.api.close()

    def wait(self) -> None:
        """Block until disconnected."""
        self.socket.wait()

    @property
    def connected(self) -> bool:
        """Check if connected."""
        return self.socket.connected

    def __enter__(self) -> ZnDraw:
        """Context manager entry."""
        if not self.connected:
            self.connect()
        return self

    def __exit__(self, *args: Any) -> None:
        """Context manager exit."""
        self.disconnect()

    # -------------------------------------------------------------------------
    # Source Mount
    # -------------------------------------------------------------------------

    def mount(self, source: FrameSource) -> None:
        """Mount a virtual frame source on this room.

        The room must be empty (``len(vis) == 0``) and have no existing mount.
        After mounting, the room becomes read-only — frames are served on demand
        from the source via the provider system.

        Registers two providers:
        - ``FrameSourceRead`` (category "frames") — serves individual frames.
        - ``FrameSourceLength`` (category "frames_meta") — serves source length.

        Then PATCHes the room with the frame count so ``get_length()`` returns
        the correct value immediately.

        Parameters
        ----------
        source
            Any object satisfying the ``FrameSource`` protocol
            (e.g. ``znh5md.IO``, ``list[ase.Atoms]``).
        """
        from zndraw.providers.frame_source import FrameSourceLength, FrameSourceRead

        try:
            length = len(source)
        except TypeError:
            raise TypeError(
                "Frame source does not support len(). "
                "For file-based backends (e.g. asebytes.ASEIO with XYZ), "
                "call db._backend.count_frames() before mounting."
            ) from None

        if self._mount is not None:
            raise RuntimeError("Room already has a mount")
        if len(self) != 0:
            raise RuntimeError("Room must be empty before mounting")
        if self.room is None:
            raise NotConnectedError("Cannot mount: no room set")

        self._mount = source
        # UUID ensures each mount gets fresh provider cache keys —
        # prevents stale results from a previous mount being served.
        mount_name = f"mount-{uuid.uuid4().hex[:12]}"
        self._mount_name = mount_name
        self.register_provider(
            FrameSourceRead, name=mount_name, handler=source, room=self.room
        )
        self.register_provider(
            FrameSourceLength, name=mount_name, handler=source, room=self.room
        )
        self.api.update_room({"frame_count": length})
        self._cached_length = length

    def unmount(self) -> None:
        """Unmount the current source. Room returns to empty.

        Unregisters both providers and clears the external frame count.
        """
        if self._mount is None:
            raise RuntimeError("No mount to remove")

        if self.room is None:
            raise NotConnectedError("Cannot unmount: no room set")
        mount_name = self._mount_name
        self.jobs.unregister_provider(f"{self.room}:frames:{mount_name}")
        self.jobs.unregister_provider(f"{self.room}:frames_meta:{mount_name}")
        self.api.update_room({"frame_count": 0})
        self._mount = None
        self._mount_name = None
        self._cached_length = 0

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------

    @property
    def step(self) -> int:
        """Current frame index."""
        return self.api.get_step()["step"]

    @step.setter
    def step(self, value: int) -> None:
        """Set current frame index."""
        if value < 0:
            value = len(self) + value
        self.api.update_step(value)

    @property
    def selections(self) -> Selections:
        """Access selections by geometry name."""
        if self._selections is None:
            self._selections = Selections(self)
        return self._selections

    @property
    def selection_groups(self) -> SelectionGroups:
        """Access named selection groups."""
        if self._selection_groups is None:
            self._selection_groups = SelectionGroups(self)
        return self._selection_groups

    @property
    def selection(self) -> tuple[int, ...]:
        """Get selection for 'particles' geometry (convenience property)."""
        return self.selections.get("particles", ())

    @selection.setter
    def selection(self, value: Iterable[int] | None) -> None:
        """Set selection for 'particles' geometry."""
        self.selections["particles"] = [] if value is None else list(value)

    @property
    def frame_selection(self) -> tuple[int, ...]:
        """Selected frame indices (sorted)."""
        indices = self.api.get_frame_selection()
        return tuple(sorted(indices)) if indices else ()

    @frame_selection.setter
    def frame_selection(self, value: Iterable[int] | None) -> None:
        """Set selected frame indices."""
        self.api.update_frame_selection(sorted(value) if value else [])

    @property
    def atoms(self) -> ase.Atoms:
        """Get the current frame as an ase.Atoms object."""
        return self[self.step]

    @atoms.setter
    def atoms(self, value: ase.Atoms) -> None:
        """Set the current frame from an ase.Atoms object."""
        self[self.step] = value

    @property
    def bookmarks(self) -> Bookmarks:
        """Access frame bookmarks."""
        if self._bookmarks is None:
            self._bookmarks = Bookmarks(self)
        return self._bookmarks

    @property
    def geometries(self) -> Geometries:
        """Access custom geometries."""
        if self._geometries is None:
            self._geometries = Geometries(self)
        return self._geometries

    @property
    def figures(self) -> Figures:
        """Access Plotly figures."""
        if self._figures is None:
            self._figures = Figures(self)
        return self._figures

    @property
    def metadata(self) -> RoomMetadata:
        """Access room metadata."""
        if self._metadata is None:
            self._metadata = RoomMetadata(self)
        return self._metadata

    @property
    def sessions(self) -> Sessions:
        """Active frontend browser sessions for this user."""
        return Sessions(_api=self.api)

    def log(self, message: str) -> None:
        """Send a chat message to the room.

        Parameters
        ----------
        message
            The message content (supports Markdown).
        """
        self.api.create_chat_message(message)

    @property
    def jobs(self) -> JobManager:
        """Access the job manager for registering jobs and submitting tasks."""
        return self._jobs

    def _resolve_room(self, room: str | None) -> str:
        """Resolve room argument, defaulting to self.room."""
        resolved = room or self.room
        if resolved is None:
            raise NotConnectedError("No room set")
        return resolved

    def register_job(self, cls: type, *, room: str | None = None) -> None:
        """Register an extension as a job.

        Parameters
        ----------
        cls
            Extension subclass to register.
        room
            Room scope. Defaults to ``self.room``.
        """
        self.jobs.register(cls, room=self._resolve_room(room))

    def register_provider(
        self,
        provider_cls: type,
        *,
        name: str,
        handler: Any,
        room: str | None = None,
    ) -> uuid.UUID:
        """Register a provider for serving read requests.

        Parameters
        ----------
        provider_cls
            Provider subclass defining category and read schema.
        name
            Unique name for this provider instance.
        handler
            Object passed to ``provider.read(handler)`` on dispatch.
        room
            Room scope. Defaults to ``self.room``.
        """
        return self.jobs.register_provider(
            provider_cls, name=name, handler=handler, room=self._resolve_room(room)
        )

    def register_fs(self, fs: Any, *, name: str, room: str | None = None) -> None:
        """Register an fsspec filesystem as a provider with file loading.

        Parameters
        ----------
        fs
            An fsspec filesystem instance (e.g. ``fsspec.filesystem("file")``).
        name
            Unique name for this filesystem (e.g. "local", "s3-bucket").
        room
            Room scope. Defaults to ``self.room``.
        """
        from zndraw.extensions.filesystem import LoadFile
        from zndraw.providers.filesystem import FilesystemRead

        r = self._resolve_room(room)
        self.register_provider(FilesystemRead, name=name, handler=fs, room=r)
        self.register_job(LoadFile, room=r)

    # -------------------------------------------------------------------------
    # Locking
    # -------------------------------------------------------------------------

    def get_lock(self, msg: str | None = None) -> ZnDrawLock:
        """Get an edit lock context manager.

        Parameters
        ----------
        msg : str | None
            Optional message describing the lock purpose.

        Returns
        -------
        ZnDrawLock
            Lock context manager.

        Examples
        --------
        >>> with vis.get_lock(msg="Uploading frames"):
        ...     vis.extend(frames)
        """
        if not self.connected:
            raise NotConnectedError("Cannot acquire lock when not connected")
        return ZnDrawLock(api=self.api, msg=msg)

    # -------------------------------------------------------------------------
    # MutableSequence Implementation (Frame Operations)
    # -------------------------------------------------------------------------

    def __len__(self) -> int:
        """Return number of frames.

        Uses a cached value when available (seeded on connect, updated by
        ``FramesInvalidate`` socket events and local mutations). Falls back
        to an HTTP request on cache miss.
        """
        if self._cached_length is not None:
            return self._cached_length
        info = self.api.get_room_info()
        length = info.get("frame_count", 0)
        self._cached_length = length
        return length

    @overload
    def __getitem__(self, index: int) -> ase.Atoms: ...

    @overload
    def __getitem__(self, index: slice) -> list[ase.Atoms]: ...

    def __getitem__(self, index: int | slice) -> ase.Atoms | list[ase.Atoms]:
        """Get frame(s) by index or slice as ase.Atoms objects."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            frame = self.api.get_frame(index)
            return raw_frame_to_atoms(frame)

        if isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            if step == 1:
                frames = self.api.get_frames(start=start, stop=stop)
            else:
                indices = list(range(start, stop, step))
                frames = self.api.get_frames(indices=indices) if indices else []
            return [raw_frame_to_atoms(f) for f in frames]

        raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    @overload
    def __setitem__(self, index: int, value: ase.Atoms) -> None: ...

    @overload
    def __setitem__(self, index: slice, value: Iterable[ase.Atoms]) -> None: ...

    def __setitem__(  # type: ignore[override]
        self, index: int | slice, value: ase.Atoms | Iterable[ase.Atoms]
    ) -> None:
        """Set frame(s) at given index from ase.Atoms objects."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            if not isinstance(value, ase.Atoms):
                raise TypeError("Value must be an ase.Atoms object")
            self.api.update_frame(index, atoms_to_json_dict(value))

        elif isinstance(index, slice):
            if isinstance(value, ase.Atoms):
                raise TypeError(
                    "Value must be iterable of ase.Atoms for slice assignment"
                )
            value_list = list(value)
            if not all(isinstance(v, ase.Atoms) for v in value_list):
                raise TypeError("All values must be ase.Atoms objects")
            length = len(self)
            start, stop, step = index.indices(length)
            indices = list(range(start, stop, step))

            if step != 1 and len(value_list) != len(indices):
                raise ValueError(
                    f"attempt to assign sequence of size {len(value_list)} "
                    f"to extended slice of size {len(indices)}"
                )

            for i, idx in enumerate(indices):
                if i < len(value_list):
                    self.api.update_frame(idx, atoms_to_json_dict(value_list[i]))

        else:
            raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    def __delitem__(self, index: int | slice) -> None:
        """Delete frame(s) at given index."""
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            self.api.delete_frame(index)
            self._cached_length = None

        elif isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            indices = sorted(range(start, stop, step), reverse=True)
            for idx in indices:
                self.api.delete_frame(idx)
            self._cached_length = None

        else:
            raise TypeError(f"Index must be int or slice, not {type(index).__name__}")

    def insert(self, index: int, value: ase.Atoms) -> None:
        """Insert an ase.Atoms frame at the given index.

        Note: Due to API limitations, this appends and then reorders.
        """
        if not isinstance(value, ase.Atoms):
            raise TypeError("Value must be an ase.Atoms object")
        # For simplicity, just append (true insert requires server support)
        result = self.api.append_frames([atoms_to_json_dict(value)])
        self._cached_length = result.get("total")

    def append(self, value: ase.Atoms) -> None:
        """Append an ase.Atoms frame to the end."""
        if not isinstance(value, ase.Atoms):
            raise TypeError("Value must be an ase.Atoms object")
        result = self.api.append_frames([atoms_to_json_dict(value)])
        self._cached_length = result.get("total")

    @contextlib.contextmanager
    @typing_extensions.deprecated(
        "Use ZnDrawTqdm directly: "
        "for x in ZnDrawTqdm(items, vis=vis, description='...', unit='it'): ..."
    )
    def progress_bar(
        self,
        iterable: Iterable[Any] | None = None,
        *,
        total: int | None = None,
        description: str = "Processing...",
        unit: str = "it",
    ) -> Generator[ZnDrawTqdm, None, None]:
        """Context manager yielding a ZnDrawTqdm progress bar.

        .. deprecated::
            Use :class:`ZnDrawTqdm` directly instead.

        Parameters
        ----------
        iterable : Iterable[Any] | None
            Optional iterable to wrap.
        total : int | None
            Total expected iterations.
        description : str
            Label shown in the UI.
        unit : str
            Unit label (e.g. ``"frames"``).

        Yields
        ------
        ZnDrawTqdm
            A tqdm-compatible progress bar.
        """
        from zndraw.tqdm import ZnDrawTqdm

        pbar = ZnDrawTqdm(
            iterable, total=total, vis=self, description=description, unit=unit
        )
        try:
            yield pbar
        finally:
            pbar.close()

    def extend(self, values: Iterable[ase.Atoms]) -> None:
        """Extend with multiple ase.Atoms frames.

        Streams frames in size-targeted chunks (~2 MB, max 1000 frames each).
        Shows a tqdm progress bar in the terminal and broadcasts progress to
        the ZnDraw UI via ``ZnDrawTqdm``. Accepts any iterable including
        generators.
        """
        from zndraw.tqdm import ZnDrawTqdm

        try:
            total_frames: int | None = len(values)  # type: ignore[arg-type]
        except TypeError:
            total_frames = None

        chunk: list[dict[str, Any]] = []
        chunk_size = 0
        result: dict[str, Any] = {}

        progress = ZnDrawTqdm(
            total=total_frames,
            vis=self,
            description="Uploading frames",
            unit="frames",
        )

        try:
            for atoms in values:
                if not isinstance(atoms, ase.Atoms):
                    raise TypeError("All values must be ase.Atoms objects")
                frame = atoms_to_json_dict(atoms)
                frame_size = _estimate_frame_size(frame)

                if chunk_size > 0 and (
                    chunk_size + frame_size > _TARGET_CHUNK_BYTES
                    or len(chunk) >= _MAX_CHUNK_FRAMES
                ):
                    result = self.api.append_frames(chunk)
                    progress.update(len(chunk))
                    chunk = []
                    chunk_size = 0

                chunk.append(frame)
                chunk_size += frame_size

            # Flush remaining frames
            if chunk:
                result = self.api.append_frames(chunk)
                progress.update(len(chunk))
        finally:
            progress.close()

        if result:
            self._cached_length = result.get("total")

    @staticmethod
    def _decode_raw_frame(frame: dict[bytes, bytes]) -> dict[str, Any]:
        """Decode a raw msgpack frame to a dict with string keys and decoded values."""
        import msgpack_numpy

        return {
            k.decode(): msgpack.unpackb(v, object_hook=msgpack_numpy.decode)
            for k, v in frame.items()
        }

    def get(
        self,
        index: int | list[int] | slice,
        keys: list[str] | None = None,
    ) -> dict[str, Any] | list[dict[str, Any]]:
        """Get decoded frame data with optional key filtering.

        Returns decoded dictionaries with string keys and Python/numpy values,
        transferring only the requested keys from the server.

        Parameters
        ----------
        index : int | list[int] | slice
            Frame index, list of indices, or slice.
        keys : list[str] | None
            If provided, only return these keys from each frame
            (e.g. ``["info.energy", "arrays.positions"]``).

        Returns
        -------
        dict[str, Any] | list[dict[str, Any]]
            Single dict for int index, list of dicts for slice/list index.
        """
        if isinstance(index, int):
            if index < 0:
                index = len(self) + index
            frames = self.api.get_frames(indices=[index], keys=keys)
            return self._decode_raw_frame(frames[0]) if frames else {}

        if isinstance(index, list):
            length = len(self)
            normalized = [i if i >= 0 else length + i for i in index]
            frames = self.api.get_frames(indices=normalized, keys=keys)
        elif isinstance(index, slice):
            length = len(self)
            start, stop, step = index.indices(length)
            if step == 1:
                frames = self.api.get_frames(start=start, stop=stop, keys=keys)
            else:
                indices = list(range(start, stop, step))
                frames = self.api.get_frames(indices=indices, keys=keys)
        else:
            raise TypeError(
                f"Index must be int, list, or slice, not {type(index).__name__}"
            )

        return [self._decode_raw_frame(f) for f in frames]

    def set_frames(
        self,
        index: int | slice | list[int],
        value: ase.Atoms | list[ase.Atoms],
    ) -> None:
        """Set frame(s) at given index from ase.Atoms objects.

        This is an alternative to __setitem__ that supports list indices.
        """
        if not self.connected:
            raise NotConnectedError("Client is not connected")

        if isinstance(index, int):
            if not isinstance(value, ase.Atoms):
                raise TypeError("Value must be ase.Atoms for single index")
            self[index] = value
        elif isinstance(index, slice):
            if isinstance(value, ase.Atoms):
                raise TypeError("Value must be list of ase.Atoms for slice index")
            self[index] = value
        elif isinstance(index, list):
            if not isinstance(value, list):
                raise TypeError("Value must be list of ase.Atoms for list index")
            if not all(isinstance(v, ase.Atoms) for v in value):
                raise TypeError("All values must be ase.Atoms objects")
            # For list indices, update each frame individually
            if len(index) != len(value):
                raise ValueError("Index and value lists must have same length")
            for idx, atoms in zip(index, value, strict=True):
                self.api.update_frame(idx, atoms_to_json_dict(atoms))
        else:
            raise TypeError(
                f"Index must be int, slice, or list, not {type(index).__name__}"
            )

    # -------------------------------------------------------------------------
    # Extension Execution
    # -------------------------------------------------------------------------

    def run(self, extension: Extension, *, job_room: str | None = None) -> str:
        """Submit an extension for execution via the job system.

        Parameters
        ----------
        extension : Extension
            The extension instance to run.
        job_room : str | None
            The room where the job is registered. If None, auto-detects:
            ``@internal`` for built-in extensions, ``@global`` otherwise.

        Returns
        -------
        str
            The task ID for tracking progress.
        """
        if job_room is None:
            module = extension.__class__.__module__
            job_room = (
                "@internal" if module.startswith("zndraw.extensions") else "@global"
            )

        if self.room is None:
            raise NotConnectedError("Cannot run: no room set")
        return self.jobs.submit(
            cast("JoblibExtension", extension), room=self.room, job_room=job_room
        )

    # -------------------------------------------------------------------------
    # Worker Serve Loop
    # -------------------------------------------------------------------------

    def _execute_task(self, task: ClaimedTask) -> None:
        """Execute a claimed task's extension logic.

        Lifecycle management (start/complete/fail) is handled by
        ``JobManager._claim_loop`` -- this callback only runs the extension.
        """
        task_vis = ZnDraw(
            url=self.url,
            room=task.room_id,
            user=self.user,
            password=self.password,
        )
        try:
            task.extension.run(task_vis, providers=self.jobs.handlers)
        finally:
            task_vis.disconnect()
