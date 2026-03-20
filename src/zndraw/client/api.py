"""REST API manager for the ZnDraw client."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import httpx
import msgpack

from zndraw.exceptions import PROBLEM_TYPES, ProblemDetail, RoomLockedError, ZnDrawError

if TYPE_CHECKING:
    from zndraw.schemas import SessionItem

log = logging.getLogger(__name__)


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
    lock_token: str | None = field(default=None, init=False)

    http: httpx.Client = field(init=False)

    def __post_init__(self) -> None:
        """Initialize HTTP client."""
        self.http = httpx.Client(base_url=self.url, timeout=30.0)

    @property
    def base_url(self) -> str:
        """Base URL (satisfies ApiManager protocol)."""
        return self.url

    def get_headers(self) -> dict[str, str]:
        """Build request headers."""
        headers: dict[str, str] = {}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        if self.session_id:
            headers["X-Session-ID"] = self.session_id
        if self.lock_token:
            headers["Lock-Token"] = self.lock_token
        return headers

    def close(self) -> None:
        """Close the HTTP client."""
        self.http.close()

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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_room_info(self) -> dict[str, Any]:
        """Get room information."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
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

        from zndraw_joblib.exceptions import ProviderTimeoutError

        while True:
            response = self.http.get(
                f"/v1/rooms/{self.room_id}/frames/{index}",
                headers=self.get_headers(),
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

        from zndraw_joblib.exceptions import ProviderTimeoutError

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
                headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_frame(self, index: int, data: dict[str, Any]) -> dict[str, Any]:
        """Update a single frame."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/frames/{index}",
            json={"data": data},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_frame(self, index: int) -> None:
        """Delete a single frame."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/frames/{index}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Step Operations
    # -------------------------------------------------------------------------

    def get_step(self) -> dict[str, Any]:
        """Get current step."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/step",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_step(self, step: int) -> dict[str, Any]:
        """Update current step."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/step",
            json={"step": step},
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["frame_selection"]

    def update_frame_selection(self, indices: list[int]) -> dict[str, Any]:
        """Update selected frame indices."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/frame-selection",
            json={"indices": indices},
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_geometry(self, key: str) -> dict[str, Any] | None:
        """Get a specific geometry."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/geometries/{key}",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_geometry(self, key: str) -> None:
        """Delete a geometry."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/geometries/{key}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Default Camera
    # -------------------------------------------------------------------------

    def get_default_camera(self) -> str | None:
        """Get the default camera key for the room."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/default-camera",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["default_camera"]

    def set_default_camera(self, camera_key: str | None) -> None:
        """Set or unset the default camera."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/default-camera",
            json={"default_camera": camera_key},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Selection Operations
    # -------------------------------------------------------------------------

    def get_selection(self, geometry: str) -> dict[str, Any]:
        """Get selection for a specific geometry."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/geometries/{geometry}/selection",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_selection(self, geometry: str, indices: list[int]) -> dict[str, Any]:
        """Update selection for a geometry."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/geometries/{geometry}/selection",
            json={"indices": indices},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def list_selection_groups(self) -> dict[str, dict[str, list[int]]]:
        """List all selection groups."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/selection-groups",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_selection_group(self, group_name: str) -> dict[str, Any]:
        """Get a selection group."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/selection-groups/{group_name}",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_selection_group(self, group_name: str) -> None:
        """Delete a selection group."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/selection-groups/{group_name}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Bookmark Operations
    # -------------------------------------------------------------------------

    def get_all_bookmarks(self) -> dict[str, str]:
        """Get all bookmarks."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/bookmarks",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_bookmark(self, index: int) -> dict[str, Any]:
        """Get a specific bookmark."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def set_bookmark(self, index: int, label: str) -> dict[str, Any]:
        """Set a bookmark."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            json={"label": label},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_bookmark(self, index: int) -> None:
        """Delete a bookmark."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/bookmarks/{index}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Figure Operations
    # -------------------------------------------------------------------------

    def list_figures(self) -> list[str]:
        """List all figure keys."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/figures",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_figure(self, key: str) -> dict[str, Any] | None:
        """Get a specific figure."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/figures/{key}",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_figure(self, key: str) -> None:
        """Delete a figure."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/figures/{key}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Preset Operations
    # -------------------------------------------------------------------------

    def list_presets(self) -> list[dict[str, Any]]:
        """List all presets."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/presets",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    def get_preset(self, name: str) -> dict[str, Any] | None:
        """Get a specific preset."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/presets/{name}",
            headers=self.get_headers(),
        )
        if response.status_code == 404:
            return None
        self.raise_for_status(response)
        return response.json()

    def create_preset(self, data: dict[str, Any]) -> dict[str, Any]:
        """Create a new preset."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/presets",
            json=data,
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def update_preset(self, name: str, data: dict[str, Any]) -> dict[str, Any]:
        """Create or update a preset."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/presets/{name}",
            json=data,
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def delete_preset(self, name: str) -> None:
        """Delete a preset."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/presets/{name}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    def apply_preset(self, name: str) -> dict[str, Any]:
        """Apply a preset to all matching geometries."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/presets/{name}/apply",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Edit Lock Operations
    # -------------------------------------------------------------------------

    def edit_lock_acquire(self, msg: str | None = None) -> dict[str, Any]:
        """Acquire the room edit lock."""
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/edit-lock",
            json={"msg": msg},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def edit_lock_refresh(
        self, lock_token: str, msg: str | None = None
    ) -> dict[str, Any]:
        """Refresh an existing edit lock using its token."""
        headers = {**self.get_headers(), "Lock-Token": lock_token}
        response = self.http.put(
            f"/v1/rooms/{self.room_id}/edit-lock",
            json={"msg": msg},
            headers=headers,
        )
        self.raise_for_status(response)
        return response.json()

    def edit_lock_release(self, lock_token: str | None = None) -> None:
        """Release the room edit lock."""
        headers = self.get_headers()
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Auth Operations
    # -------------------------------------------------------------------------

    def get_me(self) -> dict[str, Any]:
        """Get the current authenticated user's profile."""
        response = self.http.get("/v1/auth/users/me", headers=self.get_headers())
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Session Operations
    # -------------------------------------------------------------------------

    def list_sessions(self) -> list[SessionItem]:
        """List all active frontend sessions in the room."""
        from zndraw.schemas import SessionsListResponse

        response = self.http.get(
            f"/v1/rooms/{self.room_id}/sessions",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return SessionsListResponse.model_validate(response.json()).items

    def get_active_camera(self, sid: str) -> str:
        """Get active camera key for a session."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/sessions/{sid}/active-camera",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
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
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_screenshot(self, screenshot_id: int) -> dict[str, Any]:
        """Get a screenshot by ID."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/screenshots/{screenshot_id}",
            headers=self.get_headers(),
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
            headers=self.get_headers(),
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
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def progress_complete(self, progress_id: str) -> None:
        """Complete and remove a progress tracker."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/progress/{progress_id}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)

    # -------------------------------------------------------------------------
    # Room Listing (no room_id needed)
    # -------------------------------------------------------------------------

    def list_rooms(self, search: str | None = None) -> list[dict[str, Any]]:
        """List all rooms, optionally filtered by search query."""
        params: dict[str, str] = {}
        if search is not None:
            params["search"] = search
        response = self.http.get(
            "/v1/rooms",
            params=params,
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    # -------------------------------------------------------------------------
    # Chat Listing
    # -------------------------------------------------------------------------

    def list_chat_messages(
        self,
        limit: int | None = None,
        before: str | None = None,
    ) -> dict[str, Any]:
        """List chat messages with optional pagination."""
        params: dict[str, Any] = {}
        if limit is not None:
            params["limit"] = limit
        if before is not None:
            params["before"] = before
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/chat/messages",
            params=params,
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Screenshot Listing
    # -------------------------------------------------------------------------

    def list_screenshots(self) -> list[dict[str, Any]]:
        """List screenshots for the room."""
        response = self.http.get(
            f"/v1/rooms/{self.room_id}/screenshots",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()["items"]

    # -------------------------------------------------------------------------
    # Extensions (joblib jobs)
    # -------------------------------------------------------------------------

    def list_extensions(self, room: str | None = None) -> dict[str, Any]:
        """List extensions for a job room."""
        job_room = room or "@internal"
        response = self.http.get(
            f"/v1/joblib/rooms/{job_room}/jobs",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_extension(self, full_name: str) -> dict[str, Any]:
        """Get extension details by full_name."""
        # full_name is like "@internal:modifiers:Delete"
        # The job_room is the first segment
        parts = full_name.split(":", 1)
        if len(parts) != 2:
            raise ValueError(f"Invalid extension name: {full_name}")
        job_room = parts[0]
        response = self.http.get(
            f"/v1/joblib/rooms/{job_room}/jobs/{full_name}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def submit_task(self, full_name: str, payload: dict[str, Any]) -> dict[str, Any]:
        """Submit a task for an extension."""
        response = self.http.post(
            f"/v1/joblib/rooms/{self.room_id}/tasks/{full_name}",
            json={"payload": payload},
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    # -------------------------------------------------------------------------
    # Tasks
    # -------------------------------------------------------------------------

    def list_tasks(self, status: str | None = None) -> dict[str, Any]:
        """List tasks for the room."""
        params: dict[str, Any] = {}
        if status is not None:
            params["status"] = status
        response = self.http.get(
            f"/v1/joblib/rooms/{self.room_id}/tasks",
            params=params,
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def get_task(self, task_id: str) -> dict[str, Any]:
        """Get task details by ID."""
        response = self.http.get(
            f"/v1/joblib/tasks/{task_id}",
            headers=self.get_headers(),
        )
        self.raise_for_status(response)
        return response.json()
