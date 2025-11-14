import dataclasses
import typing as t

import msgpack
import msgpack_numpy as m
import requests

from zndraw.lock import ZnDrawLock


@dataclasses.dataclass
class APIManager:
    url: str
    room: str
    user_name: str | None = None
    jwt_token: str | None = None
    session_id: str | None = None  # Session ID from /join response

    def login(self, user_name: str | None, password: str | None = None) -> dict:
        """Authenticate and get JWT token.

        Parameters
        ----------
        user_name : str | None
            Username for authentication. If None, server assigns a guest username.
        password : str | None
            Optional password for admin authentication (deployment mode only)

        Returns
        -------
        dict
            {"status": str, "token": str, "userName": str, "role": str}

        Raises
        ------
        RuntimeError
            If login fails
        """
        payload = {"userName": user_name}
        if password is not None:
            payload["password"] = password

        response = requests.post(f"{self.url}/api/login", json=payload)

        if response.status_code != 200:
            raise RuntimeError(f"Login failed: {response.text}")

        data = response.json()
        # Update internal state
        self.jwt_token = data["token"]
        self.user_name = data["userName"]
        return data

    def _get_headers(self) -> dict:
        """Build headers for API requests.

        Returns
        -------
        dict
            Headers dictionary with Authorization and X-Session-ID
        """
        headers = {}
        if self.jwt_token:
            headers["Authorization"] = f"Bearer {self.jwt_token}"
        if self.session_id:
            headers["X-Session-ID"] = self.session_id
        return headers

    def _raise_for_error_type(self, response: requests.Response) -> None:
        """Convert API error responses to appropriate Python exceptions.

        Parameters
        ----------
        response : requests.Response
            The HTTP response object

        Raises
        ------
        KeyError, IndexError, ValueError, TypeError, PermissionError
            Based on the error type returned by the API
        """
        try:
            error_data = response.json()
            error_type = error_data.get("type", "")
            error_msg = error_data.get("error", response.text)

            exception_map = {
                "KeyError": KeyError,
                "IndexError": IndexError,
                "ValueError": ValueError,
                "TypeError": TypeError,
                "PermissionError": PermissionError,
            }

            if error_type in exception_map:
                raise exception_map[error_type](error_msg)
        except ValueError:
            # JSON decode failed, let raise_for_status handle it
            pass

    def get_version(self) -> str:
        """Get the server version.

        Returns
        -------
        str
            The server version string.

        Raises
        ------
        RuntimeError
            If the version endpoint is not accessible.
        """
        response = requests.get(f"{self.url}/api/version", timeout=5.0)
        response.raise_for_status()
        return response.json()["version"]

    def join_room(
        self,
        description: str | None = None,
        copy_from: str | None = None,
    ) -> dict:
        """Join a room, optionally creating it with a description or copying from an existing room.

        Args:
            description: Optional description for the room (only used if room is created)
            copy_from: Optional room ID to copy frames and settings from (only used if room is created)

        Returns:
            Dict containing room information (userName comes from JWT token)
        """
        payload = {}

        if description is not None:
            payload["description"] = description
        if copy_from is not None:
            payload["copyFrom"] = copy_from

        # Send JWT token in Authorization header
        headers = {}
        if self.jwt_token:
            headers["Authorization"] = f"Bearer {self.jwt_token}"

        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/join", json=payload, headers=headers
        )

        if response.status_code != 200:
            raise RuntimeError(
                f"Failed to join room '{self.room}': {response.status_code} {response.text}"
            )

        data = response.json()
        # Store session ID for subsequent requests
        self.session_id = data.get("sessionId")
        return data

    def get_job(self, job_id: str) -> dict:
        """Get job details via REST API.

        Parameters
        ----------
        job_id : str
            Job ID

        Returns
        -------
        dict
            Job details with keys: jobId, room, category, extension, data, public, status, etc.

        Raises
        ------
        requests.HTTPError
            If job not found (404) or request fails
        """
        headers = self._get_headers()
        url = f"{self.url}/api/jobs/{job_id}"
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        return response.json()

    def update_job_status(
        self,
        job_id: str,
        status: str,
        worker_id: str,
        error: str | None = None,
    ) -> None:
        """Update job status via REST API.

        Parameters
        ----------
        job_id : str
            Job ID
        status : str
            New status (processing, completed, failed)
        worker_id : str
            Worker ID
        error : str | None
            Error message (for failed status)

        Raises
        ------
        requests.HTTPError
            If status update fails
        """
        headers = self._get_headers()

        payload = {"status": status, "workerId": worker_id}
        if error is not None:
            payload["error"] = error

        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/jobs/{job_id}/status",
            json=payload,
            headers=headers,
        )
        response.raise_for_status()

    def register_extension(
        self, name: str, category: str, schema: dict, socket_manager, public: bool = False
    ) -> None:
        """Register extension via REST endpoint.

        This uses the sessionId to identify the worker, which the server
        maps to the socket sid for consistent worker tracking.

        Parameters
        ----------
        name : str
            Extension name
        category : str
            Extension category
        schema : dict
            JSON schema for the extension
        socket_manager : SocketManager
            Socket manager (unused, kept for API compatibility)
        public : bool
            If True, register as global extension (requires admin privileges)
        """
        headers = self._get_headers()

        payload = {
            "sessionId": self.session_id,
            "roomId": self.room,
            "name": name,
            "category": category,
            "schema": schema,
            "public": public,
        }

        response = requests.post(
            f"{self.url}/api/workers/register",
            json=payload,
            headers=headers,
            timeout=10,
        )

        if response.status_code != 200:
            try:
                error_data = response.json()
                error = error_data.get("error", response.text)
            except Exception:
                error = response.text
            raise RuntimeError(f"Extension registration failed: {error}")

        data = response.json()
        if not data.get("success"):
            raise RuntimeError(f"Extension registration failed: {data.get('error', 'Unknown error')}")

        # Return the worker_id assigned by server so caller can store it
        return data.get("workerId")

    def register_filesystem(
        self,
        name: str,
        fs_type: str,
        socket_manager,
        public: bool = False,
    ) -> None:
        """Register filesystem via REST endpoint.

        This uses the sessionId to identify the worker, which the server
        maps to the socket sid for consistent worker tracking.

        Parameters
        ----------
        name : str
            Filesystem name
        fs_type : str
            Filesystem class name (e.g., 'LocalFileSystem', 'S3FileSystem')
        socket_manager : SocketManager
            Socket manager (unused, kept for API compatibility)
        public : bool
            If True, register as global filesystem (requires admin privileges)
        """
        headers = self._get_headers()

        payload = {
            "sessionId": self.session_id,
            "roomId": self.room,
            "name": name,
            "fsType": fs_type,
            "public": public,
        }

        response = requests.post(
            f"{self.url}/api/workers/filesystem/register",
            json=payload,
            headers=headers,
            timeout=10,
        )

        if response.status_code != 200:
            try:
                error_data = response.json()
                error = error_data.get("error", response.text)
            except Exception:
                error = response.text
            raise RuntimeError(f"Filesystem registration failed: {error}")

        data = response.json()
        if not data.get("success"):
            raise RuntimeError(
                f"Filesystem registration failed: {data.get('error', 'Unknown error')}"
            )

        # Return the worker_id assigned by server so caller can store it
        return data.get("workerId")

    def get_frames(self, indices_or_slice, keys: list[str] | None = None) -> list[dict[bytes, bytes]]:
        """Get frames as raw dict[bytes, bytes] format.

        Returns frames in their serialized form - decoding to ase.Atoms
        should be done on-demand by the caller.
        """
        if isinstance(indices_or_slice, list):
            payload = {"indices": ",".join(str(i) for i in indices_or_slice)}
        elif isinstance(indices_or_slice, slice):
            payload = {}
            if indices_or_slice.start is not None:
                payload["start"] = indices_or_slice.start
            if indices_or_slice.stop is not None:
                payload["stop"] = indices_or_slice.stop
            if indices_or_slice.step is not None:
                payload["step"] = indices_or_slice.step
        else:
            raise ValueError(
                "indices_or_slice must be either a list of integers or a slice object"
            )

        if keys is not None:
            payload["keys"] = ",".join(keys)

        headers = self._get_headers()
        full_url = f"{self.url}/api/rooms/{self.room}/frames"
        response = requests.get(full_url, params=payload, headers=headers, timeout=30)

        if response.status_code == 404:
            try:
                error_data = response.json()
                error_type = error_data.get("type", "")
                error_msg = error_data.get("error", response.text)

                if error_type == "KeyError":
                    raise KeyError(error_msg)
                elif error_type == "IndexError":
                    raise IndexError(error_msg)
            except ValueError:
                pass

        response.raise_for_status()

        # Deserialize msgpack response - returns list[dict[bytes, bytes]]
        # The server sends msgpack-encoded list of frames, where each frame
        # is dict[bytes, bytes] with msgpack-encoded values
        serialized_frames = msgpack.unpackb(
            response.content, strict_map_key=False
        )

        # Return raw dict[bytes, bytes] format - no decoding
        return serialized_frames

    def upload_frames(self, action: str, data, lock_token: str | None = None, **kwargs) -> dict:
        try:
            # Serialize data with msgpack-numpy support for numpy arrays
            # This handles both dict[bytes, bytes] from encode() and raw dicts with numpy arrays
            packed_data = msgpack.packb(
                data if isinstance(data, dict) else data, default=m.encode
            )

            upload_url = f"{self.url}/api/rooms/{self.room}/frames"
            params = {"action": action}
            params.update(kwargs)

            headers = self._get_headers()
            http_response = requests.post(
                upload_url, data=packed_data, params=params, headers=headers, timeout=30
            )

            if http_response.status_code == 404:
                try:
                    error_data = http_response.json()
                    error_type = error_data.get("type", "")
                    error_msg = error_data.get("error", http_response.text)

                    if error_type == "IndexError":
                        raise IndexError(error_msg)
                except ValueError:
                    pass

            http_response.raise_for_status()
            result = http_response.json()

            if not result.get("success"):
                raise RuntimeError(f"Server reported failure: {result.get('error')}")

            return result
        except requests.exceptions.RequestException as e:
            raise RuntimeError(f"Error uploading frame data: {e}") from e

    def delete_frames(self, index: int | slice | list[int]):
        delete_url = f"{self.url}/api/rooms/{self.room}/frames"
        params = {"action": "delete"}

        if isinstance(index, int):
            params["frame_id"] = index
        elif isinstance(index, list):
            params["indices"] = ",".join(str(i) for i in index)
        elif isinstance(index, slice):
            if index.start is not None:
                params["start"] = index.start
            if index.stop is not None:
                params["stop"] = index.stop
            if index.step is not None:
                params["step"] = index.step

        # Add userName if available (for lock checking)
        if self.user_name:
            params["userName"] = self.user_name

        headers = self._get_headers()
        response = requests.delete(delete_url, params=params, headers=headers, timeout=30)

        if response.status_code == 404:
            try:
                error_data = response.json()
                error_type = error_data.get("type", "")
                error_msg = error_data.get("error", response.text)

                if error_type == "IndexError":
                    raise IndexError(error_msg)
                elif error_type == "PermissionError":
                    raise PermissionError(error_msg)
            except ValueError:
                pass

        response.raise_for_status()
        return response.json()

    def bulk_patch_frames(
        self, data: list, start: int = None, stop: int = None, indices: list[int] = None
    ):
        # Serialize data with msgpack-numpy support for numpy arrays
        # This handles both dict[bytes, bytes] from encode() and raw dicts with numpy arrays
        packed_data = msgpack.packb(data, default=m.encode)

        bulk_url = f"{self.url}/api/rooms/{self.room}/frames/bulk"
        params = {}
        if start is not None and stop is not None:
            params = {"start": start, "stop": stop}
        elif indices is not None:
            params = {"indices": ",".join(str(i) for i in indices)}
        else:
            raise ValueError("Either start/stop or indices must be provided.")

        headers = self._get_headers()
        response = requests.patch(bulk_url, data=packed_data, params=params, headers=headers, timeout=30)

        if response.status_code == 404:
            try:
                error_data = response.json()
                error_type = error_data.get("type", "")
                error_msg = error_data.get("error", response.text)

                if error_type == "IndexError":
                    raise IndexError(error_msg)
            except ValueError:
                pass

        response.raise_for_status()

    def get_extension_settings(self, extension_name: str) -> dict:
        headers = self._get_headers()

        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/extensions/settings/{extension_name}/data",
            headers=headers,
        )
        response.raise_for_status()
        return response.json()

    def submit_extension_settings(self, extension_name: str, data: dict) -> None:
        headers = self._get_headers()

        # Settings extensions are always server-side (Celery) extensions,
        # so they use the public endpoint
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/public/settings/{extension_name}/submit",
            json=data,
            headers=headers,
        )
        response.raise_for_status()

    def run_extension(self, category: str, name: str, data: dict, public: bool = False, auto_retry: bool = False) -> dict:
        headers = self._get_headers()

        # Route to correct endpoint based on public flag
        scope = "public" if public else "private"
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/{scope}/{category}/{name}/submit",
            json={"data": data},
            headers=headers,
        )
        if response.status_code != 200:
            error_data = response.json()
            error_type = error_data.get("type", "")
            error_msg = error_data.get("error", response.text)
            error_code = error_data.get("code", "")

            # Auto-retry with public=True if this is a server-side extension
            # Only retry if auto_retry=True (user didn't explicitly specify public parameter)
            if error_code == "WRONG_ENDPOINT" and not public and auto_retry:
                return self.run_extension(category, name, data, public=True, auto_retry=False)

            # Convert EXTENSION_NOT_FOUND to ValueError with clearer message
            if error_code == "EXTENSION_NOT_FOUND":
                namespace = error_data.get("namespace", "unknown")
                extension_name = error_data.get("extension", name)
                # Use "public" and "private" in error messages to match test expectations
                namespace_display = "public" if namespace == "global" else "private"
                raise ValueError(
                    f"Extension '{extension_name}' is not registered in {namespace_display} namespace. "
                    f"Register it with vis.register_extension({extension_name}, public={public})"
                )

            if error_type == "KeyError":
                raise KeyError(error_msg)
            elif error_type == "ValueError":
                raise ValueError(error_msg)
            elif error_type == "TypeError":
                raise TypeError(error_msg)
            else:
                raise RuntimeError(f"Error running extension: {error_msg}")
        response.raise_for_status()
        return response.json()

    def get_messages(self, limit: int, before: int | None, after: int | None) -> dict:
        params = {"limit": limit}
        if before:
            params["before"] = before
        if after:
            params["after"] = after

        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/chat/messages", params=params, headers=headers
        )
        response.raise_for_status()
        return response.json()

    def set_geometry(self, data: dict, key: str, geometry_type: str) -> None:
        """Create or update a geometry.

        Args:
            data: The geometry data dict
            key: The geometry key/name
            geometry_type: The type of geometry (Sphere, Arrow, Bond, Curve)
        """
        with ZnDrawLock(api=self, target="trajectory:meta", msg=None) as lock:
            headers = self._get_headers()

            response = requests.post(
                f"{self.url}/api/rooms/{self.room}/geometries",
                json={"key": key, "data": data, "type": geometry_type},
                headers=headers,
            )
            response.raise_for_status()

    def get_geometry(self, key: str) -> dict | None:
        """Get a specific geometry by key.

        Args:
            key: The geometry key/name

        Returns:
            The geometry dict with 'type' and 'data' keys, or None if not found
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/geometries/{key}",
            headers=headers,
        )
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json().get("geometry", None)

    def del_geometry(self, key: str) -> None:
        """Delete a geometry.

        Args:
            key: The geometry key/name
        """
        with ZnDrawLock(api=self, target="trajectory:meta", msg=None) as lock:
            headers = self._get_headers()

            response = requests.delete(
                f"{self.url}/api/rooms/{self.room}/geometries/{key}", headers=headers
            )
            response_data = response.json()
            if response.status_code != 200:
                error = response_data.get("error", None)
                error_type = response_data.get("type", None)
                if error_type == "KeyError":
                    raise KeyError(error)
            response.raise_for_status()

    def list_geometries(self) -> dict:
        """Get all geometries with their full data.

        Returns:
            Dict mapping geometry keys to their data {key: {type, data}, ...}
            Example: {"particles": {"type": "Sphere", "data": {...}}}
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/geometries",
            headers=headers,
        )
        response.raise_for_status()
        return response.json().get("geometries", {})

    def get_geometries(self) -> dict | None:
        """Get all geometries (DEPRECATED: use list_geometries).

        This is now just an alias for list_geometries() since the endpoint
        was updated to return full geometry data in a single request.

        Returns:
            Dict mapping geometry keys to their data {key: {type, data}, ...}
        """
        return self.list_geometries()

    def add_figure(self, key: str, figure: dict) -> None:
        headers = self._get_headers()
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/figures",
            json={"key": key, "figure": figure},
            headers=headers,
        )
        response.raise_for_status()

    def delete_figure(self, key: str) -> None:
        headers = self._get_headers()
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/figures/{key}",
            headers=headers,
        )
        response_data = response.json()
        if response.status_code != 200:
            error = response_data.get("error", None)
            error_type = response_data.get("type", None)
            if error_type == "KeyError":
                raise KeyError(error)

        response.raise_for_status()

    def get_figure(self, key: str) -> dict | None:
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/figures/{key}",
            headers=headers,
        )
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json().get("figure", None)

    def list_figures(self) -> list[str]:
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/figures",
            headers=headers,
        )
        response.raise_for_status()
        return response.json().get("figures", [])

    def get_metadata(self) -> dict[str, str]:
        """Get all metadata for the room.

        Returns
        -------
        dict[str, str]
            Dictionary of metadata fields.
        """
        headers = self._get_headers()
        response = requests.get(f"{self.url}/api/rooms/{self.room}/metadata", headers=headers)
        response.raise_for_status()
        return response.json().get("metadata", {})

    def set_metadata(self, data: dict[str, str]) -> None:
        """Update room metadata.

        Parameters
        ----------
        data : dict[str, str]
            Dictionary of metadata fields to update.

        Raises
        ------
        requests.HTTPError
            If the request fails (e.g., room is locked).
        """
        headers = self._get_headers()
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/metadata", json=data, headers=headers
        )
        response.raise_for_status()

    def delete_metadata_field(self, field: str) -> None:
        """Delete a metadata field.

        Parameters
        ----------
        field : str
            The metadata field name to delete.

        Raises
        ------
        requests.HTTPError
            If the request fails (e.g., room is locked).
        """
        headers = self._get_headers()
        response = requests.delete(f"{self.url}/api/rooms/{self.room}/metadata/{field}", headers=headers)
        response.raise_for_status()

    # ========================================================================
    # Selection API Methods
    # ========================================================================

    def get_all_selections(self) -> dict:
        """Get all current selections and groups.

        Returns
        -------
        dict
            {
                "selections": {"particles": [1,2,3], "forces": [2,3]},
                "groups": {"group1": {"particles": [1,3], "forces": [1,3]}},
                "activeGroup": "group1" | None
            }
        """
        headers = self._get_headers()
        response = requests.get(f"{self.url}/api/rooms/{self.room}/selections", headers=headers)
        response.raise_for_status()
        return response.json()

    def get_selection(self, geometry: str) -> dict:
        """Get selection for a specific geometry.

        Parameters
        ----------
        geometry : str
            The geometry name (e.g., "particles", "forces")

        Returns
        -------
        dict
            {"selection": [1, 2, 3]}
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/selections/{geometry}",
            headers=headers,
        )
        response.raise_for_status()
        return response.json()

    def update_selection(self, geometry: str, indices: list[int]) -> None:
        """Update selection for a specific geometry.

        Parameters
        ----------
        geometry : str
            The geometry name (e.g., "particles", "forces")
        indices : list[int]
            List of indices to select

        Raises
        ------
        requests.HTTPError
            If the request fails (e.g., invalid indices).
        """
        headers = self._get_headers()
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/selections/{geometry}",
            json={"indices": indices},
            headers=headers,
        )
        response.raise_for_status()

    def get_selection_group(self, group_name: str) -> dict:
        """Get a specific selection group.

        Parameters
        ----------
        group_name : str
            The group name

        Returns
        -------
        dict
            {"group": {"particles": [1, 3], "forces": [1, 3]}}

        Raises
        ------
        requests.HTTPError
            If group not found (404).
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}",
            headers=headers,
        )
        response.raise_for_status()
        return response.json()

    def create_update_selection_group(
        self, group_name: str, group_data: dict[str, list[int]]
    ) -> None:
        """Create or update a selection group.

        Parameters
        ----------
        group_name : str
            The group name
        group_data : dict[str, list[int]]
            Dict mapping geometry names to index lists
            e.g., {"particles": [1, 3], "forces": [1, 3]}

        Raises
        ------
        requests.HTTPError
            If the request fails (e.g., invalid data).
        """
        headers = self._get_headers()
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}",
            json=group_data,
            headers=headers,
        )
        response.raise_for_status()

    def delete_selection_group(self, group_name: str) -> None:
        """Delete a selection group.

        Parameters
        ----------
        group_name : str
            The group name

        Raises
        ------
        requests.HTTPError
            If group not found (404).
        """
        headers = self._get_headers()
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}",
            headers=headers,
        )
        response.raise_for_status()

    def load_selection_group(self, group_name: str) -> None:
        """Load a selection group (apply it to current selections).

        This sets the active group and updates all selections to match the group.

        Parameters
        ----------
        group_name : str
            The group name

        Raises
        ------
        requests.HTTPError
            If group not found (404).
        """
        headers = self._get_headers()
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}/load",
            headers=headers,
        )
        response.raise_for_status()

    # ========================================================================
    # Bookmarks API Methods
    # ========================================================================

    def get_all_bookmarks(self) -> dict[int, str]:
        """Get all bookmarks for the room.

        Returns
        -------
        dict[int, str]
            Dictionary mapping frame indices to bookmark labels.
        """
        headers = self._get_headers()
        response = requests.get(f"{self.url}/api/rooms/{self.room}/bookmarks", headers=headers)
        response.raise_for_status()
        bookmarks = response.json().get("bookmarks", {})
        # Convert string keys from JSON back to integers
        return {int(k): v for k, v in bookmarks.items()}

    def get_bookmark(self, index: int) -> dict:
        """Get a specific bookmark by frame index.

        Parameters
        ----------
        index : int
            Frame index

        Returns
        -------
        dict
            {"index": 0, "label": "First Frame"}

        Raises
        ------
        KeyError
            If bookmark not found at the given index
        IndexError
            If index is out of range
        """
        headers = self._get_headers()
        response = requests.get(f"{self.url}/api/rooms/{self.room}/bookmarks/{index}", headers=headers)
        if response.status_code == 404:
            self._raise_for_error_type(response)

        response.raise_for_status()
        return response.json()

    def set_bookmark(self, index: int, label: str) -> None:
        """Set or update a bookmark at a specific frame index.

        Parameters
        ----------
        index : int
            Frame index
        label : str
            Bookmark label

        Raises
        ------
        ValueError
            If label is invalid (empty or not a string)
        IndexError
            If index is out of range
        """
        headers = self._get_headers()
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/bookmarks/{index}",
            json={"label": label},
            headers=headers,
        )
        if response.status_code == 400:
            self._raise_for_error_type(response)

        response.raise_for_status()

    def delete_bookmark(self, index: int) -> None:
        """Delete a bookmark at a specific frame index.

        Parameters
        ----------
        index : int
            Frame index

        Raises
        ------
        KeyError
            If bookmark not found at the given index
        """
        headers = self._get_headers()
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/bookmarks/{index}",
            headers=headers,
        )
        if response.status_code == 404:
            self._raise_for_error_type(response)

        response.raise_for_status()

    # ========================================================================
    # Screenshot API Methods
    # ========================================================================

    def list_screenshots(self, limit: int = 20, offset: int = 0) -> dict:
        """List all screenshots for the room.

        Parameters
        ----------
        limit : int
            Maximum number of screenshots to return (default 20, max 100).
        offset : int
            Number of screenshots to skip (default 0).

        Returns
        -------
        dict
            {
                "screenshots": [
                    {
                        "id": 1,
                        "format": "png",
                        "size": 123456,
                        "width": 1920,
                        "height": 1080,
                        "url": "/api/rooms/{room}/screenshots/1"
                    }
                ],
                "total": 45,
                "limit": 20,
                "offset": 0
            }
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/screenshots",
            params={"limit": limit, "offset": offset},
            headers=headers,
        )
        response.raise_for_status()
        return response.json()

    def get_screenshot_metadata(self, screenshot_id: int) -> dict:
        """Get metadata for a specific screenshot.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier

        Returns
        -------
        dict
            {
                "id": 1,
                "format": "png",
                "size": 123456,
                "width": 1920,
                "height": 1080,
                "url": "/api/rooms/{room}/screenshots/1"
            }

        Raises
        ------
        requests.HTTPError
            If screenshot not found (404).
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/screenshots/{screenshot_id}/metadata",
            headers=headers,
        )
        response.raise_for_status()
        return response.json()

    def download_screenshot(self, screenshot_id: int) -> bytes:
        """Download a screenshot as bytes.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier

        Returns
        -------
        bytes
            The image file data

        Raises
        ------
        requests.HTTPError
            If screenshot not found (404).
        """
        headers = self._get_headers()
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/screenshots/{screenshot_id}",
            headers=headers,
        )
        response.raise_for_status()
        return response.content

    def delete_screenshot(self, screenshot_id: int) -> bool:
        """Delete a screenshot.

        Parameters
        ----------
        screenshot_id : int
            Screenshot identifier

        Returns
        -------
        bool
            True if deleted successfully

        Raises
        ------
        requests.HTTPError
            If screenshot not found (404).
        """
        headers = self._get_headers()
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/screenshots/{screenshot_id}",
            headers=headers,
        )
        response.raise_for_status()
        return response.json().get("success", False)

    # ========================================================================
    # Lock API Methods
    # ========================================================================

    def lock_acquire(self, target: str, msg: str | None = None) -> dict:
        """Acquire a lock for a specific target.

        Parameters
        ----------
        target : str
            Lock target identifier (e.g., "trajectory:meta")
        msg : str | None
            Optional message describing the lock purpose

        Returns
        -------
        dict
            {"success": true, "lockToken": "...", "ttl": 60, "refreshInterval": 30}

        Raises
        ------
        RuntimeError
            If the lock is already held by another user (423 Locked)
        requests.HTTPError
            If the request fails
        """
        # Include session ID header for lock acquisition
        headers = self._get_headers()

        payload = {}
        if msg is not None:
            payload["msg"] = msg

        url = f"{self.url}/api/rooms/{self.room}/locks/{target}/acquire"
        response = requests.post(url, json=payload, headers=headers, timeout=30)

        if response.status_code == 423:  # Locked
            error_data = response.json()
            raise RuntimeError(f"Failed to acquire lock: {error_data.get('error')}")

        response.raise_for_status()
        return response.json()

    def lock_refresh(self, target: str, lock_token: str, msg: str | None = None) -> dict:
        """Refresh lock TTL and optionally update message.

        Parameters
        ----------
        target : str
            Lock target identifier (e.g., "trajectory:meta")
        lock_token : str
            Lock token from acquire response
        msg : str | None
            Optional updated message (if None, keeps existing message)

        Returns
        -------
        dict
            {"success": true}

        Raises
        ------
        PermissionError
            If not the lock holder (403 Forbidden)
        requests.HTTPError
            If the request fails
        """
        # Include session ID header for lock refresh
        headers = self._get_headers()

        payload = {"lockToken": lock_token}
        if msg is not None:
            payload["msg"] = msg

        url = f"{self.url}/api/rooms/{self.room}/locks/{target}/refresh"
        response = requests.post(url, json=payload, headers=headers, timeout=30)

        if response.status_code == 403:
            error_data = response.json()
            raise PermissionError(f"Failed to refresh lock: {error_data.get('error')}")

        response.raise_for_status()
        return response.json()

    def lock_release(self, target: str, lock_token: str) -> dict:
        """Release a lock.

        Parameters
        ----------
        target : str
            Lock target identifier (e.g., "trajectory:meta")
        lock_token : str
            Lock token from acquire response

        Returns
        -------
        dict
            {"success": true}

        Raises
        ------
        PermissionError
            If not the lock holder (403 Forbidden)
        requests.HTTPError
            If the request fails
        """
        # Include session ID header for lock release
        headers = self._get_headers()

        payload = {"lockToken": lock_token}
        url = f"{self.url}/api/rooms/{self.room}/locks/{target}/release"
        response = requests.post(url, json=payload, headers=headers, timeout=30)

        if response.status_code == 403:
            error_data = response.json()
            raise PermissionError(f"Failed to release lock: {error_data.get('error')}")

        response.raise_for_status()
        return response.json()

    def get_room_info(self) -> dict:
        """Get room information including frameCount.

        Returns
        -------
        dict
            {
                "id": str,
                "description": str | None,
                "frameCount": int,
                "locked": bool,
                "hidden": bool,
                "metadata": dict
            }

        Raises
        ------
        requests.HTTPError
            If the request fails (e.g., room not found)
        """
        headers = self._get_headers()
        url = f"{self.url}/api/rooms/{self.room}"
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        return response.json()

    def get_frame_selection(self) -> list[int] | None:
        """Get frame selection for the room.

        Returns
        -------
        list[int] | None
            List of selected frame indices, or None if no selection

        Raises
        ------
        requests.HTTPError
            If the request fails
        """
        headers = self._get_headers()
        url = f"{self.url}/api/rooms/{self.room}/frame-selection"
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        data = response.json()
        return data.get("frameSelection")

    def get_step(self) -> dict:
        """Get current step/frame for the room.

        Returns
        -------
        dict
            {"step": int, "totalFrames": int}

        Raises
        ------
        RuntimeError
            If the GET request fails
        """
        headers = self._get_headers()
        url = f"{self.url}/api/rooms/{self.room}/step"
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
        return response.json()

    def update_step(self, step: int) -> dict:
        """Set current step/frame for the room.

        Requires holding the 'step' lock for the room.

        Parameters
        ----------
        step : int
            Frame index to set (non-negative integer)

        Returns
        -------
        dict
            {"success": bool, "step": int}

        Raises
        ------
        ValueError
            If step is negative or non-integer
        RuntimeError
            If the PUT request fails (401, 403, 423)
        """
        if not isinstance(step, int) or step < 0:
            raise ValueError("Step must be a non-negative integer")

        headers = self._get_headers()
        url = f"{self.url}/api/rooms/{self.room}/step"
        response = requests.put(url, json={"step": step}, headers=headers, timeout=10)

        if response.status_code == 400:
            error_data = response.json()
            raise ValueError(f"Invalid step value: {error_data.get('error')}")
        elif response.status_code == 401:
            error_data = response.json()
            raise PermissionError(f"Authentication failed: {error_data.get('error')}")
        elif response.status_code == 403:
            error_data = response.json()
            raise PermissionError(f"Access denied: {error_data.get('error')}")
        elif response.status_code == 423:
            from zndraw.exceptions import LockError
            error_data = response.json()
            raise LockError(f"Lock not held for target='step': {error_data.get('error')}")

        response.raise_for_status()
        return response.json()

    def progress_start(self, progress_id: str, description: str):
        """Start tracking progress for a long-running operation via REST API.

        Parameters
        ----------
        progress_id : str
            Unique progress identifier
        description : str
            Progress description

        Returns
        -------
        dict
            Response from server
        """
        url = f"{self.url}/api/rooms/{self.room}/progress"
        response = requests.post(
            url,
            json={"progressId": progress_id, "description": description},
            headers=self._get_headers(),
        )
        response.raise_for_status()
        return response.json()

    def progress_update(
        self, progress_id: str, description: str | None = None, progress: float | None = None
    ):
        """Update progress tracking for an ongoing operation via REST API.

        Parameters
        ----------
        progress_id : str
            Progress identifier
        description : str | None
            Updated progress description (optional)
        progress : float | None
            Progress percentage 0-100 (optional)

        Returns
        -------
        dict
            Response from server
        """
        data = {}
        if description is not None:
            data["description"] = description
        if progress is not None:
            data["progress"] = progress

        url = f"{self.url}/api/rooms/{self.room}/progress/{progress_id}"
        response = requests.put(url, json=data, headers=self._get_headers())
        response.raise_for_status()
        return response.json()

    def progress_complete(self, progress_id: str):
        """Complete and remove progress tracking for an operation via REST API.

        Parameters
        ----------
        progress_id : str
            Progress identifier

        Returns
        -------
        dict
            Response from server
        """
        url = f"{self.url}/api/rooms/{self.room}/progress/{progress_id}"
        response = requests.delete(url, headers=self._get_headers())
        response.raise_for_status()
        return response.json()
