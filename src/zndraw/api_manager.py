import dataclasses
import typing as t

import msgpack
import requests

from zndraw.storage import decode_data, encode_data


@dataclasses.dataclass
class APIManager:
    url: str
    room: str
    client_id: str | None = None

    def join_room(
        self,
        description: str | None = None,
        copy_from: str | None = None,
        user_id: str = "guest"
    ) -> dict:
        """Join a room, optionally creating it with a description or copying from an existing room.
        
        Args:
            description: Optional description for the room (only used if room is created)
            copy_from: Optional room ID to copy frames and settings from (only used if room is created)
            user_id: User identifier
        
        Returns:
            Dict containing room information and join token
        """
        payload = {"userId": user_id}
        
        if description is not None:
            payload["description"] = description
        if copy_from is not None:
            payload["copyFrom"] = copy_from
        if self.client_id:
            payload["clientId"] = self.client_id

        response = requests.post(f"{self.url}/api/rooms/{self.room}/join", json=payload)

        if response.status_code != 200:
            raise RuntimeError(
                f"Failed to join room '{self.room}': {response.status_code} {response.text}"
            )
        return response.json()

    def get_next_job(self, worker_id: str) -> dict | None:
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/jobs/next", json={"workerId": worker_id}
        )
        if response.status_code == 200:
            return response.json()
        elif response.status_code == 400:
            print("No jobs available.")
            return None
        else:
            response.raise_for_status()
            return None

    def update_job_status(
        self, job_id: str, status: str, worker_id: str, error: str | None = None
    ) -> None:
        payload = {"status": status, "workerId": worker_id}
        if error:
            payload["error"] = error
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/jobs/{job_id}/status",
            json=payload,
        )
        response.raise_for_status()

    def register_extension(
        self, name: str, category: str, schema: dict, client_id: str
    ) -> None:
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/register",
            json={
                "name": name,
                "category": category,
                "schema": schema,
                "clientId": client_id,
            },
        )
        response.raise_for_status()

    def get_frames(self, indices_or_slice, keys: list[str] | None = None) -> list[dict]:
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

        full_url = f"{self.url}/api/rooms/{self.room}/frames"
        response = requests.get(full_url, params=payload, timeout=30)

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

        serialized_frames = msgpack.unpackb(response.content, strict_map_key=False)
        return [decode_data(frame) for frame in serialized_frames]

    def upload_frames(self, action: str, data, **kwargs) -> dict:
        try:
            serialized_data = (
                encode_data(data)
                if isinstance(data, dict)
                else [encode_data(frame) for frame in data]
            )
            packed_data = msgpack.packb(serialized_data)

            upload_url = f"{self.url}/api/rooms/{self.room}/frames"
            params = {"action": action}
            params.update(kwargs)
            
            # Add client_id if available (for lock checking)
            if self.client_id:
                params["client_id"] = self.client_id

            http_response = requests.post(
                upload_url, data=packed_data, params=params, timeout=30
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
        
        # Add client_id if available (for lock checking)
        if self.client_id:
            params["client_id"] = self.client_id

        response = requests.delete(delete_url, params=params, timeout=30)

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
        serialized_data = [encode_data(value) for value in data]
        packed_data = msgpack.packb(serialized_data)

        bulk_url = f"{self.url}/api/rooms/{self.room}/frames/bulk"
        params = {}
        if start is not None and stop is not None:
            params = {"start": start, "stop": stop}
        elif indices is not None:
            params = {"indices": ",".join(str(i) for i in indices)}
        else:
            raise ValueError("Either start/stop or indices must be provided.")

        response = requests.patch(bulk_url, data=packed_data, params=params, timeout=30)

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

    def get_extension_settings(self, extension_name: str, user_id: str) -> dict:
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/extensions/settings/{extension_name}/data?userId={user_id}"
        )
        response.raise_for_status()
        return response.json()

    def submit_extension_settings(
        self, extension_name: str, user_id: str, data: dict
    ) -> None:
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/settings/{extension_name}/submit?userId={user_id}",
            json=data,
        )
        response.raise_for_status()

    def run_extension(self, category: str, name: str, user_id: str, data: dict) -> dict:
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/extensions/{category}/{name}/submit",
            json={"data": data, "userId": user_id},
        )
        if response.status_code != 200:
            error_data = response.json()
            error_type = error_data.get("type", "")
            error_msg = error_data.get("error", response.text)

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

        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/chat/messages", params=params
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
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/geometries",
            json={"key": key, "data": data, "type": geometry_type},
        )
        response.raise_for_status()

    def get_geometry(self, key: str) -> dict | None:
        """Get a specific geometry by key.
        
        Args:
            key: The geometry key/name
            
        Returns:
            The geometry dict with 'type' and 'data' keys, or None if not found
        """
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/geometries/{key}",
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
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/geometries/{key}",
        )
        response_data = response.json()
        if response.status_code != 200:
            error = response_data.get("error", None)
            error_type = response_data.get("type", None)
            if error_type == "KeyError":
                raise KeyError(error)
        response.raise_for_status()

    def list_geometries(self) -> list[str]:
        """List all geometry keys.
        
        Returns:
            List of geometry keys
        """
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/geometries",
        )
        response.raise_for_status()
        return response.json().get("geometries", [])
    
    def get_geometries(self) -> dict | None:
        """Get all geometries (DEPRECATED: use list_geometries and get_geometry).
        
        This method fetches all geometry keys and then retrieves each geometry.
        For better performance, use list_geometries() and get_geometry(key) as needed.
        
        Returns:
            Dict mapping geometry keys to their data {key: {type, data}, ...}
        """
        keys = self.list_geometries()
        if not keys:
            return {}
        
        result = {}
        for key in keys:
            geometry = self.get_geometry(key)
            if geometry:
                result[key] = geometry
        return result

    def add_figure(self, key: str, figure: dict) -> None:
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/figures",
            json={"key": key, "figure": figure},
        )
        response.raise_for_status()

    def delete_figure(self, key: str) -> None:
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/figures/{key}",
        )
        response_data = response.json()
        if response.status_code != 200:
            error = response_data.get("error", None)
            error_type = response_data.get("type", None)
            if error_type == "KeyError":
                raise KeyError(error)

        response.raise_for_status()

    def get_figure(self, key: str) -> dict | None:
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/figures/{key}",
        )
        if response.status_code == 404:
            return None
        response.raise_for_status()
        return response.json().get("figure", None)

    def list_figures(self) -> list[str]:
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/figures",
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
        response = requests.get(f"{self.url}/api/rooms/{self.room}/metadata")
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
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/metadata",
            json=data
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
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/metadata/{field}"
        )
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
        response = requests.get(f"{self.url}/api/rooms/{self.room}/selections")
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
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/selections/{geometry}"
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
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/selections/{geometry}",
            json={"indices": indices},
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
        response = requests.get(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}"
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
        response = requests.put(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}",
            json=group_data,
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
        response = requests.delete(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}"
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
        response = requests.post(
            f"{self.url}/api/rooms/{self.room}/selections/groups/{group_name}/load"
        )
        response.raise_for_status()
