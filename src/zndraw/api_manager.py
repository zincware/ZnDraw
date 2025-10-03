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

    def join_room(self, template: t.Any, user_id: str = "guest") -> dict:
        # Import here to avoid circular import
        from zndraw.zndraw import _TemplateValue

        # If template is _TemplateValue (default), don't send it (use server default)
        # If template is None, send {"template": None} (create blank room)
        # If template is a string, send {"template": template_name}
        payload = {"userId": user_id}
        if template is not _TemplateValue:
            payload["template"] = template
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
