import logging
import re
import traceback
import typing as t

import socketio

from zndraw.app.constants import SocketEvents

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

log = logging.getLogger(__name__)


class SocketManager:
    def __init__(self, zndraw_instance: "ZnDraw"):
        self.zndraw = zndraw_instance
        self.sio = socketio.Client()
        self._initial_connect_done = False
        self._register_handlers()

    def _register_handlers(self):
        self.sio.on("connect", self._on_connect)
        self.sio.on(SocketEvents.FRAME_UPDATE, self._on_frame_update)
        self.sio.on(SocketEvents.SELECTION_UPDATE, self._on_selection_update)
        self.sio.on(SocketEvents.ROOM_UPDATE, self._on_room_update)
        self.sio.on(SocketEvents.JOB_ASSIGNED, self._on_job_assigned)
        self.sio.on(
            SocketEvents.FRAME_SELECTION_UPDATE, self._on_frame_selection_update
        )
        self.sio.on(SocketEvents.INVALIDATE_BOOKMARK, self._on_bookmarks_invalidate)
        self.sio.on(SocketEvents.INVALIDATE_FRAMES, self._on_frames_invalidate)
        self.sio.on(SocketEvents.INVALIDATE_GEOMETRY, self._on_geometry_invalidate)
        self.sio.on(SocketEvents.INVALIDATE_FIGURE, self._on_figure_invalidate)
        self.sio.on(SocketEvents.FILESYSTEM_LIST, self._on_filesystem_list)
        self.sio.on(SocketEvents.FILESYSTEM_METADATA, self._on_filesystem_metadata)
        self.sio.on(SocketEvents.FILESYSTEM_LOAD, self._on_filesystem_load)

    def connect(self):
        """Connect to server with JWT authentication and join room."""
        if self.sio.connected:
            log.debug("Already connected.")
            return
        # Connect with JWT token only (sessionId is created by room:join)
        self.sio.connect(
            self.zndraw.url,
            auth={"token": self.zndraw.api.jwt_token},
            wait=True,
        )

        # Join room and get sessionId + room data
        room_id = self.zndraw.room
        response = self.sio.call(
            "room:join", {"roomId": room_id, "clientType": "python"}
        )

        # Handle room not found - create it first
        if response.get("code") == 404:
            self.zndraw.api.create_room(
                description=self.zndraw.description,
                copy_from=self.zndraw.copy_from,
                template="none",  # Python client will upload its own data
            )
            # Retry join after creation
            response = self.sio.call(
                "room:join", {"roomId": room_id, "clientType": "python"}
            )

        if response.get("status") != "ok":
            raise RuntimeError(f"Failed to join room: {response.get('message')}")

        # Store sessionId for subsequent operations
        self.zndraw.api.session_id = response["sessionId"]

        # Initialize minimal room data from socket response
        self.zndraw._len = response.get("frameCount", 0)
        self.zndraw._step = response.get("step", 0)

        # Fetch additional data via REST
        # Geometries (includes schemas/defaults)
        geometries = self.zndraw.api.get_geometries()
        if geometries is not None:
            self.zndraw._geometries = geometries

        # Frame selection
        frame_selection = self.zndraw.api.get_frame_selection()
        self.zndraw._frame_selection = frozenset(frame_selection or [])

        # Bookmarks
        bookmarks = self.zndraw.api.get_all_bookmarks()
        self.zndraw._bookmarks = bookmarks

        # Re-register extensions after session is established
        self._register_extensions_after_join()

        # Mark initial connection as done - subsequent connects are reconnects
        self._initial_connect_done = True

    def disconnect(self):
        if self.sio.connected:
            self.sio.disconnect()
            # Reset flag so next connect() goes through full flow
            self._initial_connect_done = False
            log.debug("Disconnected.")

    @property
    def connected(self) -> bool:
        return self.sio.connected

    def _on_connect(self):
        """Handle socket connection event.

        For initial connection, connect() handles room:join and extension
        registration. For auto-reconnects (e.g., after network issues),
        we need to re-join the room and re-register extensions here.
        """
        log.debug("Connected to server")

        # On initial connect, connect() handles everything
        if not self._initial_connect_done:
            return

        # This is an auto-reconnect - need to re-join room and re-register
        log.debug(
            "Auto-reconnect detected, re-joining room and re-registering extensions"
        )

        # Re-join room to get new session
        room_id = self.zndraw.room
        response = self.sio.call(
            "room:join", {"roomId": room_id, "clientType": "python"}
        )

        if response.get("status") != "ok":
            log.error(
                f"Failed to re-join room after reconnect: {response.get('message')}"
            )
            return

        # Update session ID and minimal state
        self.zndraw.api.session_id = response["sessionId"]
        self.zndraw._len = response.get("frameCount", 0)
        self.zndraw._step = response.get("step", 0)

        # Re-register extensions
        self._register_extensions_after_join()

    def _register_extensions_after_join(self):
        """Re-register extensions after session is established.

        Called from connect() after room:join succeeds and session_id is set.
        This ensures the session exists before attempting to register extensions.
        """
        # Re-register any extensions that were registered before connection
        # Process public extensions
        for name, ext in self.zndraw._public_extensions.items():
            worker_id = self.zndraw.api.register_extension(
                name=name,
                category=ext["extension"].category,
                schema=ext["extension"].model_json_schema(),
                socket_manager=self,
                public=True,
            )
            # Store the worker_id assigned by server
            if worker_id:
                self.zndraw._worker_id = worker_id

        # Process private extensions
        for name, ext in self.zndraw._private_extensions.items():
            worker_id = self.zndraw.api.register_extension(
                name=name,
                category=ext["extension"].category,
                schema=ext["extension"].model_json_schema(),
                socket_manager=self,
                public=False,
            )
            # Store the worker_id assigned by server
            if worker_id:
                self.zndraw._worker_id = worker_id

        # Re-register any filesystems that were registered before connection
        for name, fs in self.zndraw._filesystems.items():
            provider = fs["provider"]
            worker_id = self.zndraw.api.register_filesystem(
                name=name,
                provider_type=provider.__class__.__name__,
                root_path=provider.root_path,
                socket_manager=self,
                public=fs["public"],
            )
            # Store the worker_id assigned by server
            if worker_id:
                self.zndraw._worker_id = worker_id

    def _on_frame_update(self, data):
        if "frame" in data:
            self.zndraw._step = data["frame"]

    def _on_room_update(self, data):
        """Handle room:update events (consolidated room metadata updates).

        Note: We invalidate the cache instead of setting _len directly to avoid
        race conditions. Socket events are asynchronous and can arrive after
        a more recent HTTP response has already set the correct _len value.
        Invalidating forces the next len() call to fetch fresh data.
        """
        if "frameCount" in data:
            self.zndraw._len = None

    def _on_selection_update(self, data):
        if "indices" in data:
            self.zndraw._selection = frozenset(data["indices"])

    def _on_frame_selection_update(self, data):
        if "indices" in data:
            self.zndraw._frame_selection = frozenset(data["indices"])

    def _on_bookmarks_invalidate(self, data):
        """Handle bookmark invalidation by updating specific entry or refetching."""
        index = data.get("index")
        operation = data.get("operation")

        if operation == "set" and index is not None:
            # Targeted update - only fetch and update the specific bookmark
            try:
                response = self.zndraw.api.get_bookmark(index)
                label = response.get("label")
                if label is not None:
                    self.zndraw._bookmarks[index] = label
            except Exception as e:
                # Check for 404 - bookmark was deleted
                response = getattr(e, "response", None)
                if response is not None and response.status_code == 404:
                    self.zndraw._bookmarks.pop(index, None)
                else:
                    log.error("Failed to fetch bookmark %s: %s", index, e)
                    # Fallback to full refresh
                    self._refresh_all_bookmarks()
        elif operation == "delete" and index is not None:
            # Remove the specific bookmark from cache
            self.zndraw._bookmarks.pop(index, None)
        else:
            # Full refresh for bulk operations (clear, shift, etc.)
            self._refresh_all_bookmarks()

    def _refresh_all_bookmarks(self):
        """Refresh all bookmarks from server, retaining cache on failure."""
        try:
            bookmarks = self.zndraw.api.get_all_bookmarks()
            self.zndraw._bookmarks = bookmarks
        except Exception as e:
            log.error("Failed to refresh bookmarks: %s", e)

    def _on_geometry_invalidate(self, data):
        if data.get("key"):
            # Refresh geometries from server to get updated state
            response = self.zndraw.api.get_geometries()
            if response is not None:
                self.zndraw._geometries = response

    def _on_figure_invalidate(self, data):
        if key := data.get("key"):
            self.zndraw._figures.pop(key, None)

    def _on_job_assigned(self, data: dict):
        """Handle job:assign socket event.

        New workflow (REST-based):
        1. Receive socket notification with jobId
        2. Fetch full job details via GET /api/jobs/{jobId}
        3. Execute job using shared executor
        4. Status updates handled by executor via PUT endpoints

        Parameters
        ----------
        data : dict
            Event payload containing:
            - jobId: Job identifier to fetch and execute
        """
        job_id = data.get("jobId")

        if not job_id:
            log.error("Received job:assign event without jobId")
            return

        if not self.zndraw.auto_pickup_jobs:
            log.debug(f"Auto-pickup disabled, ignoring job assignment {job_id}")
            return

        log.debug(f"Worker {self.zndraw.sid} received job assignment: {job_id}")

        # Use shared job executor (same code as Celery workers)
        from zndraw.job_executor import execute_job_for_worker

        # Build extension registry from parent's public and private extensions
        # Store entire _ExtensionStore objects to preserve run_kwargs
        extension_registry = {}
        for ext_name, ext_store in self.zndraw._public_extensions.items():
            # Use actual category from extension, not all categories
            category = ext_store["extension"].category.value
            extension_registry[f"{category}:{ext_name}"] = ext_store
        for ext_name, ext_store in self.zndraw._private_extensions.items():
            category = ext_store["extension"].category.value
            extension_registry[f"{category}:{ext_name}"] = ext_store

        try:
            execute_job_for_worker(
                job_id=job_id,
                server_url=self.zndraw.url,
                worker_id=self.zndraw.sid,
                extension_registry=extension_registry,
            )
        except Exception as e:
            log.error(f"Error executing assigned job {job_id}: {e}", exc_info=True)
            # Error handling is done inside execute_job_for_worker

    def _on_task_run(
        self,
        data: dict,
        extension: str,
        category: str,
        room: str,
        public: str = "false",
    ):
        """Execute an extension job with proper room context.

        Creates a temporary ZnDraw instance connected to the job's target room,
        ensuring the extension operates in the correct room context.

        Parameters
        ----------
        data : dict
            Job input parameters
        extension : str
            Extension name
        category : str
            Extension category
        room : str
            Room where the job was triggered (may differ from worker's room)
        public : str
            Whether this is a public/global extension ("true" or "false")
        """
        from zndraw import ZnDraw

        # Create temporary ZnDraw instance for the job's target room
        temp_vis = ZnDraw.for_job_execution(
            url=self.zndraw.url,
            room=room,
            user=self.zndraw.user,
            password=getattr(self.zndraw, "password", None),
        )

        try:
            # Redis stores booleans as strings, so check for "true"
            public_bool = public == "true"

            # Find the extension in the correct namespace
            extensions_dict = (
                self.zndraw._public_extensions
                if public_bool
                else self.zndraw._private_extensions
            )
            ext_data = extensions_dict.get(extension)

            if ext_data is None:
                namespace = "public" if public_bool else "private"
                raise ValueError(
                    f"Extension '{extension}' not found in {namespace} namespace. "
                    f"Available extensions in {namespace}: {list(extensions_dict.keys())}"
                )

            ext = ext_data["extension"]
            instance = ext(**(data))
            instance.run(temp_vis, **(ext_data["run_kwargs"] or {}))
        finally:
            # Cleanup: disconnect temporary instance
            if hasattr(temp_vis, "socket") and temp_vis.socket:
                temp_vis.socket.disconnect()

    def _on_frames_invalidate(self, data: dict):
        log.debug(f"Received cache invalidation event: {data}")
        if self.zndraw.cache is None:
            return

        operation = data.get("operation")
        affected_keys = data.get("affectedKeys")

        # Log key-level invalidation info if present
        if affected_keys:
            log.warning(
                "Key-level invalidation received, but frame cache does not support key-level invalidation. Clearing entire frame cache."
            )

        if operation == "replace":
            idx = data.get("affectedIndex")
            if idx is not None:
                self.zndraw.cache.pop(idx, None)
        elif operation in ("insert", "delete", "bulk_replace"):
            from_idx = data.get("affectedFrom")
            if from_idx is not None:
                self.zndraw.cache.invalidate_from(from_idx)
        elif operation == "clear_all":
            self.zndraw.cache.clear()
        else:
            log.warning(
                "Unknown or broad invalidation event received. Clearing entire frame cache."
            )
            self.zndraw.cache.clear()

    def _on_filesystem_list(self, data: dict):
        """Handle filesystem list request from server.

        Server sends: {requestId, fsName, public, path, recursive, filterExtensions, search}
        Client responds: {requestId, success, files, error}
        """
        request_id = data.get("requestId")
        fs_name = data.get("fsName")
        public = data.get("public", False)
        path = data.get("path", "")
        recursive = data.get("recursive", False)
        filter_extensions = data.get("filterExtensions")
        search = data.get("search")

        try:
            # Get the filesystem instance using composite key
            fs_key = (
                f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
            )
            if fs_key not in self.zndraw._filesystems:
                available_keys = list(self.zndraw._filesystems.keys())
                raise ValueError(
                    f"Filesystem '{fs_name}' (key: '{fs_key}') not found. "
                    f"Available: {available_keys}. "
                    f"Public flag: {public}, Room: {self.zndraw.room}"
                )

            fs = self.zndraw._filesystems[fs_key]["fs"]

            # Compile search pattern if provided
            search_pattern = None
            if search:
                try:
                    search_pattern = re.compile(search, re.IGNORECASE)
                except re.error as e:
                    raise ValueError(f"Invalid search pattern: {e}")

            # List files using fsspec API
            if recursive:
                # Use glob for recursive listing
                pattern = f"{path}/**" if path else "**"
                all_paths = fs.glob(pattern, detail=True)
            else:
                # Use ls for non-recursive listing
                all_paths = fs.ls(path or ".", detail=True)

            # Convert fsspec info dicts to file metadata format
            files_data = []
            for item in all_paths:
                # Handle both dict and tuple returns from fs.ls()
                if isinstance(item, dict):
                    info = item
                    item_path = info.get("name", "")
                else:
                    # Fallback for simple path strings
                    info = fs.info(item)
                    item_path = item

                item_name = item_path.split("/")[-1]

                # Filter by extension if specified
                if filter_extensions:
                    if not any(item_path.endswith(ext) for ext in filter_extensions):
                        continue

                # Filter by search pattern if specified
                if search_pattern:
                    if not search_pattern.search(item_name):
                        continue

                # Build file metadata
                file_data = {
                    "name": item_name,
                    "path": item_path,
                    "size": info.get("size", 0),
                    "type": info.get("type", "file"),
                    "modified": info.get("mtime"),
                }
                files_data.append(file_data)

            # Return response (for socketio.call() on server)
            return {
                "requestId": request_id,
                "success": True,
                "files": files_data,
            }

        except Exception as e:
            log.error(f"Error listing files from filesystem '{fs_name}': {e}")
            traceback.print_exc()
            return {
                "requestId": request_id,
                "success": False,
                "error": str(e),
            }

    def _on_filesystem_load(self, data: dict):
        """Handle filesystem load request from server.

        Server sends: {
            requestId, fsName, public, path, room,
            batchSize (optional), start (optional), stop (optional), step (optional)
        }
        Client: Reads file, uploads to room, sends completion response
        """
        import itertools
        from pathlib import Path

        import ase.db
        import ase.io
        import h5py
        import znh5md

        from zndraw.app.file_browser import FORMAT_BACKENDS
        from zndraw.app.tasks import apply_adaptive_resolution, upload_frames_to_room

        request_id = data.get("requestId")
        fs_name = data.get("fsName")
        public = data.get("public", False)
        path = data.get("path")
        target_room = data.get("room")
        batch_size = data.get("batchSize", 10)
        start = data.get("start")
        stop = data.get("stop")
        step = data.get("step")

        temp_vis = None
        try:
            # Get the filesystem instance
            fs_key = (
                f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
            )
            if fs_key not in self.zndraw._filesystems:
                raise ValueError(f"Filesystem '{fs_name}' not found")

            fs = self.zndraw._filesystems[fs_key]["fs"]

            from zndraw import ZnDraw

            temp_vis = ZnDraw(
                url=self.zndraw.url,
                room=target_room,
                user=self.zndraw.user,
                password=getattr(self.zndraw, "password", None),
            )

            temp_vis.log(f"Loading {path} from {fs_name}...")

            # Determine backend and file mode
            ext = Path(path).suffix.lstrip(".").lower()
            backends = FORMAT_BACKENDS.get(ext, ["ASE"])
            backend_name = backends[0]
            binary_formats = {"h5", "h5md", "traj", "db", "nc", "hdf5"}
            mode = "rb" if ext in binary_formats else "r"

            with fs.open(path, mode) as f:
                frame_iterator = None
                total_expected_frames = None

                if backend_name == "ZnH5MD":
                    if step is not None and step <= 0:
                        raise ValueError("Step must be a positive integer")
                    h5_file = h5py.File(f, "r")
                    io = znh5md.IO(file_handle=h5_file)
                    n_frames = len(io)

                    if start is not None and start >= n_frames:
                        raise ValueError(f"Start frame {start} exceeds file length")

                    s = slice(start, stop, step)
                    _start, _stop, _step = s.indices(n_frames)
                    total_expected_frames = (_stop - _start + _step - 1) // _step
                    frame_iterator = itertools.islice(io, _start, _stop, _step)

                elif backend_name == "ASE-DB":
                    temp_vis.log("⚠️ ASE database may not support remote filesystems")
                    db = ase.db.connect(f)
                    n_rows = db.count()

                    if n_rows == 0:
                        raise ValueError("Database is empty")
                    if step is not None and step <= 0:
                        raise ValueError("Step must be a positive integer")

                    _start = start if start is not None else 0
                    _stop = min(stop if stop is not None else n_rows, n_rows)
                    _step = step if step is not None else 1

                    if _start >= n_rows:
                        raise ValueError(f"Start row {_start} exceeds database size")

                    selected_ids = list(range(_start + 1, _stop + 1, _step))
                    total_expected_frames = len(selected_ids)

                    def db_iterator():
                        for row_id in selected_ids:
                            try:
                                yield db.get(id=row_id).toatoms()
                            except KeyError:
                                continue

                    frame_iterator = db_iterator()

                else:
                    start_str = str(start) if start is not None else ""
                    stop_str = str(stop) if stop is not None else ""
                    step_str = str(step) if step is not None else ""
                    index_str = f"{start_str}:{stop_str}:{step_str}"
                    frame_iterator = ase.io.iread(
                        f, index=index_str, format=ext or None
                    )

                # Upload frames using shared helper
                if frame_iterator:
                    loaded_frame_count, max_particles = upload_frames_to_room(
                        temp_vis,
                        frame_iterator,
                        batch_size=batch_size,
                        total_expected_frames=total_expected_frames,
                        progress_description=f"Loading {path}",
                    )
                else:
                    loaded_frame_count, max_particles = 0, 0

            temp_vis.log(
                f"✓ Successfully loaded {loaded_frame_count} frames from {path}"
            )

            if loaded_frame_count > 0:
                apply_adaptive_resolution(temp_vis, max_particles)

            return {
                "requestId": request_id,
                "success": True,
                "frameCount": loaded_frame_count,
            }

        except Exception as e:
            log.error(f"Error loading file from filesystem '{fs_name}': {e}")
            traceback.print_exc()
            return {
                "requestId": request_id,
                "success": False,
                "error": str(e),
            }
        finally:
            if temp_vis is not None:
                try:
                    temp_vis.disconnect()
                except Exception as e:
                    log.error(f"Error disconnecting temp_vis: {e}")

    def _on_filesystem_metadata(self, data: dict):
        """Handle filesystem metadata request from server.

        Server sends: {requestId, fsName, public, path}
        Client responds: {requestId, success, metadata, error}
        """
        request_id = data.get("requestId")
        fs_name = data.get("fsName")
        public = data.get("public", False)
        path = data.get("path")

        try:
            # Get the filesystem instance using composite key
            fs_key = (
                f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
            )
            if fs_key not in self.zndraw._filesystems:
                available_keys = list(self.zndraw._filesystems.keys())
                raise ValueError(
                    f"Filesystem '{fs_name}' (key: '{fs_key}') not found. "
                    f"Available: {available_keys}. "
                    f"Public flag: {public}, Room: {self.zndraw.room}"
                )

            fs = self.zndraw._filesystems[fs_key]["fs"]

            # Get metadata using fsspec API
            info = fs.info(path)

            # Convert to standard metadata format
            metadata = {
                "name": path.split("/")[-1],
                "path": path,
                "size": info.get("size", 0),
                "type": info.get("type", "file"),
                "modified": info.get("mtime"),
                "created": info.get("created"),
            }

            # Return response (for socketio.call() on server)
            return {
                "requestId": request_id,
                "success": True,
                "metadata": metadata,
            }

        except Exception as e:
            log.error(f"Error getting metadata from filesystem '{fs_name}': {e}")
            traceback.print_exc()
            return {
                "requestId": request_id,
                "success": False,
                "error": str(e),
            }
