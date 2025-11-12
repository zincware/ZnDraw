import dataclasses
import logging
import threading
import time
import traceback
import typing as t
import warnings

import socketio

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

log = logging.getLogger(__name__)



class SocketManager:
    def __init__(self, zndraw_instance: "ZnDraw"):
        self.zndraw = zndraw_instance
        self.sio = socketio.Client()
        self._register_handlers()

    def _register_handlers(self):
        self.sio.on("connect", self._on_connect)
        self.sio.on("frame_update", self._on_frame_update)
        self.sio.on("selection:update", self._on_selection_update)
        self.sio.on("room:update", self._on_room_update)
        self.sio.on("invalidate", self._on_invalidate)
        self.sio.on("job:assigned", self._on_job_assigned)
        self.sio.on("frame_selection:update", self._on_frame_selection_update)
        self.sio.on("bookmarks:invalidate", self._on_bookmarks_invalidate)
        self.sio.on("frames:invalidate", self._on_frames_invalidate)
        self.sio.on("invalidate:geometry", self._on_geometry_invalidate)
        self.sio.on("invalidate:figure", self._on_figure_invalidate)
        self.sio.on("filesystem:list", self._on_filesystem_list)
        self.sio.on("filesystem:metadata", self._on_filesystem_metadata)
        self.sio.on("filesystem:load", self._on_filesystem_load)

    def connect(self):
        """Connect to server with JWT authentication."""
        if self.sio.connected:
            print("Already connected.")
            return
        # Connect with JWT token and sessionId for authentication
        self.sio.connect(
            self.zndraw.url,
            auth={
                "token": self.zndraw.api.jwt_token,
                "sessionId": self.zndraw.api.session_id,
            },
            wait=True,
        )

    def disconnect(self):
        if self.sio.connected:
            self.sio.disconnect()
            print("Disconnected.")

    @property
    def connected(self) -> bool:
        return self.sio.connected

    def _on_connect(self):
        """Handle connection to server."""
        log.debug(f"Connected to server")

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
        """Handle room:update events (consolidated room metadata updates)."""
        if "frameCount" in data:
            self.zndraw._len = data["frameCount"]

    def _on_selection_update(self, data):
        if "indices" in data:
            self.zndraw._selection = frozenset(data["indices"])

    def _on_frame_selection_update(self, data):
        if "indices" in data:
            self.zndraw._frame_selection = frozenset(data["indices"])

    def _on_bookmarks_invalidate(self, data):
        """Handle bookmark invalidation by refetching from server."""
        # Refetch all bookmarks from server to update local cache
        bookmarks = self.zndraw.api.get_all_bookmarks()
        self.zndraw._bookmarks = bookmarks

    def _on_geometry_invalidate(self, data):
        if key := data.get("key"):
            # Refresh geometries from server to get updated state
            response = self.zndraw.api.get_geometries()
            if response is not None:
                self.zndraw._geometries = response

    def _on_figure_invalidate(self, data):
        if key := data.get("key"):
            self.zndraw._figures.pop(key, None)

    def _on_job_assigned(self, data: dict):
        """Handle job:assigned socket event.

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
            log.error("Received job:assigned event without jobId")
            return

        if not self.zndraw.auto_pickup_jobs:
            log.info(f"Auto-pickup disabled, ignoring job assignment {job_id}")
            return

        log.info(f"Worker {self.zndraw.sid} received job assignment: {job_id}")

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

    def _on_task_run(self, data: dict, extension: str, category: str, room: str, public: str = "false"):
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
            extensions_dict = self.zndraw._public_extensions if public_bool else self.zndraw._private_extensions
            ext_data = extensions_dict.get(extension)

            if ext_data is None:
                namespace = "public" if public_bool else "private"
                raise ValueError(
                    f"Extension '{extension}' not found in {namespace} namespace. "
                    f"Available extensions in {namespace}: {list(extensions_dict.keys())}"
                )

            ext = ext_data["extension"]
            instance = ext(**(data))
            instance.run(
                temp_vis, **(ext_data["run_kwargs"] or {})
            )
        finally:
            # Cleanup: disconnect temporary instance
            if hasattr(temp_vis, "socket") and temp_vis.socket:
                temp_vis.socket.disconnect()

    def _on_invalidate(self, data: dict):
        if data["category"] == "settings":
            self.zndraw._settings.pop(data["extension"], None)

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
        import re

        request_id = data.get("requestId")
        fs_name = data.get("fsName")
        public = data.get("public", False)
        path = data.get("path", "")
        recursive = data.get("recursive", False)
        filter_extensions = data.get("filterExtensions")
        search = data.get("search")

        try:
            # Get the filesystem instance using composite key
            fs_key = f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
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

        This triggers the client to read a file from the filesystem and upload
        it to a target room using normal vis.extend() logic.

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
            # Get the filesystem instance using composite key
            fs_key = f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
            if fs_key not in self.zndraw._filesystems:
                available_keys = list(self.zndraw._filesystems.keys())
                raise ValueError(
                    f"Filesystem '{fs_name}' (key: '{fs_key}') not found. "
                    f"Available: {available_keys}. "
                    f"Public flag: {public}, Room: {self.zndraw.room}"
                )

            fs = self.zndraw._filesystems[fs_key]["fs"]

            # Create a temporary ZnDraw instance for the target room
            from zndraw import ZnDraw

            temp_vis = ZnDraw(
                url=self.zndraw.url,
                room=target_room,
                user=self.zndraw.user,
                password=getattr(self.zndraw, "password", None),
            )

            # Use progress tracking for loading progress
            with temp_vis.progress_tracker(f"Loading {path} from {fs_name}") as task:
                temp_vis.log(f"Loading {path} from {fs_name}...")

                # Import utilities from tasks.py
                from zndraw.app.tasks import (
                    FORMAT_BACKENDS,
                    batch_generator,
                    calculate_adaptive_resolution,
                )

                # Determine backend based on file extension
                ext = Path(path).suffix.lstrip(".").lower()
                backends = FORMAT_BACKENDS.get(ext, ["ASE"])
                backend_name = backends[0]

                loaded_frame_count = 0
                max_particles = 0
                total_expected_frames = None  # Will be set if we know file size

                # Note: extend() internally acquires lock, so we can't wrap this in get_lock()
                # Determine file mode based on format
                # Binary formats need "rb", text formats need "r"
                binary_formats = {"h5", "h5md", "traj", "db", "nc", "hdf5"}
                mode = "rb" if ext in binary_formats else "r"

                # Open file from fsspec filesystem
                with fs.open(path, mode) as f:
                    frame_iterator = None

                    if backend_name == "ZnH5MD":
                        # Use ZnH5MD for H5/H5MD files
                        if step is not None and step <= 0:
                            raise ValueError("Step must be a positive integer for H5MD files.")

                        # Wrap fsspec file object with h5py.File for znh5md compatibility
                        h5_file = h5py.File(f, "r")
                        io = znh5md.IO(file_handle=h5_file)
                        n_frames = len(io)
                        total_expected_frames = (_stop - _start) // _step if step else n_frames

                        if start is not None and start >= n_frames:
                            raise ValueError(f"Start frame {start} exceeds file length")

                        s = slice(start, stop, step)
                        _start, _stop, _step = s.indices(n_frames)
                        frame_iterator = itertools.islice(io, _start, _stop, _step)

                    elif backend_name == "ASE-DB":
                        # ASE database - note: this may not work with file objects
                        temp_vis.log("⚠️ Warning: ASE database format may not support remote filesystems")
                        # Try anyway in case fsspec filesystem supports it
                        db = ase.db.connect(f)
                        n_rows = db.count()

                        if n_rows == 0:
                            temp_vis.log("⚠️ Warning: Database is empty")
                            raise ValueError("Database is empty")

                        if step is not None and step <= 0:
                            raise ValueError("Step must be a positive integer for database files.")

                        _start = start if start is not None else 0
                        _stop = stop if stop is not None else n_rows
                        _step = step if step is not None else 1

                        if _start >= n_rows:
                            raise ValueError(f"Start row {_start} exceeds database size")

                        if _stop > n_rows:
                            _stop = n_rows

                        selected_ids = list(range(_start + 1, _stop + 1, _step))
                        total_expected_frames = len(selected_ids)

                        def db_iterator():
                            for row_id in selected_ids:
                                try:
                                    row = db.get(id=row_id)
                                    yield row.toatoms()
                                except KeyError:
                                    temp_vis.log(
                                        f"⚠️ Warning: Row ID {row_id} not found, skipping"
                                    )
                                    continue

                        frame_iterator = db_iterator()

                    else:
                        # Use ASE for all other formats
                        start_str = str(start) if start is not None else ""
                        stop_str = str(stop) if stop is not None else ""
                        step_str = str(step) if step is not None else ""
                        index_str = f"{start_str}:{stop_str}:{step_str}"

                        frame_iterator = ase.io.iread(f, index=index_str, format=ext or None)

                    # Process frames in batches
                    if frame_iterator:
                        for batch in batch_generator(frame_iterator, batch_size):
                            # Track max particle count
                            for atoms in batch:
                                max_particles = max(max_particles, len(atoms))

                            temp_vis.extend(batch)
                            loaded_frame_count += len(batch)

                            # Update task progress if we know total frames
                            if total_expected_frames and total_expected_frames > 0:
                                progress = (loaded_frame_count / total_expected_frames) * 100
                                task.update(progress=min(progress, 99))  # Cap at 99 until complete

                temp_vis.log(f"✓ Successfully loaded {loaded_frame_count} frames from {path}")

                # Apply adaptive resolution if needed
                if loaded_frame_count > 0 and max_particles > 0:
                    adaptive_resolution = calculate_adaptive_resolution(max_particles)
                    if adaptive_resolution < 16:
                        try:
                            from zndraw.geometries import Sphere

                            sphere = temp_vis.geometries.get("particles")
                            if sphere and isinstance(sphere, Sphere):
                                sphere.resolution = adaptive_resolution
                                temp_vis.geometries["particles"] = sphere
                                temp_vis.log(
                                    f"ℹ️ Reduced particle resolution to {adaptive_resolution} for {max_particles} particles"
                                )
                        except Exception as e:
                            log.warning(f"Failed to apply adaptive resolution: {e}")

                # Task will complete automatically when exiting context manager
                # Return success response (for socketio.call() on server)
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
            # Always disconnect temp_vis to prevent resource leaks
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
            fs_key = f"global:{fs_name}" if public else f"room:{self.zndraw.room}:{fs_name}"
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
