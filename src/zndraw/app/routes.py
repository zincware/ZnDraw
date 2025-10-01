import json
import logging

import msgpack
import zarr
from flask import Response, current_app, request

from zndraw.storage import ZarrStorageSequence, decode_data, encode_data
from zndraw.server import socketio
from .constants import SocketEvents
from .redis_keys import ExtensionKeys
from .worker_stats import WorkerStats
from .queue_manager import emit_queue_update
import traceback


from . import main

# --- Logging Setup ---
log = logging.getLogger(__name__)


# TODO: move to utils
def get_lock_key(room: str, target: str) -> str:
    """Constructs a standardized Redis key for a lock."""
    return f"room:{room}:lock:{target}"


def get_zarr_store_path(room_id: str) -> str:
    """Returns the path to the Zarr store for a given room."""
    # In a real app, this path might come from a config file.
    return f"data/{room_id}.zarr"

@main.route("/internal/emit", methods=["POST"])
def internal_emit():
    """Internal endpoint to emit Socket.IO events. Secured via a shared secret."""
    data = request.get_json()

    event = data.get("event")
    sid = data.get("sid")
    payload = data.get("data", {})

    if not event or not sid:
        return {"error": "Event and sid are required"}, 400

    socketio.emit(event, payload, to=sid)
    return {"success": True}


@main.route("/api/frames/<string:room_id>", methods=["POST"])
def get_frames(room_id):
    """Serves multiple frames' data from the room's Zarr store using either indices or slice parameters."""
    r = current_app.extensions["redis"]
    try:
        # Parse the request data
        request_data = request.get_json()
        if request_data is None:
            return {"error": "Request body required"}, 400

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)
        storage = ZarrStorageSequence(root)

        # Get logical-to-physical mapping from Redis
        indices_key = f"room:{room_id}:trajectory:indices"
        frame_mapping = r.zrange(indices_key, 0, -1)

        if not frame_mapping:
            return {"error": "No frames found in room"}, 404

        max_frame = len(frame_mapping) - 1

        # Determine frame indices based on request parameters
        if "indices" in request_data:
            # Direct list of indices
            frame_indices = request_data["indices"]
            if not isinstance(frame_indices, list):
                return {"error": "Indices must be a list"}, 400

            # Validate frame indices
            for frame_id in frame_indices:
                if (
                    not isinstance(frame_id, int)
                    or frame_id < 0
                    or frame_id > max_frame
                ):
                    error_data = {
                        "error": f"Invalid frame index {frame_id}, valid range: 0-{max_frame}",
                        "type": "IndexError",
                    }
                    return Response(
                        json.dumps(error_data),
                        status=404,
                        content_type="application/json",
                    )

        else:
            # Default to slice behavior for any remaining cases (including empty payload)
            # This handles slice parameters and slice(None, None, None) which sends empty payload
            start = request_data.get("start", 0)
            stop = request_data.get("stop", len(frame_mapping))
            step = request_data.get("step", 1)

            # Validate slice parameters
            if not all(isinstance(x, int) for x in [start, stop, step]):
                return {"error": "start, stop, and step must be integers"}, 400

            if step == 0:
                return {"error": "step cannot be zero"}, 400

            # Generate frame indices from slice
            try:
                frame_indices = list(range(start, stop, step))
                # Filter out invalid indices
                frame_indices = [i for i in frame_indices if 0 <= i <= max_frame]
            except ValueError as e:
                return {"error": f"Invalid slice parameters: {e}"}, 400

        # Get keys parameter if specified
        requested_keys = request_data.get("keys")
        print(f"get_frames called with keys: {requested_keys}")
        # error_data = {"error": f"Key(s) not found: {', '.join(sorted(missing_keys))}", "type": "KeyError"}
        # return Response(json.dumps(error_data), status=404, content_type='application/json')

        # TODO: requested keys and KeyError handling
        # TODO: instead of iterate load at once
        try:
            frames_data = []
            for frame_id in frame_indices:
                # Get the physical index for this logical frame
                physical_index = int(frame_mapping[frame_id])
                frame_data = storage.get(physical_index, keys=requested_keys)
                frames_data.append(encode_data(frame_data))
        except KeyError as e:
            error_data = {
                "error": f"Key(s) not found: {e}",
                "type": "KeyError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )
        except IndexError as e:
            error_data = {
                "error": f"Frame index out of range: {e}",
                "type": "IndexError",
            }
            return Response(
                json.dumps(error_data), status=404, content_type="application/json"
            )

        packed_data = msgpack.packb(frames_data)
        return Response(packed_data, content_type="application/octet-stream")
    except Exception as e:
        error_data = {
            "error": f"Server error: {e}",
            "type": type(e).__name__,
            "success": False,
        }
        print(traceback.format_exc())
        return Response(
            json.dumps(error_data), status=500, content_type="application/json"
        )


@main.route("/api/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    """Appends a new frame. Authorized via a short-lived Bearer token."""
    r = current_app.extensions["redis"]
    auth_header = request.headers.get("Authorization")
    if not auth_header or not auth_header.startswith("Bearer "):
        return {"error": "Authorization token is missing or invalid"}, 401

    token = auth_header.split(" ")[1]
    token_key = f"room:{room_id}:upload_token:{token}"

    # Get token metadata
    token_data = r.hgetall(token_key)
    if not token_data:
        return {"error": "Token is invalid or has expired"}, 403

    sid_from_token = token_data.get("sid")
    action = token_data.get("action", "append")
    target_frame_id = (
        int(token_data.get("frame_id", -1))
        if token_data.get("frame_id") != "-1"
        else None
    )

    r.delete(token_key)  # Invalidate the token after first use

    lock_key = get_lock_key(room_id, "trajectory:meta")
    if r.get(lock_key) != sid_from_token:
        return {"error": "Client does not hold the trajectory lock"}, 403

    try:
        # Unpack the msgpack data
        serialized_data = msgpack.unpackb(request.data, strict_map_key=False)

        store_path = get_zarr_store_path(room_id)
        root = zarr.group(store_path)
        storage = ZarrStorageSequence(root)

        indices_key = f"room:{room_id}:trajectory:indices"

        if action == "replace":
            frame_mapping = r.zrange(indices_key, 0, -1)

            # Validate target_frame_id
            if target_frame_id is None or not (
                0 <= target_frame_id < len(frame_mapping)
            ):
                return {
                    "error": f"Invalid or missing frame_id for replace. Valid range: 0-{len(frame_mapping) - 1}"
                }, 404

            new_physical_index = len(storage)
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            old_physical_index_to_unmap = frame_mapping[target_frame_id]

            pipeline = r.pipeline()
            pipeline.zrem(indices_key, old_physical_index_to_unmap)
            pipeline.zadd(indices_key, {str(new_physical_index): target_frame_id})
            pipeline.execute()

            log.info(
                f"Replaced frame {target_frame_id} (old physical: {old_physical_index_to_unmap}, new physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "replaced_frame": target_frame_id}

        elif action == "extend":
            # Extend operation: add multiple frames in one go
            if not isinstance(serialized_data, list):
                return {
                    "error": "For extend action, data must be a list of frame dictionaries"
                }, 400

            # 1. Determine starting logical and physical positions
            start_logical_pos = r.zcard(indices_key)
            start_physical_pos = len(storage)
            num_frames = len(serialized_data)

            # 2. Decode all frames and extend the physical storage
            decoded_frames = [decode_data(frame) for frame in serialized_data]
            storage.extend(decoded_frames)

            # 3. Create the new logical-to-physical mapping for all new frames
            new_mapping = {
                str(start_physical_pos + i): start_logical_pos + i
                for i in range(num_frames)
            }
            if new_mapping:
                r.zadd(indices_key, new_mapping)

            # 4. Prepare response data
            new_indices = list(range(start_logical_pos, start_logical_pos + num_frames))

            log.info(
                f"Extended trajectory with {num_frames} frames (physical: {start_physical_pos}-{start_physical_pos + num_frames - 1}) to room '{room_id}'"
            )
            return {"success": True, "new_indices": new_indices}

        elif action == "insert":
            # Insert operation: add a new physical frame and insert it into the logical sequence
            insert_position = int(token_data.get("insert_position", 0))
            current_length = r.zcard(indices_key)

            if not (0 <= insert_position <= current_length):
                return {
                    "error": f"Insert position {insert_position} out of range [0, {current_length}]"
                }, 400

            # 1. Determine the new physical index (always at the end of the physical store)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data to the physical store
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Atomically shift existing logical indices and insert the new one
            pipeline = r.pipeline()
            # Get members whose scores (logical positions) need to be incremented
            members_to_shift = r.zrangebyscore(indices_key, insert_position, "+inf")
            if members_to_shift:
                for member in members_to_shift:
                    pipeline.zincrby(indices_key, 1, member)

            # Add the new frame at the correct logical position
            pipeline.zadd(indices_key, {str(new_physical_index): insert_position})
            pipeline.execute()

            log.info(
                f"Inserted frame at position {insert_position} (physical: {new_physical_index}) in room '{room_id}'"
            )
            return {"success": True, "inserted_position": insert_position}

        elif action == "append":
            # Append operation: add a new frame to the end of the logical sequence
            # 1. Determine the new logical and physical positions
            logical_position = r.zcard(indices_key)
            new_physical_index = len(storage)

            # 2. Decode and append the new frame data
            decoded_data = decode_data(serialized_data)
            storage.append(decoded_data)

            # 3. Add the new physical index to the logical sequence
            r.zadd(indices_key, {str(new_physical_index): logical_position})

            log.info(
                f"Appended frame {logical_position} (physical: {new_physical_index}) to room '{room_id}'"
            )
            return {"success": True, "new_index": logical_position}

        else:
            # Default case for any unknown actions
            return {"error": f"The requested action '{action}' is not supported."}, 400
    except Exception as e:
        log.error(f"Failed to write to Zarr store: {e}")
        print(traceback.format_exc())
        return {"error": "Failed to write to data store"}, 500


@main.route("/api/exit")
def exit_app():
    """Endpoint to gracefully shut down the server. Secured via a shared secret."""
    socketio.stop()
    return {"success": True}

@main.route("/api/rooms/<string:room_id>/schema/<string:category>", methods=["GET"])
def get_room_schema(room_id: str, category: str):
    """Get the schema for a specific room with worker and queue statistics.

    Returns schema along with metadata about each extension:
    - provider: "celery" for server-side extensions, or count of registered workers
    - queueLength: number of queued tasks for this extension
    - idleWorkers: number of idle workers available
    - progressingWorkers: number of workers currently processing tasks
    """
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    # Map category strings to the corresponding imported objects
    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
    }

    if category not in category_map:
        return {"error": f"Unknown schema category '{category}'"}

    redis_client = current_app.extensions["redis"]
    schema = {}

    # Add server-provided extensions (Celery-based)
    for name, cls in category_map[category].items():
        schema[name] = {
            "schema": cls.model_json_schema(),
            "provider": "celery",
            "queueLength": 0,
            "idleWorkers": 0,
            "progressingWorkers": 0
        }

    # Add client-provided extensions from Redis
    schema_key = ExtensionKeys.schema_key(room_id, category)
    redis_schema = redis_client.hgetall(schema_key)

    for name, sch_str in redis_schema.items():
        sch = json.loads(sch_str)

        # Get worker statistics for this extension
        keys = ExtensionKeys.for_extension(room_id, category, name)
        stats = WorkerStats.fetch(redis_client, keys)

        if name in schema:
            if schema[name]["schema"] != sch:
                print(
                    f"Warning: {category.capitalize()} extension '{name}' schema "
                    "in Redis differs from server schema."
                )
        else:
            schema[name] = {
                "schema": sch,
                "provider": stats.total_workers,  # Number of workers for client extensions
                **stats.to_dict()
            }

    return schema


@main.route("/api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>", methods=["POST"])
def log_room_extension(room_id: str, category: str, extension: str):
    """Logs a user extension action in the room's action log."""
    user_id = request.args.get("userId")

    if user_id is None:
        return {"error": "User ID is required"}, 400

    data = request.json
    if data is None:
        return {"error": "Request body required"}, 400

    print(f"Logging extension for room {room_id}: category={category}, extension={extension}, data={json.dumps(data)}")

    # store in redis
    redis_client = current_app.extensions["redis"]
    # Store the entire extension data as a JSON string
    redis_client.hset(
        f"room:{room_id}:user:{user_id}:{category}", extension, json.dumps(data)
    )

    # Check if this is a server-side (Celery) extension
    from zndraw.extensions.modifiers import modifiers
    from zndraw.extensions.selections import selections
    from zndraw.settings import settings

    category_map = {
        "selections": selections,
        "modifiers": modifiers,
        "settings": settings,
    }

    is_celery_extension = (
        category in category_map and extension in category_map[category]
    )

    queue_position = 0

    if is_celery_extension:
        # Dispatch to Celery task for server-side execution
        from zndraw.app.tasks import run_extension_task

        log.info(
            f"Dispatching Celery task for {category}/{extension} in room {room_id}"
        )
        run_extension_task.delay(room_id, category, extension, data, user_id)
        print(
            f"Dispatched Celery task for user {user_id}, category {category}, extension {extension}, room {room_id}"
        )
    else:
        # send data via `task:run` emit to client workers
        keys = ExtensionKeys.for_extension(room_id, category, extension)

        worker_sid = redis_client.spop(keys.idle_workers)
        if worker_sid is not None:
            socketio.emit(
                SocketEvents.TASK_RUN,
                {"data": data, "extension": extension, "category": category, "room": room_id},
                to=worker_sid,
            )
            print(
                f"Emitted task:run to worker {worker_sid} for user {user_id}, category {category}, extension {extension}, room {room_id}"
            )
            redis_client.sadd(keys.progressing_workers, worker_sid)

            # Notify all clients in room about worker state change
            emit_queue_update(redis_client, room_id, category, extension, socketio)
        else:
            # add to queue
            redis_client.rpush(
                keys.queue, json.dumps({"user_id": user_id, "data": data, "room": room_id})
            )
            print(
                f"No idle workers, queued task for user {user_id}, category {category}, extension {extension}, room {room_id}"
            )
            queue_position = redis_client.llen(keys.queue)

            # Notify all clients in room about queue update
            emit_queue_update(redis_client, room_id, category, extension, socketio)

    print(
        f"Emitting invalidate for user {user_id}, category {category}, extension {extension}, room {room_id} to user:{user_id}"
    )
    socketio.emit(
        SocketEvents.INVALIDATE,
        {"userId": user_id, "category": category, "extension": extension, "roomId": room_id},
        to=f"user:{user_id}",
    )
    return {"status": "success", "queuePosition": queue_position}, 200

@main.route("/api/rooms/<string:room_id>/extension-data/<string:category>/<string:extension>", methods=["GET"])
def get_extension_data(room_id: str, category: str, extension: str):
    user_id = request.args.get("userId")
    print(f"get_extension_data called with userId={user_id}, category={category}, extension={extension} for room {room_id}")

    if not user_id:
        return {"error": "User ID is required"}, 400

    redis_client = current_app.extensions["redis"]
    extension_data = redis_client.hget(f"room:{room_id}:user:{user_id}:{category}", extension)
    if extension_data is None:
        return {"data": None}, 200
    extension_data = json.loads(extension_data)
    return {"data": extension_data}, 200

