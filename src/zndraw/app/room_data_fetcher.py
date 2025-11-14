"""Batched room data fetcher for efficient room listing.

This module provides utilities to fetch data for multiple rooms efficiently
using Redis pipelining to avoid N+1 query patterns.
"""

from zndraw.app.metadata_manager import RoomMetadataManager
from zndraw.app.redis_keys import RoomKeys
from zndraw.app.route_utils import get_lock_key


class BatchedRoomDataFetcher:
    """Efficiently fetch data for multiple rooms using Redis pipeline.

    This class batches Redis queries to avoid N+1 query patterns when
    listing multiple rooms. Instead of making 6+ queries per room,
    it makes a small constant number of queries total.
    """

    def fetch_rooms_data(
        self, redis_client, room_service, room_ids: list[str]
    ) -> dict[str, dict]:
        """Fetch all data for multiple rooms in batched queries.

        Parameters
        ----------
        redis_client : Redis
            Redis client instance
        room_service : RoomService
            Room service for frame counts
        room_ids : list[str]
            List of room IDs to fetch data for

        Returns
        -------
        dict[str, dict]
            Mapping of room_id to room data dict with keys:
            - frameCount: int
            - description: str | None
            - locked: bool
            - hidden: bool
            - metadataLocked: bool
            - metadata: dict[str, str]

        Notes
        -----
        This method uses Redis pipelining to batch queries and reduce
        network round trips from O(n * queries_per_room) to O(1).

        Performance: For 100 rooms, reduces ~600 queries to ~5-10 queries.
        """
        if not room_ids:
            return {}

        # Step 1: Get frame counts (LMDB access, cannot batch Redis but is local)
        frame_counts = room_service.get_frame_counts_batch(room_ids)

        # Step 2: Build pipeline for all Redis queries
        pipe = redis_client.pipeline()

        # Batch simple keys using MGET (3 queries total instead of n*3)
        desc_keys = [RoomKeys(rid).description() for rid in room_ids]
        lock_keys = [RoomKeys(rid).locked() for rid in room_ids]
        hide_keys = [RoomKeys(rid).hidden() for rid in room_ids]

        pipe.mget(desc_keys)
        pipe.mget(lock_keys)
        pipe.mget(hide_keys)

        # Batch metadata locks (n queries but in single pipeline)
        for room_id in room_ids:
            lock_key = get_lock_key(room_id, "trajectory:meta")
            pipe.get(lock_key)

        # Batch metadata hashes (n queries but in single pipeline)
        for room_id in room_ids:
            metadata_key = f"room:{room_id}:metadata"
            pipe.hgetall(metadata_key)

        # Execute all queries at once
        results = pipe.execute()

        # Parse results
        descriptions = results[0]
        lockeds = results[1]
        hiddens = results[2]
        metadata_locks = results[3 : 3 + len(room_ids)]
        metadatas = results[3 + len(room_ids) :]

        # Build structured response
        rooms_data = {}
        for i, room_id in enumerate(room_ids):
            # Handle both bytes and string responses from Redis
            locked_val = lockeds[i]
            hidden_val = hiddens[i]

            if isinstance(locked_val, bytes):
                locked = locked_val == b"1"
            else:
                locked = locked_val == "1"

            if isinstance(hidden_val, bytes):
                hidden = hidden_val == b"1"
            else:
                hidden = hidden_val == "1"

            rooms_data[room_id] = {
                "frameCount": frame_counts.get(room_id, 0),
                "description": descriptions[i],
                "locked": locked,
                "hidden": hidden,
                "metadataLocked": metadata_locks[i] is not None,
                "metadata": metadatas[i] if metadatas[i] else {},
            }

        return rooms_data
