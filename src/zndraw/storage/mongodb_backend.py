"""MongoDB storage backend for distributed frame storage."""

import logging

import msgpack
import numpy as np
from bson import Binary
from pymongo import MongoClient, ReturnDocument
from pymongo.collection import Collection

from .base import StorageBackend

log = logging.getLogger(__name__)


class MongoDBStorageBackend(StorageBackend):
    """MongoDB storage backend for distributed frame storage.

    This backend stores frames in MongoDB with:
    - One collection per room
    - Frame index as document _id for efficient lookups
    - Raw msgpack bytes stored as BSON Binary

    Document schema:
        {
            "_id": <frame_index>,
            "data": Binary(<msgpack_bytes>)
        }
    """

    def __init__(self, uri: str, database: str, room_id: str):
        """Initialize MongoDB storage backend.

        Parameters
        ----------
        uri : str
            MongoDB connection URI (e.g., "mongodb://root:example@localhost:27017/")
        database : str
            Database name
        room_id : str
            Room identifier (used as collection name)
        """
        self.uri = uri
        self.database_name = database
        self.room_id = room_id

        self.client: MongoClient = MongoClient(
            uri,
            serverSelectionTimeoutMS=5000,
            connectTimeoutMS=10000,
            maxPoolSize=50,
        )
        try:
            # Validate connection
            self.client.admin.command("ping")

            self.db = self.client[database]
            self.collection: Collection = self.db[room_id]

            log.debug(
                f"Initialized MongoDBStorageBackend: "
                f"database='{database}', collection='{room_id}'"
            )
        except Exception:
            self.client.close()
            raise

    def close(self) -> None:
        """Close MongoDB connection and release resources."""
        self.client.close()
        log.debug(f"Closed MongoDB connection for room '{self.room_id}'")

    def get(
        self,
        index: int | list[int] | slice | np.ndarray,
        keys: list[str] | None = None,
    ) -> dict[bytes, bytes] | list[dict[bytes, bytes]]:
        """Get frame(s) with optional key filtering.

        Parameters
        ----------
        index : int | list[int] | slice | np.ndarray
            Frame index/indices to retrieve
        keys : list[str] | None
            Optional list of keys to filter (e.g., ["arrays.positions", "info.energy"])

        Returns
        -------
        dict[bytes, bytes] | list[dict[bytes, bytes]]
            For single index: dict[bytes, bytes] with msgpack-serialized values
            For multiple indices: list of dicts
        """
        # Handle numpy arrays and scalars
        if isinstance(index, np.ndarray):
            if index.ndim == 0:
                index = int(index.item())
            else:
                index = index.tolist()
        elif isinstance(index, np.integer):
            index = int(index)

        if isinstance(index, slice):
            index = list(range(*index.indices(len(self))))

        is_single = isinstance(index, int)
        if is_single:
            index = [index]

        # Validate bounds
        length = len(self)
        for i in index:
            if i < -length or i >= length:
                raise IndexError(
                    f"Index {i} is out of bounds for storage of length {length}"
                )

        # Normalize negative indices
        normalized_indices = [i if i >= 0 else length + i for i in index]

        # Query MongoDB - use $in for batch retrieval
        docs = list(
            self.collection.find(
                {"_id": {"$in": normalized_indices}},
                {"_id": 1, "data": 1},
            )
        )

        # Create a lookup for fast access
        doc_map = {doc["_id"]: doc["data"] for doc in docs}

        # Build results in order of requested indices
        results = []
        keys_bytes = None
        if keys is not None:
            keys_bytes = set(k.encode() for k in keys)

        for i in normalized_indices:
            if i not in doc_map:
                raise IndexError(f"Frame {i} not found in storage")

            # Unpack the msgpack dict from BSON Binary
            msgpack_dict = msgpack.unpackb(bytes(doc_map[i]))

            # Filter by keys if requested
            if keys_bytes is not None:
                msgpack_dict = {
                    k: v for k, v in msgpack_dict.items() if k in keys_bytes
                }

            results.append(msgpack_dict)

        if is_single:
            return results[0]

        return results

    def extend(self, values: list[dict[bytes, bytes]]) -> None:
        """Extend storage with multiple frames (batch write).

        Parameters
        ----------
        values : list[dict[bytes, bytes]]
            List of frame dictionaries with msgpack-serialized values
            Keys should be like b"arrays.positions", b"info.energy"
        """
        if not values:
            return

        # Use atomic counter for sequential IDs to prevent race conditions
        counters = self.db["_counters"]
        result = counters.find_one_and_update(
            {"_id": f"{self.room_id}_frame_counter"},
            {"$inc": {"value": len(values)}},
            upsert=True,
            return_document=ReturnDocument.AFTER,
        )
        start_idx = result["value"] - len(values)

        # Prepare documents for bulk insert
        documents = []
        for i, frame_data in enumerate(values):
            # Pack the entire frame dict as msgpack bytes
            packed = msgpack.packb(frame_data)
            documents.append({"_id": start_idx + i, "data": Binary(packed)})

        # Bulk insert
        self.collection.insert_many(documents)

        log.debug(
            f"Extended storage with {len(values)} frames (indices {start_idx}-{start_idx + len(values) - 1})"
        )

    def get_available_keys(self, index: int) -> list[str]:
        """List all keys available for a frame.

        Parameters
        ----------
        index : int
            Frame index

        Returns
        -------
        list[str]
            List of available keys (e.g., ["arrays.positions", "info.energy"])
        """
        length = len(self)
        if index < -length or index >= length:
            raise IndexError(
                f"Index {index} is out of bounds for storage of length {length}"
            )

        # Normalize negative index
        if index < 0:
            index = length + index

        doc = self.collection.find_one({"_id": index}, {"data": 1})
        if doc is None:
            raise IndexError(f"Frame {index} not found in storage")

        # Unpack and get keys
        msgpack_dict = msgpack.unpackb(bytes(doc["data"]))
        return [k.decode() for k in msgpack_dict.keys()]

    def __len__(self) -> int:
        """Return the number of frames in storage."""
        return self.collection.count_documents({})
