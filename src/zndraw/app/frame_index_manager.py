"""Frame indexing manager using gap-based sorted set scores.

This module implements a scalable indexing strategy for Redis sorted sets
that avoids O(N) reindexing operations by using floating-point scores with gaps.
"""

import typing as t

from redis import Redis


class FrameIndexManager:
    """Manages frame indices in Redis sorted sets using gap-based scoring.

    Instead of contiguous integer scores (0, 1, 2, ...), this uses floating-point
    scores with gaps to enable O(log N) insertions without shifting all subsequent frames.

    Operations:
    - Append: O(log N) - get last score and add 1.0
    - Insert: O(log N) - average of surrounding scores
    - Delete: O(log N) - simple ZREM
    - Renormalize: O(N) - reassign to contiguous integers (infrequent)
    """

    def __init__(self, redis_client: Redis, key: str):
        """Initialize the frame index manager.

        Args:
            redis_client: Redis connection instance
            key: Redis sorted set key for frame indices
        """
        self.redis = redis_client
        self.key = key

    def append(self, member: str) -> float:
        """Append a frame to the end of the sequence.

        Args:
            member: The member to add to the sorted set

        Returns:
            The score assigned to the new member
        """
        # Get the last element's score
        last_elements = self.redis.zrange(self.key, -1, -1, withscores=True)

        if not last_elements:
            # First frame
            new_score = 1.0
        else:
            # Add 1.0 to the last score
            new_score = last_elements[0][1] + 1.0

        self.redis.zadd(self.key, {member: new_score})
        return new_score

    def append_batch(self, members: t.List[str]) -> t.List[float]:
        """Append multiple frames to the end of the sequence in a single batch.

        This is much more efficient than calling append() multiple times,
        as it reduces N Redis operations to just 2 (one ZRANGE and one ZADD).

        Args:
            members: List of members to add to the sorted set

        Returns:
            List of scores assigned to the new members
        """
        if not members:
            return []

        # Get the last element's score once
        last_elements = self.redis.zrange(self.key, -1, -1, withscores=True)

        if not last_elements:
            # First batch of frames
            start_score = 1.0
        else:
            # Add 1.0 to the last score
            start_score = last_elements[0][1] + 1.0

        # Build the mapping for all members
        mapping = {member: start_score + i for i, member in enumerate(members)}

        # Add all members in a single Redis call
        self.redis.zadd(self.key, mapping)

        # Return the list of scores
        return [start_score + i for i in range(len(members))]

    def insert(self, position: int, member: str) -> float:
        """Insert a frame at a specific logical position.

        Args:
            position: Logical position (0-indexed) where to insert
            member: The member to add to the sorted set

        Returns:
            The score assigned to the new member
        """
        total_count = self.redis.zcard(self.key)

        if position < 0:
            position = max(0, total_count + position + 1)

        # Get surrounding elements
        prev_element = None
        next_element = None

        if position > 0:
            prev_element = self.redis.zrange(
                self.key, position - 1, position - 1, withscores=True
            )

        if position < total_count:
            next_element = self.redis.zrange(
                self.key, position, position, withscores=True
            )

        # Calculate new score based on surrounding elements
        new_score = self._calculate_insert_score(prev_element, next_element)

        self.redis.zadd(self.key, {member: new_score})
        return new_score

    def _calculate_insert_score(
        self,
        prev_element: t.Optional[t.List[t.Tuple[bytes, float]]],
        next_element: t.Optional[t.List[t.Tuple[bytes, float]]],
    ) -> float:
        """Calculate the score for a new element based on its neighbors.

        Args:
            prev_element: The element before the insertion point (if any)
            next_element: The element after the insertion point (if any)

        Returns:
            The calculated score for the new element
        """
        if not prev_element and not next_element:
            # First frame
            return 1.0
        elif not prev_element:
            # Insert before first element
            return next_element[0][1] - 1.0
        elif not next_element:
            # Append after last element
            return prev_element[0][1] + 1.0
        else:
            # Insert between two elements - use average
            return (prev_element[0][1] + next_element[0][1]) / 2.0

    def delete(self, member: str) -> int:
        """Delete a frame from the sequence.

        Args:
            member: The member to remove from the sorted set

        Returns:
            Number of elements removed (0 or 1)
        """
        return self.redis.zrem(self.key, member)

    def delete_at_position(self, position: int) -> int:
        """Delete a frame at a specific logical position.

        Args:
            position: Logical position (0-indexed) to delete

        Returns:
            Number of elements removed (0 or 1)
        """
        elements = self.redis.zrange(self.key, position, position)
        if not elements:
            return 0
        return self.redis.zrem(self.key, elements[0])

    def renormalize(self) -> int:
        """Reassign all scores to contiguous integers starting from 0.

        This is an O(N) operation that should be called infrequently when
        precision drift becomes an issue (e.g., scores too close together).

        Returns:
            Number of frames renormalized
        """
        # Get all elements with their current scores
        all_elements = self.redis.zrange(self.key, 0, -1, withscores=True)

        if not all_elements:
            return 0

        # Create new mapping with contiguous integer scores
        new_mapping = {
            member: float(idx) for idx, (member, _) in enumerate(all_elements)
        }

        # Replace the entire sorted set atomically using a pipeline
        pipeline = self.redis.pipeline()
        pipeline.delete(self.key)
        pipeline.zadd(self.key, new_mapping)
        pipeline.execute()

        return len(all_elements)

    def __len__(self) -> int:
        """Get the total number of frames.

        Returns:
            Number of frames in the sequence
        """
        return self.redis.zcard(self.key)

    def __getitem__(self, index: int) -> str:
        """Get frame at a specific index.

        Args:
            index: Frame index (0-indexed, negative indices supported)

        Returns:
            Frame key at the specified index

        Raises:
            IndexError: If index is out of range
        """
        result = self.redis.zrange(self.key, index, index)
        if not result:
            raise IndexError(f"Frame index {index} out of range")
        # Decode bytes to string if necessary
        frame_key = result[0]
        if isinstance(frame_key, bytes):
            frame_key = frame_key.decode("utf-8")
        return frame_key

    def get_count(self) -> int:
        """Get the total number of frames.

        Returns:
            Number of frames in the sequence
        """
        return self.redis.zcard(self.key)

    def get_all(self, withscores: bool = False) -> t.List:
        """Get all frames in order.

        Args:
            withscores: Whether to include scores in the result

        Returns:
            List of frames (with scores if requested)
        """
        return self.redis.zrange(self.key, 0, -1, withscores=withscores)

    def get_range(self, start: int, end: int, withscores: bool = False) -> t.List:
        """Get a range of frames by logical position.

        Args:
            start: Start position (inclusive, 0-indexed)
            end: End position (inclusive, -1 for last)
            withscores: Whether to include scores in the result

        Returns:
            List of frames in the specified range
        """
        return self.redis.zrange(self.key, start, end, withscores=withscores)

    def get_by_indices(self, indices: t.List[int], withscores: bool = False) -> t.List:
        """
        Get a list of frames by their specific logical indices (ranks).
        This is highly efficient for fetching non-contiguous frames.
        O(M * log(N)) complexity using a pipeline.

        Args:
            indices: A list of 0-indexed logical positions to fetch.
            withscores: Whether to include scores in the result.

        Returns:
            List of frames (and scores if requested) in the same order as the input indices.
        """

        if not indices:
            return []

        with self.redis.pipeline() as pipe:
            for idx in indices:
                pipe.zrange(self.key, idx, idx, withscores=withscores)
            results = pipe.execute()

        # The pipeline returns a list of lists. Flatten it.
        # e.g., [[b'item1'], [b'item5'], []] -> [b'item1', b'item5']
        return [item[0] for item in results if item]
