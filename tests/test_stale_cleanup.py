"""Tests for stale worker detection and cleanup after server restart.

Tests cover:
1. ClusterHeartbeat - initializing, refreshing, and using heartbeat timestamps
2. Stale extension worker detection and cleanup
3. Stale filesystem worker detection and cleanup
4. WorkerStats filtering of stale workers
5. Job dispatcher filtering of stale workers
"""

import time
from unittest.mock import MagicMock

from zndraw.app.cluster_heartbeat import ClusterHeartbeat
from zndraw.app.redis_keys import ClusterKeys, ExtensionKeys, FilesystemKeys
from zndraw.app.stale_cleanup import (
    cleanup_all_stale_workers_for_category,
    cleanup_stale_extension_workers,
    cleanup_stale_filesystem_worker,
)
from zndraw.app.worker_stats import WorkerStats


class TestClusterHeartbeat:
    """Tests for the ClusterHeartbeat class."""

    def test_start_creates_new_heartbeat(self, redis_client):
        """Test that start() creates a new heartbeat when none exists."""
        heartbeat = ClusterHeartbeat(redis_client, ttl_seconds=60)

        before = time.time()
        ts = heartbeat.start()
        after = time.time()

        try:
            # Timestamp should be close to current time
            assert before <= ts <= after

            # Heartbeat key should exist in Redis
            stored = redis_client.get(ClusterKeys.HEARTBEAT)
            assert stored is not None
            assert float(stored) == ts
        finally:
            heartbeat.stop()

    def test_start_reuses_existing_heartbeat(self, redis_client):
        """Test that start() reuses existing heartbeat from another instance."""
        # Set up existing heartbeat
        existing_ts = time.time() - 30  # 30 seconds ago
        redis_client.set(ClusterKeys.HEARTBEAT, existing_ts, ex=60)

        heartbeat = ClusterHeartbeat(redis_client, ttl_seconds=60)
        ts = heartbeat.start()

        try:
            # Should return the existing timestamp
            assert ts == existing_ts
        finally:
            heartbeat.stop()

    def test_get_heartbeat_timestamp(self, redis_client):
        """Test getting heartbeat timestamp from Redis."""
        expected = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, expected)

        actual = ClusterHeartbeat.get_heartbeat_timestamp(redis_client)

        assert actual == expected

    def test_get_heartbeat_timestamp_returns_none_when_missing(self, redis_client):
        """Test that get_heartbeat_timestamp returns None when no heartbeat."""
        actual = ClusterHeartbeat.get_heartbeat_timestamp(redis_client)

        assert actual is None

    def test_is_worker_stale_returns_true_for_old_timestamp(self, redis_client):
        """Test that is_worker_stale returns True for worker registered before heartbeat."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        worker_ts = heartbeat_ts - 100  # Registered before heartbeat

        assert ClusterHeartbeat.is_worker_stale(redis_client, worker_ts) is True
        assert ClusterHeartbeat.is_worker_stale(redis_client, str(worker_ts)) is True

    def test_is_worker_stale_returns_false_for_new_timestamp(self, redis_client):
        """Test that is_worker_stale returns False for worker registered after heartbeat."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        worker_ts = heartbeat_ts + 100  # Registered after heartbeat

        assert ClusterHeartbeat.is_worker_stale(redis_client, worker_ts) is False

    def test_is_worker_stale_returns_false_when_no_heartbeat(self, redis_client):
        """Test that is_worker_stale returns False when no heartbeat exists."""
        worker_ts = time.time()

        assert ClusterHeartbeat.is_worker_stale(redis_client, worker_ts) is False

    def test_stop_stops_thread(self, redis_client):
        """Test that stop() properly stops the refresh thread."""
        heartbeat = ClusterHeartbeat(redis_client, ttl_seconds=60)
        heartbeat.start()

        assert heartbeat._thread is not None
        assert heartbeat._thread.is_alive()

        heartbeat.stop()

        assert heartbeat._thread is None or not heartbeat._thread.is_alive()


class TestStaleExtensionCleanup:
    """Tests for stale extension worker cleanup."""

    def test_cleanup_removes_stale_worker(self, redis_client):
        """Test that cleanup_stale_extension_workers removes stale workers."""
        # Set up heartbeat (current time)
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        # Register a stale worker (timestamp before heartbeat)
        keys = ExtensionKeys.for_extension(room_id, category, extension)
        stale_worker_id = "stale-worker-123"
        stale_ts = heartbeat_ts - 100

        redis_client.hset(keys.workers, stale_worker_id, stale_ts)
        redis_client.hset(keys.schema, extension, '{"type": "object"}')

        # Mock socketio
        mock_socketio = MagicMock()

        # Run cleanup
        cleaned = cleanup_stale_extension_workers(
            redis_client, mock_socketio, room_id, category, extension
        )

        # Worker should be removed
        assert stale_worker_id in cleaned
        assert redis_client.hexists(keys.workers, stale_worker_id) is False

        # Schema should be deleted (no more workers, no pending jobs)
        assert redis_client.hexists(keys.schema, extension) is False

        # Should emit schema invalidation
        mock_socketio.emit.assert_called()

    def test_cleanup_keeps_valid_worker(self, redis_client):
        """Test that cleanup_stale_extension_workers keeps valid workers."""
        # Set up heartbeat
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        # Register a valid worker (timestamp after heartbeat)
        keys = ExtensionKeys.for_extension(room_id, category, extension)
        valid_worker_id = "valid-worker-456"
        valid_ts = heartbeat_ts + 10

        redis_client.hset(keys.workers, valid_worker_id, valid_ts)
        redis_client.hset(keys.schema, extension, '{"type": "object"}')

        mock_socketio = MagicMock()

        # Run cleanup
        cleaned = cleanup_stale_extension_workers(
            redis_client, mock_socketio, room_id, category, extension
        )

        # No workers should be cleaned
        assert cleaned == []
        assert redis_client.hexists(keys.workers, valid_worker_id) is True
        assert redis_client.hexists(keys.schema, extension) is True

    def test_cleanup_mixed_stale_and_valid_workers(self, redis_client):
        """Test cleanup with both stale and valid workers."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add stale worker
        stale_worker = "stale-worker"
        redis_client.hset(keys.workers, stale_worker, heartbeat_ts - 100)

        # Add valid worker
        valid_worker = "valid-worker"
        redis_client.hset(keys.workers, valid_worker, heartbeat_ts + 10)

        redis_client.hset(keys.schema, extension, '{"type": "object"}')

        mock_socketio = MagicMock()

        cleaned = cleanup_stale_extension_workers(
            redis_client, mock_socketio, room_id, category, extension
        )

        # Only stale worker should be cleaned
        assert stale_worker in cleaned
        assert valid_worker not in cleaned

        # Valid worker should remain
        assert redis_client.hexists(keys.workers, valid_worker) is True
        assert redis_client.hexists(keys.workers, stale_worker) is False

        # Schema should remain (still has valid worker)
        assert redis_client.hexists(keys.schema, extension) is True

    def test_cleanup_global_extension(self, redis_client):
        """Test cleanup for global extensions (room_id=None)."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        category = "modifiers"
        extension = "GlobalExtension"

        keys = ExtensionKeys.for_global_extension(category, extension)
        stale_worker = "stale-global-worker"

        redis_client.hset(keys.workers, stale_worker, heartbeat_ts - 100)
        redis_client.hset(keys.schema, extension, '{"type": "object"}')

        mock_socketio = MagicMock()

        cleaned = cleanup_stale_extension_workers(
            redis_client, mock_socketio, None, category, extension
        )

        assert stale_worker in cleaned
        assert redis_client.hexists(keys.workers, stale_worker) is False


class TestStaleFilesystemCleanup:
    """Tests for stale filesystem worker cleanup."""

    def test_cleanup_removes_stale_filesystem(self, redis_client):
        """Test that cleanup_stale_filesystem_worker removes stale filesystem."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        fs_name = "test-fs"

        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        worker_id = "stale-fs-worker"
        stale_ts = heartbeat_ts - 100

        # Set up filesystem metadata with registration timestamp
        redis_client.hset(
            keys.metadata,
            mapping={
                "name": fs_name,
                "fsType": "local",
                "registration_timestamp": str(stale_ts),
            },
        )
        redis_client.set(keys.worker, worker_id)

        mock_socketio = MagicMock()

        cleaned_worker = cleanup_stale_filesystem_worker(
            redis_client, mock_socketio, room_id, fs_name
        )

        assert cleaned_worker == worker_id
        assert redis_client.exists(keys.metadata) == 0
        assert redis_client.get(keys.worker) is None

        mock_socketio.emit.assert_called()

    def test_cleanup_keeps_valid_filesystem(self, redis_client):
        """Test that cleanup_stale_filesystem_worker keeps valid filesystem."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        fs_name = "valid-fs"

        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        worker_id = "valid-fs-worker"
        valid_ts = heartbeat_ts + 10

        redis_client.hset(
            keys.metadata,
            mapping={
                "name": fs_name,
                "fsType": "local",
                "registration_timestamp": str(valid_ts),
            },
        )
        redis_client.set(keys.worker, worker_id)

        mock_socketio = MagicMock()

        cleaned_worker = cleanup_stale_filesystem_worker(
            redis_client, mock_socketio, room_id, fs_name
        )

        assert cleaned_worker is None
        assert redis_client.exists(keys.metadata) == 1
        assert redis_client.get(keys.worker) == worker_id

    def test_cleanup_global_filesystem(self, redis_client):
        """Test cleanup for global filesystems (room_id=None)."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        fs_name = "global-fs"

        keys = FilesystemKeys.for_global_filesystem(fs_name)
        worker_id = "stale-global-fs-worker"

        redis_client.hset(
            keys.metadata,
            mapping={
                "name": fs_name,
                "fsType": "s3",
                "registration_timestamp": str(heartbeat_ts - 100),
            },
        )
        redis_client.set(keys.worker, worker_id)

        mock_socketio = MagicMock()

        cleaned_worker = cleanup_stale_filesystem_worker(
            redis_client, mock_socketio, None, fs_name
        )

        assert cleaned_worker == worker_id
        assert redis_client.exists(keys.metadata) == 0

    def test_cleanup_skips_filesystem_without_timestamp(self, redis_client):
        """Test that filesystem without registration_timestamp is skipped."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        fs_name = "legacy-fs"

        keys = FilesystemKeys.for_filesystem(room_id, fs_name)
        worker_id = "legacy-worker"

        # No registration_timestamp in metadata
        redis_client.hset(
            keys.metadata,
            mapping={
                "name": fs_name,
                "fsType": "local",
            },
        )
        redis_client.set(keys.worker, worker_id)

        mock_socketio = MagicMock()

        cleaned_worker = cleanup_stale_filesystem_worker(
            redis_client, mock_socketio, room_id, fs_name
        )

        # Should not clean up (can't determine staleness)
        assert cleaned_worker is None
        assert redis_client.exists(keys.metadata) == 1


class TestCleanupAllStaleWorkersForCategory:
    """Tests for cleanup_all_stale_workers_for_category."""

    def test_cleanup_multiple_extensions(self, redis_client):
        """Test cleaning up stale workers across multiple extensions."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"

        # Set up two extensions with stale workers
        for ext_name in ["Extension1", "Extension2"]:
            keys = ExtensionKeys.for_extension(room_id, category, ext_name)
            redis_client.hset(keys.workers, f"stale-{ext_name}", heartbeat_ts - 100)
            redis_client.hset(keys.schema, ext_name, '{"type": "object"}')

        mock_socketio = MagicMock()

        cleaned = cleanup_all_stale_workers_for_category(
            redis_client, mock_socketio, room_id, category
        )

        assert "Extension1" in cleaned
        assert "Extension2" in cleaned
        assert "stale-Extension1" in cleaned["Extension1"]
        assert "stale-Extension2" in cleaned["Extension2"]


class TestWorkerStatsFiltersStaleWorkers:
    """Tests that WorkerStats correctly filters out stale workers."""

    def test_worker_stats_excludes_stale_workers(self, redis_client):
        """Test that WorkerStats.fetch() excludes stale workers from counts."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add stale worker
        stale_worker = "stale-worker"
        redis_client.hset(keys.workers, stale_worker, heartbeat_ts - 100)
        capacity_key = ExtensionKeys.worker_capacity_key(stale_worker)
        redis_client.set(capacity_key, 1)  # Idle

        # Add valid worker
        valid_worker = "valid-worker"
        redis_client.hset(keys.workers, valid_worker, heartbeat_ts + 10)
        capacity_key = ExtensionKeys.worker_capacity_key(valid_worker)
        redis_client.set(capacity_key, 1)  # Idle

        stats = WorkerStats.fetch(redis_client, keys)

        # Only valid worker should be counted
        assert stats.idle_count == 1
        assert stats.total_workers == 1

    def test_worker_stats_counts_only_valid_workers(self, redis_client):
        """Test WorkerStats counts with mix of idle and busy valid workers."""
        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add valid idle worker
        idle_worker = "idle-worker"
        redis_client.hset(keys.workers, idle_worker, heartbeat_ts + 10)
        redis_client.set(ExtensionKeys.worker_capacity_key(idle_worker), 1)

        # Add valid busy worker
        busy_worker = "busy-worker"
        redis_client.hset(keys.workers, busy_worker, heartbeat_ts + 10)
        redis_client.set(ExtensionKeys.worker_capacity_key(busy_worker), 0)

        # Add stale worker (should be ignored)
        stale_worker = "stale-worker"
        redis_client.hset(keys.workers, stale_worker, heartbeat_ts - 100)
        redis_client.set(ExtensionKeys.worker_capacity_key(stale_worker), 1)

        stats = WorkerStats.fetch(redis_client, keys)

        assert stats.idle_count == 1
        assert stats.busy_count == 1
        assert stats.total_workers == 2


class TestJobDispatcherFiltersStaleWorkers:
    """Tests that job dispatcher correctly filters out stale workers."""

    def test_get_available_workers_excludes_stale(self, redis_client):
        """Test that get_available_workers excludes stale workers."""
        from zndraw.app.job_dispatcher import get_available_workers

        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add stale worker
        stale_worker = "stale-worker"
        redis_client.hset(keys.workers, stale_worker, heartbeat_ts - 100)
        redis_client.set(ExtensionKeys.worker_capacity_key(stale_worker), 1)

        # Add valid worker
        valid_worker = "valid-worker"
        redis_client.hset(keys.workers, valid_worker, heartbeat_ts + 10)
        redis_client.set(ExtensionKeys.worker_capacity_key(valid_worker), 1)

        available = get_available_workers(redis_client, keys)

        assert valid_worker in available
        assert stale_worker not in available

    def test_get_available_workers_returns_empty_if_all_stale(self, redis_client):
        """Test that get_available_workers returns empty if all workers are stale."""
        from zndraw.app.job_dispatcher import get_available_workers

        heartbeat_ts = time.time()
        redis_client.set(ClusterKeys.HEARTBEAT, heartbeat_ts)

        room_id = "test-room"
        category = "modifiers"
        extension = "TestExtension"

        keys = ExtensionKeys.for_extension(room_id, category, extension)

        # Add only stale workers
        for i in range(3):
            worker = f"stale-worker-{i}"
            redis_client.hset(keys.workers, worker, heartbeat_ts - 100)
            redis_client.set(ExtensionKeys.worker_capacity_key(worker), 1)

        available = get_available_workers(redis_client, keys)

        assert available == []
