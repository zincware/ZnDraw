"""Redis key patterns for ephemeral state.

All persistent room state (geometries, bookmarks, figures, selections, step)
is stored in SQL. Redis is used only for:
- Session cameras (frontend presence derived from camera hash)
- Ephemeral locks (TTL-based with metadata)
- Socket.IO pub/sub adapter
"""


class RedisKey:
    """Redis key patterns - avoids magic strings."""

    # =========================================================================
    # Session Camera Keys
    # =========================================================================

    @staticmethod
    def room_cameras(room_id: str) -> str:
        """Hash key for all session cameras in a room.

        Each hash field = camera key (``cam:{email}:{uuid_short}``),
        hash value = JSON blob with sid, user_id, email, data.
        """
        return f"room:{room_id}:cameras"

    # =========================================================================
    # Active Camera Keys
    # =========================================================================

    @staticmethod
    def active_cameras(room_id: str) -> str:
        """Hash: session_id -> camera_key for all sessions in a room."""
        return f"room:{room_id}:active-cameras"

    # =========================================================================
    # Edit Lock Keys
    # =========================================================================

    @staticmethod
    def edit_lock(room_id: str) -> str:
        """Room-level edit lock key (JSON value with 10s TTL)."""
        return f"room:{room_id}:edit-lock"

    # =========================================================================
    # Download Token Keys
    # =========================================================================

    @staticmethod
    def download_token(token: str) -> str:
        """Key for a temporary download token. Value is the room_id."""
        return f"download-token:{token}"

    # =========================================================================
    # Provider Frame Count Keys
    # =========================================================================

    @staticmethod
    def provider_frame_count(room_id: str) -> str:
        """External frame count set by a provider mount (int value)."""
        return f"room:{room_id}:provider_frame_count"

    # =========================================================================
    # Provider Dispatch Keys
    # =========================================================================

    @staticmethod
    def provider_result(full_name: str, rhash: str) -> str:
        """Cache key for a completed provider result (msgpack bytes)."""
        return f"provider-result:{full_name}:{rhash}"

    @staticmethod
    def provider_inflight(full_name: str, rhash: str) -> str:
        """Inflight lock key preventing duplicate provider dispatches."""
        return f"provider-inflight:{full_name}:{rhash}"

    # =========================================================================
    # Progress Tracker Keys
    # =========================================================================

    @staticmethod
    def room_progress(room_id: str) -> str:
        """Hash of active progress trackers in a room.

        Field: progress_id, Value: JSON ProgressResponse.
        """
        return f"room:{room_id}:progress"
