"""Redis key patterns for presence tracking and ephemeral locks.

All persistent room state (geometries, bookmarks, figures, selections, step)
is stored in SQL. Redis is used only for:
- Presence tracking (TTL-based per-session keys)
- Ephemeral locks (TTL-based with metadata)
- Socket.IO pub/sub adapter
"""


class RedisKey:
    """Redis key patterns - avoids magic strings."""

    # =========================================================================
    # Presence Keys
    # =========================================================================

    @staticmethod
    def presence_sid(room_id: str, sid: str) -> str:
        """Key for session presence in a room (per-sid).

        Each Socket.IO session (tab) gets its own presence key.
        This enables multi-tab support for the same user.
        """
        return f"presence:room:{room_id}:sid:{sid}"

    @staticmethod
    def presence_sid_pattern(room_id: str) -> str:
        """Pattern for scanning all sessions in a room."""
        return f"presence:room:{room_id}:sid:*"

    @staticmethod
    def parse_presence_sid(key: str) -> str | None:
        """Extract SID from a presence key.

        Expected format: presence:room:{room_id}:sid:{sid}
        Returns the SID or None if the key format is invalid.
        """
        parts = key.split(":")
        if (
            len(parts) == 5
            and parts[0] == "presence"
            and parts[1] == "room"
            and parts[3] == "sid"
        ):
            return parts[4]
        return None

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
    # Session Settings Keys
    # =========================================================================

    @staticmethod
    def session_settings(room_id: str, sid: str) -> str:
        """Per-session settings key (JSON-serialized RoomConfig)."""
        return f"room:{room_id}:session-settings:{sid}"

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
