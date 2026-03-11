"""Socket.IO manager for the ZnDraw client."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import socketio
from zndraw_socketio import SyncClientWrapper, wrap

from zndraw.client.exceptions import NotConnectedError
from zndraw.exceptions import ZnDrawError

if TYPE_CHECKING:
    from zndraw.client.core import ZnDraw

log = logging.getLogger(__name__)


# =============================================================================
# SocketManager - WebSocket Real-time Sync
# =============================================================================


@dataclass
class SocketManager:
    """Manages Socket.IO connection and real-time synchronization."""

    zndraw: ZnDraw
    tsio: SyncClientWrapper = field(init=False)
    _connected: bool = field(default=False, init=False)

    def __post_init__(self) -> None:
        """Initialize wrapped Socket.IO client and register handlers."""
        self.tsio = wrap(socketio.Client())
        self._register_handlers()

    def _register_handlers(self) -> None:
        """Register Socket.IO event handlers using typed models."""
        from zndraw.socket_events import FramesInvalidate

        self.tsio.on("connect", self._on_connect)
        self.tsio.on("disconnect", self._on_disconnect)
        self.tsio.on(FramesInvalidate, self._on_frames_invalidate)

    @property
    def connected(self) -> bool:
        """Check if connected."""
        return self._connected and self.tsio.connected

    def connect(self) -> None:
        """Connect to the server and join the room."""
        from zndraw.socket_events import RoomJoin, RoomJoinResponse

        if self.tsio.connected:
            log.debug("Already connected")
            return

        # Connect with JWT token
        self.tsio.connect(
            self.zndraw.url,
            auth={"token": self.zndraw.api.token},
            wait=True,
        )

        # Join room and get session info
        if self.zndraw.room is None:
            raise NotConnectedError("Cannot join: no room set")
        join_request = RoomJoin(room_id=self.zndraw.room, client_type="pyclient")
        raw_response = self.tsio.call(join_request)
        response: dict[str, Any] = raw_response if raw_response else {}

        if "type" in response and response.get("status") == 404:
            self.zndraw.api.create_room(copy_from=self.zndraw.copy_from)
            # Retry join
            raw_response = self.tsio.call(join_request)
            response = raw_response if raw_response else {}

        if "type" in response:
            error_msg = (
                response.get("detail") or response.get("title") or "Unknown error"
            )
            raise ZnDrawError(f"Failed to join room: {error_msg}")

        # Validate response
        join_response = RoomJoinResponse.model_validate(response)
        self.zndraw.api.session_id = join_response.session_id

        # Seed the frame count cache from join response
        self.zndraw.cached_length = join_response.frame_count

        self._connected = True
        log.debug("Connected to room %s", self.zndraw.room)

    def disconnect(self) -> None:
        """Disconnect from the server."""
        if self.tsio.connected:
            self.tsio.disconnect()
        self._connected = False
        log.debug("Disconnected")

    def wait(self) -> None:
        """Block until disconnected."""
        self.tsio.wait()

    def _on_connect(self) -> None:
        """Handle connection event. Re-register providers if mounted."""
        log.debug("Socket connected")
        if self.zndraw.mount is not None:
            from zndraw.providers.frame_source import FrameSourceLength, FrameSourceRead

            source = self.zndraw.mount
            room = self.zndraw.room
            if room is None:
                return
            mount_name = self.zndraw.mount_name
            if mount_name is None:
                return
            self.zndraw.register_provider(
                FrameSourceRead, name=mount_name, handler=source, room=room
            )
            self.zndraw.register_provider(
                FrameSourceLength, name=mount_name, handler=source, room=room
            )
            self.zndraw.api.update_room({"frame_count": len(source)})

    def _on_disconnect(self, reason: str = "") -> None:
        """Handle disconnection event."""
        log.debug("Socket disconnected: %s", reason)
        self._connected = False

    def _on_frames_invalidate(self, data: Any) -> None:
        """Handle FramesInvalidate broadcast -- update cached length."""
        from zndraw.socket_events import FramesInvalidate

        event = FramesInvalidate.model_validate(data)
        if event.room_id != self.zndraw.room:
            return
        if event.count is not None:
            self.zndraw.cached_length = event.count
        else:
            self.zndraw.cached_length = None
