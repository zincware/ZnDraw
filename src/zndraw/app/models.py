"""Pydantic models for room metadata shared between REST API and WebSocket events."""

from pydantic import BaseModel, Field


class LockMetadata(BaseModel):
    """Lock metadata for temporary locks (e.g., uploads).

    Attributes
    ----------
    msg : str | None
        Human-readable message describing the lock purpose
    userName : str | None
        Name of user holding the lock
    timestamp : str | None
        ISO 8601 timestamp when lock was acquired (e.g., "2025-01-24T19:27:16.369371")
    """

    msg: str | None = None
    userName: str | None = None
    timestamp: str | None = None


class RoomMetadata(BaseModel):
    """Room metadata shared between REST API and WebSocket events.

    This is the single source of truth for room state that is synchronized
    between clients via WebSocket events.

    Attributes
    ----------
    id : str
        Unique room identifier
    description : str | None
        Human-readable description
    frameCount : int
        Number of trajectory frames
    locked : bool
        Permanent lock (immutable)
    hidden : bool
        Hidden from room list
    isDefault : bool
        Set as default room
    presenterSid : str | None
        Socket ID of current presenter (for presenter mode)
    """

    id: str = Field(..., description="Unique room identifier")
    description: str | None = Field(None, description="Human-readable description")
    frameCount: int = Field(0, description="Number of trajectory frames")
    locked: bool = Field(False, description="Permanent lock (immutable)")
    hidden: bool = Field(False, description="Hidden from room list")
    isDefault: bool = Field(False, description="Set as default room")
    presenterSid: str | None = Field(None, description="Socket ID of current presenter")
