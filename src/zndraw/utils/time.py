"""Time utilities for consistent UTC timestamp handling.

Use these functions instead of deprecated datetime.utcnow().
"""

import datetime


def utc_now() -> datetime.datetime:
    """Get current UTC time (timezone-aware).

    Use this instead of deprecated datetime.utcnow().

    Returns
    -------
    datetime.datetime
        Current UTC time with timezone info
    """
    return datetime.datetime.now(datetime.timezone.utc)


def utc_now_iso() -> str:
    """Get current UTC time as ISO format string.

    Returns
    -------
    str
        Current UTC time in ISO 8601 format
    """
    return utc_now().isoformat()


def utc_now_timestamp() -> float:
    """Get current UTC time as Unix timestamp.

    Returns
    -------
    float
        Current UTC time as seconds since epoch
    """
    return utc_now().timestamp()
