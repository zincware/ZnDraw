"""Version compatibility utilities for ZnDraw client-server connections."""

import logging
import re
import warnings

log = logging.getLogger(__name__)


def parse_version(version: str) -> tuple[int, int, int] | None:
    """Parse a semantic version string into its components.

    Handles pre-release versions like "0.6.0a4"

    Parameters
    ----------
    version : str
        Version string to parse

    Returns
    -------
    tuple[int, int, int] | None
        (major, minor, patch) or None if parsing fails
    """
    # Strip any pre-release suffix (e.g., "a4", "b1", "rc2")
    clean_version = re.sub(r"[a-zA-Z].*", "", version)
    parts = clean_version.split(".")

    if len(parts) < 3:
        return None

    try:
        major, minor, patch = int(parts[0]), int(parts[1]), int(parts[2])
        return (major, minor, patch)
    except ValueError:
        return None


def check_version_compatibility(
    client_version: str, server_version: str
) -> tuple[bool, str, str]:
    """Check version compatibility between client and server.

    Parameters
    ----------
    client_version : str
        The client version string
    server_version : str
        The server version string

    Returns
    -------
    tuple[bool, str, str]
        (compatible, severity, message) where:
        - compatible: bool indicating if versions are compatible
        - severity: 'error' | 'warning' | 'info'
        - message: human-readable message
    """
    client = parse_version(client_version)
    server = parse_version(server_version)

    if not client or not server:
        return (
            False,
            "error",
            f"Invalid version format. Client: {client_version}, Server: {server_version}",
        )

    client_major, client_minor, client_patch = client
    server_major, server_minor, server_patch = server

    # Check major version mismatch
    if client_major != server_major:
        return (
            False,
            "error",
            f"Major version mismatch. Client: {client_version}, Server: {server_version}. Connection refused.",
        )

    # Check minor version mismatch
    if client_minor != server_minor:
        return (
            False,
            "error",
            f"Minor version mismatch. Client: {client_version}, Server: {server_version}. Connection refused.",
        )

    # Check patch version mismatch (warning only)
    if client_patch != server_patch:
        return (
            True,
            "warning",
            f"Patch version mismatch. Client: {client_version}, Server: {server_version}. This may cause unexpected behavior.",
        )

    return (True, "info", f"Version compatible: {client_version}")


def validate_server_version(api_manager, client_version: str) -> None:
    """Validate that the server version is compatible with the client.

    Parameters
    ----------
    api_manager : APIManager
        The API manager to use for fetching the server version
    client_version : str
        The client version string

    Raises
    ------
    RuntimeError
        If versions are incompatible (major or minor mismatch)

    Warnings
    --------
    UserWarning
        If patch versions differ
    """
    try:
        server_version = api_manager.get_version()

        log.info(f"Server version: {server_version}")
        log.info(f"Client version: {client_version}")

        compatible, severity, message = check_version_compatibility(
            client_version, server_version
        )

        if not compatible:
            log.error(message)
            raise RuntimeError(message)

        if severity == "warning":
            log.warning(message)
            warnings.warn(message, UserWarning, stacklevel=2)
        else:
            log.info(message)

    except Exception as e:
        if isinstance(e, RuntimeError):
            raise
        log.warning(f"Could not validate server version: {e}")
        # Don't fail connection if version check fails for other reasons
