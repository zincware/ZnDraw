"""ZnDraw CLI - Command-line interface for ZnDraw visualization server.

This module provides commands to:
- Start a ZnDraw server (local or detached)
- Upload trajectory files to a running server
- Check server status
- Shutdown a running server
- Connect to a remote server
"""

from __future__ import annotations

import os
import re
import secrets
import sys
import threading
import uuid
import webbrowser
from pathlib import Path
from typing import Annotated

import typer
import uvicorn

from zndraw import __version__
from zndraw.client import ZnDraw
from zndraw.server_manager import (
    DEFAULT_PORT,
    ServerInfo,
    find_running_server,
    remove_server_info,
    shutdown_server,
    wait_for_server_ready,
    write_server_info,
)

app = typer.Typer(
    name="zndraw",
    help="ZnDraw - Interactive visualization for atomistic simulations",
    no_args_is_help=False,
)

# Separate CLI for database management (registered as `zndraw-db` entry point)
db_app = typer.Typer(
    name="zndraw-db",
    help="ZnDraw database management utilities",
)


@db_app.command()
def db_init(
    database_url: Annotated[
        str | None,
        typer.Option(help="Database URL (overrides ZNDRAW_DATABASE_URL)"),
    ] = None,
) -> None:
    """Initialize database tables.

    Run once before starting multiple workers in production:

        zndraw-db
        gunicorn -w 4 -k uvicorn.workers.UvicornWorker zndraw.app:socket_app

    For development (single worker), init happens automatically on startup.
    """
    import asyncio

    from zndraw.database import init_database

    if database_url:
        os.environ["ZNDRAW_DATABASE_URL"] = database_url

    asyncio.run(init_database())
    typer.echo("Database initialized successfully")


def sanitize_room_name(name: str) -> str:
    """Convert a string to a valid room name by replacing non-alphanumeric characters.

    Preserves letters, numbers, hyphens, and underscores. Other characters are
    replaced with underscores.
    """
    return re.sub(r"[^a-zA-Z0-9\-_]", "_", name)


def path_to_room(path: str, unique: bool = True) -> str:
    """Convert a file path to a valid room name.

    Parameters
    ----------
    path : str
        The file path to convert.
    unique : bool
        If True, append a random UUID suffix to ensure uniqueness.
        If False, return the sanitized path directly (for append mode).

    Returns
    -------
    str
        A valid room name, optionally with UUID suffix.
    """
    room = sanitize_room_name(path)
    if unique:
        room = f"{room}_{uuid.uuid4().hex[:4]}"
    return room


def daemonize(log_file: str = "zndraw.log") -> None:
    """Daemonize the current process by forking and detaching from terminal.

    Parameters
    ----------
    log_file : str
        Path to log file for stdout/stderr redirection.
    """
    try:
        pid = os.fork()
        if pid > 0:
            # Parent process - exit
            typer.echo(f"Server started in background (PID: {pid})")
            typer.echo(f"  Logs will be written to: {log_file}")
            typer.echo("  Use 'zndraw --status' to check server status")
            typer.echo("  Use 'zndraw --shutdown' to stop the server")
            sys.exit(0)
    except OSError as e:
        typer.echo(f"Fork failed: {e}", err=True)
        sys.exit(1)

    # Child process continues here
    # Create a new session and become session leader
    os.setsid()

    # Redirect standard file descriptors
    sys.stdout.flush()
    sys.stderr.flush()

    # Open log file for writing
    log_fd = os.open(log_file, os.O_CREAT | os.O_WRONLY | os.O_APPEND, 0o644)

    # Redirect stdout and stderr to log file
    os.dup2(log_fd, sys.stdout.fileno())
    os.dup2(log_fd, sys.stderr.fileno())

    # Close stdin
    devnull = os.open(os.devnull, os.O_RDONLY)
    os.dup2(devnull, sys.stdin.fileno())

    # Close the temporary file descriptors
    os.close(log_fd)
    os.close(devnull)


def upload_file(
    file_path: str,
    server_url: str,
    room_id: str,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
) -> None:
    """Upload a file to the server.

    Streams frames from disk via ``ase.io.iread`` directly into
    ``client.extend()``, so only one chunk (~2 MB) is in memory at a time.

    Parameters
    ----------
    file_path : str
        Path to the file to upload.
    server_url : str
        URL of the ZnDraw server.
    room_id : str
        Room ID to upload to.
    start, stop, step : int | None
        Slice parameters for frame selection.
    """
    import ase.io

    index: str | slice = (
        slice(start, stop, step)
        if any(x is not None for x in [start, stop, step])
        else ":"
    )

    try:
        frames = ase.io.iread(file_path, index=index)  # type: ignore[arg-type]
    except Exception as e:
        typer.echo(f"Error reading file {file_path}: {e}", err=True)
        raise typer.Exit(1) from e

    client = ZnDraw(url=server_url, room=room_id)
    try:
        old_len = len(client)
        client.extend(frames)  # type: ignore[arg-type]
        uploaded = len(client) - old_len
        if uploaded == 0:
            typer.echo(f"No frames found in {file_path}", err=True)
            raise typer.Exit(1)
        typer.echo(f"  Uploaded {uploaded} frames to room {room_id}")
    finally:
        client.disconnect()


# ── Extracted pipeline helpers ───────────────────────────────────────


def validate_flags(
    detached: bool,
    status: bool,
    shutdown: bool,
    connect: str | None,
    append: bool,
    room: str | None,
    path: list[str] | None,
) -> None:
    """Validate CLI flag combinations."""
    if detached and (status or shutdown or connect):
        typer.echo(
            "Error: --detached cannot be used with --status, --shutdown, or --connect",
            err=True,
        )
        raise typer.Exit(1)

    if append and room:
        typer.echo(
            "Error: --append and --room cannot be used together. "
            "Use --append for file-derived room names or --room for explicit naming.",
            err=True,
        )
        raise typer.Exit(1)

    if (append or room) and not path:
        typer.echo(
            "Error: --append and --room require file path(s) to be specified.",
            err=True,
        )
        raise typer.Exit(1)


def validate_files_exist(paths: list[str]) -> None:
    """Validate that all specified files exist."""
    for p in paths:
        if not Path(p).exists():
            typer.echo(f"Error: File not found: {p}", err=True)
            raise typer.Exit(1)


def handle_status(port: int | None) -> None:
    """Handle --status flag. Always raises typer.Exit."""
    server_info = find_running_server(port)

    if server_info is not None:
        typer.echo(
            f"Server running (PID: {server_info.pid}, "
            f"Port: {server_info.port}, Version: {server_info.version})"
        )
        typer.echo(f"  Server URL: http://localhost:{server_info.port}")
        raise typer.Exit(0)

    if port is not None:
        typer.echo(f"No ZnDraw server running on port {port}")
    else:
        typer.echo("No local ZnDraw server is running")
    raise typer.Exit(1)


def handle_shutdown(port: int | None) -> None:
    """Handle --shutdown flag. Always raises typer.Exit."""
    server_info = find_running_server(port)

    if server_info is None:
        if port is not None:
            typer.echo(f"No server running on port {port}. Nothing to shut down.")
        else:
            typer.echo("No running server found. Nothing to shut down.")
        raise typer.Exit(0)

    typer.echo(
        f"Shutting down server (PID: {server_info.pid}, Port: {server_info.port})..."
    )
    if shutdown_server(server_info):
        typer.echo("Server shut down successfully")
        raise typer.Exit(0)

    typer.echo("Failed to shut down server")
    raise typer.Exit(1)


def get_room_names(
    paths: list[str],
    room: str | None,
    append: bool,
) -> list[str]:
    """Compute room names for given paths based on --room and --append flags."""
    if room:
        return [sanitize_room_name(room)] * len(paths)
    if append:
        return [path_to_room(p, unique=False) for p in paths]
    return [path_to_room(p, unique=True) for p in paths]


def resolve_server(
    connect: str | None,
    port: int | None,
    host: str,
    detached: bool,
    verbose: bool,
) -> tuple[str, uvicorn.Server | None, int]:
    """Resolve or start a ZnDraw server.

    Returns
    -------
    tuple[str, uvicorn.Server | None, int]
        (server_url, server, effective_port).
        server is None for remote/existing connections.
        effective_port is only meaningful when server is not None.
    """
    if connect:
        typer.echo(f"Connecting to remote server: {connect}")
        return connect, None, 0

    server_info = find_running_server(port)
    if server_info is not None:
        typer.echo(
            f"Found existing server (PID: {server_info.pid}, "
            f"Port: {server_info.port}, Version: {server_info.version})"
        )
        typer.echo(f"  Server URL: http://localhost:{server_info.port}")

        if server_info.version != __version__:
            typer.echo(
                f"Warning: Server version ({server_info.version}) "
                f"differs from CLI version ({__version__})"
            )
            typer.echo(
                "  Consider running 'zndraw --shutdown' and starting a new server."
            )

        return f"http://localhost:{server_info.port}", None, server_info.port

    # Start new server
    effective_port = port if port is not None else DEFAULT_PORT

    # Write back so Settings() in lifespan reads the actual CLI values
    os.environ["ZNDRAW_HOST"] = host
    os.environ["ZNDRAW_PORT"] = str(effective_port)

    typer.echo(f"Starting new server on port {effective_port}...")

    shutdown_token = secrets.token_urlsafe(32)

    if detached:
        daemonize()

    write_server_info(
        ServerInfo(
            pid=os.getpid(),
            port=effective_port,
            version=__version__,
            shutdown_token=shutdown_token,
        )
    )

    log_level = "debug" if verbose else "info"

    # Import here to avoid circular imports
    from zndraw.app import app as fastapi_app, socket_app

    fastapi_app.state.shutdown_token = shutdown_token

    config = uvicorn.Config(
        socket_app,
        host=host,
        port=effective_port,
        log_level=log_level,
    )
    return f"http://{host}:{effective_port}", uvicorn.Server(config), effective_port


def open_browser_to(
    server_url: str,
    room: str | None,
    browser: bool,
    *,
    copy_from: str | None = None,
) -> None:
    """Open the browser to the appropriate server URL.

    Parameters
    ----------
    server_url : str
        Base server URL.
    room : str | None
        Room to navigate to, or None for root URL.
    browser : bool
        Whether to open the browser.
    copy_from : str | None
        If set, append ``?copy_from=<value>`` query parameter.
    """
    if not browser:
        return
    url = server_url
    if room:
        url = f"{url}/rooms/{room}"
        if copy_from:
            url += f"?copy_from={copy_from}"
    typer.echo(f"Opening browser at {url}")
    webbrowser.open(url)


def upload_files(
    paths: list[str],
    server_url: str,
    room_names: list[str],
    start: int | None,
    stop: int | None,
    step: int | None,
) -> None:
    """Upload multiple files to the server."""
    if not paths:
        return
    typer.echo("Uploading files...")
    for p, room_name in zip(paths, room_names, strict=True):
        typer.echo(f"  Uploading {p} to room {room_name}")
        upload_file(p, server_url, room_name, start, stop, step)
    typer.echo("Files uploaded successfully")


@app.command()
def main(
    path: Annotated[
        list[str] | None,
        typer.Argument(help="Path to file(s) to load on startup (optional)."),
    ] = None,
    start: Annotated[
        int | None,
        typer.Option(help="Start frame index for slicing."),
    ] = None,
    stop: Annotated[
        int | None,
        typer.Option(help="Stop frame index for slicing."),
    ] = None,
    step: Annotated[
        int | None,
        typer.Option(help="Step for frame slicing."),
    ] = None,
    append: Annotated[
        bool,
        typer.Option(
            "--append",
            help="Append to existing room derived from file path. "
            "Without this flag, a unique room name is generated each time.",
        ),
    ] = False,
    room: Annotated[
        str | None,
        typer.Option(
            "--room",
            help="Explicitly specify the room name. "
            "All files will be loaded into this room.",
        ),
    ] = None,
    port: Annotated[
        int | None,
        typer.Option(
            "--port",
            help="Server port. If specified and server exists on "
            "that port, connects to it. "
            "If specified and no server on that port, starts new server. "
            "If not specified, auto-discovers running servers.",
            envvar="ZNDRAW_PORT",
        ),
    ] = None,
    host: Annotated[
        str,
        typer.Option(help="Server hostname or IP address.", envvar="ZNDRAW_HOST"),
    ] = "127.0.0.1",
    connect: Annotated[
        str | None,
        typer.Option(
            "--connect",
            help="Connect to a remote server URL, bypassing local server checks.",
        ),
    ] = None,
    status: Annotated[
        bool,
        typer.Option(
            "--status",
            help="Check if a local server is running.",
        ),
    ] = False,
    shutdown: Annotated[
        bool,
        typer.Option(
            "--shutdown",
            help="Stop the local server.",
        ),
    ] = False,
    browser: Annotated[
        bool,
        typer.Option(
            "--browser/--no-browser",
            help="Automatically open the web browser.",
        ),
    ] = True,
    detached: Annotated[
        bool,
        typer.Option(
            "--detached",
            help="Start the server as a detached background process.",
        ),
    ] = False,
    verbose: Annotated[
        bool,
        typer.Option("--verbose", "-v", help="Enable verbose logging."),
    ] = False,
    version: Annotated[
        bool,
        typer.Option("--version", "-V", help="Show version and exit."),
    ] = False,
) -> None:
    """Start or connect to a ZnDraw server.

    By default, this command will check if a local server is already running and
    connect to it. If no server is running, it will start a new one.

    Examples:

        # Start a new server
        zndraw

        # Load a file into the server
        zndraw trajectory.xyz

        # Load multiple files
        zndraw file1.xyz file2.xyz

        # Load specific frames
        zndraw trajectory.xyz --start 0 --stop 100 --step 2

        # Connect to remote server
        zndraw --connect https://zndraw.example.com trajectory.xyz

        # Check server status
        zndraw --status

        # Shutdown running server
        zndraw --shutdown
    """
    if version:
        typer.echo(f"zndraw {__version__}")
        raise typer.Exit(0)

    # ── Validate ─────────────────────────────────────────────────────
    validate_flags(detached, status, shutdown, connect, append, room, path)

    if status:
        handle_status(port)
    if shutdown:
        handle_shutdown(port)
    if path:
        validate_files_exist(path)

    if detached and browser:
        typer.echo("Note: Browser will not be opened in detached mode")
        browser = False

    # ── Compute rooms ────────────────────────────────────────────────
    room_names = get_room_names(path or [], room, append)
    first_room = room_names[0] if room_names else f"workspace-{uuid.uuid4().hex[:8]}"
    has_files = bool(path)

    if verbose:
        msg = f"Rooms: {room_names}" if has_files else "No files loaded on startup."
        typer.echo(msg)

    # ── Resolve server ───────────────────────────────────────────────
    url, server, effective_port = resolve_server(connect, port, host, detached, verbose)

    if server is None:
        # Remote or existing server — open browser, upload, done
        open_browser_to(url, first_room if has_files else None, browser)
        upload_files(path or [], url, room_names, start, stop, step)
        if not connect:
            typer.echo(f"\nServer is running at {url}")
        return

    # ── New server ───────────────────────────────────────────────────
    if not (has_files or browser):
        # Nothing to do post-startup — run server blocking
        typer.echo(f"Server starting at {url}")
        try:
            server.run()
        except KeyboardInterrupt:
            typer.echo("\nShutting down...")
        remove_server_info(effective_port)
        typer.echo("Server stopped.")
        return

    # Thread server so we can open browser / upload before it exits
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()

    if not wait_for_server_ready(url, timeout=30.0):
        typer.echo("Error: Server failed to start within timeout", err=True)
        raise typer.Exit(1)
    typer.echo(f"Server started at {url}")

    open_browser_to(url, first_room, browser, copy_from="@none" if has_files else None)
    upload_files(path or [], url, room_names, start, stop, step)

    try:
        thread.join()
    except KeyboardInterrupt:
        typer.echo("\nShutting down...")
    remove_server_info(effective_port)
    typer.echo("Server stopped.")


if __name__ == "__main__":
    app()
