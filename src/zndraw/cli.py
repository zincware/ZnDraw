import os
import secrets
import sys
import uuid
import webbrowser

import typer

from zndraw import __version__
from zndraw.app.tasks import read_file
from zndraw.config import (
    LMDBStorageConfig,
    MongoDBStorageConfig,
    ZnDrawConfig,
    set_config,
)
from zndraw.server import create_app, socketio
from zndraw.server_manager import (
    DEFAULT_PORT,
    ServerInfo,
    find_running_server,
    remove_server_info,
    shutdown_server,
    write_server_info,
)
from zndraw.start_celery import run_celery_worker
from zndraw.utils import path_to_room, sanitize_room_name

app = typer.Typer()


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
            typer.echo(f"✓ Server started in background (PID: {pid})")
            typer.echo(f"  Logs will be written to: {log_file}")
            typer.echo("  Use 'zndraw --status' to check server status")
            typer.echo("  Use 'zndraw --shutdown' to stop the server")
            sys.exit(0)
    except OSError as e:
        typer.echo(f"✗ Fork failed: {e}", err=True)
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


def _build_config(
    port: int | None,
    host: str | None,
    redis_url: str | None,
    storage_type: str | None,
    storage_path: str | None,
    storage_url: str | None,
    storage_database: str | None,
    media_path: str | None,
    simgen: bool | None,
    file_browser: bool | None,
    file_browser_root: str | None,
    celery: bool,
    verbose: bool,
) -> ZnDrawConfig:
    """Build config from CLI arguments.

    CLI arguments override environment variables. If a CLI argument is None,
    pydantic-settings loads the value from environment variables.
    """
    # Build kwargs for config - only include non-None values
    config_kwargs: dict = {}

    # Server settings - only set if explicitly provided
    if port is not None:
        config_kwargs["server_port"] = port
    if host is not None:
        config_kwargs["server_host"] = host

    # Redis
    if redis_url is not None:
        config_kwargs["redis_url"] = redis_url

    # Storage - only construct if any storage arg was provided
    if any([storage_type, storage_path, storage_url, storage_database]):
        # Determine storage type
        effective_type = storage_type or "lmdb"

        if effective_type not in ("lmdb", "mongodb"):
            raise typer.BadParameter(
                f"--storage-type must be 'lmdb' or 'mongodb', got '{effective_type}'"
            )

        if effective_type == "mongodb":
            if storage_url is None:
                raise typer.BadParameter(
                    "--storage-url is required when --storage-type=mongodb"
                )
            config_kwargs["storage"] = MongoDBStorageConfig(
                url=storage_url,
                database=storage_database or "zndraw",
            )
        else:
            config_kwargs["storage"] = LMDBStorageConfig(
                path=storage_path or "./zndraw-data",
            )

    # Media path
    if media_path is not None:
        config_kwargs["media_path"] = media_path

    # Feature flags - only set if explicitly provided via CLI
    # When None, pydantic-settings reads from env vars (e.g., ZNDRAW_SIMGEN_ENABLED)
    if simgen is not None:
        config_kwargs["simgen_enabled"] = simgen
    if file_browser is not None:
        config_kwargs["file_browser_enabled"] = file_browser

    if file_browser_root is not None:
        config_kwargs["file_browser_root"] = file_browser_root

    # Verbose mode sets DEBUG logging
    if verbose:
        config_kwargs["log_level"] = "DEBUG"

    return ZnDrawConfig(**config_kwargs)


@app.command()
def main(
    path: list[str] | None = typer.Argument(
        None, help="Path to file(s) to load on startup (optional)."
    ),
    start: int | None = typer.Option(
        None, help="Start frame (optional, only for certain file types)."
    ),
    stop: int | None = typer.Option(
        None, help="Stop frame (optional, only for certain file types)."
    ),
    step: int | None = typer.Option(
        None, help="Step frame (optional, only for certain file types)."
    ),
    append: bool = typer.Option(
        False,
        "--append",
        help="Append to existing room derived from file path (e.g., tmp/s22.xyz -> tmp_s22_xyz). "
        "Without this flag, a unique room name is generated each time.",
    ),
    room: str | None = typer.Option(
        None,
        "--room",
        help="Explicitly specify the room name. All files will be loaded into this room.",
    ),
    port: int | None = typer.Option(
        None,
        "--port",
        help="Server port. If specified and server exists on that port, connects to it. "
        "If specified and no server on that port, starts new server. "
        "If not specified, auto-discovers running servers (default port first, then smallest).",
    ),
    debug: bool = typer.Option(False, help="Enable debug mode."),
    verbose: bool = typer.Option(False, help="Enable verbose logging."),
    celery: bool = typer.Option(True, help="Enable Celery task processing."),
    storage_type: str | None = typer.Option(
        None,
        "--storage-type",
        help="Storage backend type: 'lmdb' (local) or 'mongodb' (distributed).",
    ),
    storage_path: str | None = typer.Option(
        None,
        "--storage-path",
        help="Path to storage directory for LMDB files (when storage-type=lmdb).",
    ),
    storage_url: str | None = typer.Option(
        None,
        "--storage-url",
        help="MongoDB connection URL (when storage-type=mongodb).",
    ),
    storage_database: str | None = typer.Option(
        None,
        "--storage-database",
        help="MongoDB database name (when storage-type=mongodb).",
    ),
    media_path: str | None = typer.Option(
        None,
        "--media-path",
        help="Path for local media files (screenshots, etc.).",
    ),
    redis_url: str | None = typer.Option(
        None,
        help="Redis server URL (e.g., `redis://localhost:6379`). "
        "If not provided, uses in-memory storage.",
    ),
    host: str | None = typer.Option(
        None,
        help="Server hostname or IP address (default: localhost).",
    ),
    connect: str | None = typer.Option(
        None,
        "--connect",
        help="Connect to a specified remote server (e.g., https://zndraw.myserver.com), "
        "bypassing local server checks.",
    ),
    status: bool = typer.Option(
        False,
        "--status",
        help="Check if a local server is running. Uses --port if specified, otherwise auto-discovers.",
    ),
    shutdown: bool = typer.Option(
        False,
        "--shutdown",
        help="Stop the local server. Uses --port if specified, otherwise auto-discovers.",
    ),
    browser: bool = typer.Option(
        True,
        "--browser/--no-browser",
        help="Automatically open the web browser.",
    ),
    file_browser: bool | None = typer.Option(
        None,
        "--file-browser/--no-file-browser",
        help="Enable local filesystem browser endpoint (default: disabled).",
    ),
    file_browser_root: str | None = typer.Option(
        None,
        "--file-browser-root",
        help="Root directory for file browser (defaults to current working directory).",
    ),
    detached: bool = typer.Option(
        False,
        "--detached",
        help="Start the server as a detached background process.",
    ),
    simgen: bool | None = typer.Option(
        None,
        "--simgen/--no-simgen",
        help="Enable SiMGen molecular generation features (default: disabled).",
    ),
):
    """Start or connect to a ZnDraw server.

    By default, this command will check if a local server is already running and
    connect to it. If no server is running, it will start a new one on the default
    port (5000).

    If --port is specified:
    - If a server is running on that port, connects to it
    - If no server on that port, starts a new server on that port

    Use --status to check server status or --shutdown to stop the server.
    """
    # Validate flag combinations
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
    # Handle --status flag
    if status:
        server_info = find_running_server(port)

        if server_info is not None:
            typer.echo(
                f"✓ Server running (PID: {server_info.pid}, "
                f"Port: {server_info.port}, Version: {server_info.version})"
            )
            typer.echo(f"  Server URL: http://localhost:{server_info.port}")
            raise typer.Exit(0)
        else:
            if port is not None:
                typer.echo(f"✗ No ZnDraw server running on port {port}")
            else:
                typer.echo("✗ No local ZnDraw server is running")
            raise typer.Exit(1)

    # Handle --shutdown flag
    if shutdown:
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
            typer.echo("✓ Server shut down successfully")
            raise typer.Exit(0)
        else:
            typer.echo("✗ Failed to shut down server")
            raise typer.Exit(1)

    # Helper to compute room names based on flags
    def get_room_names(paths: list[str]) -> list[str]:
        """Compute room names for given paths based on --room and --append flags."""
        if room:
            # All files go to the same explicitly named room
            return [sanitize_room_name(room)] * len(paths)
        elif append:
            # Each file goes to a deterministic room derived from its path
            return [path_to_room(p, unique=False) for p in paths]
        else:
            # Default: each file gets a unique room name
            return [path_to_room(p, unique=True) for p in paths]

    # Handle remote connection
    if connect:
        typer.echo(f"Connecting to remote server: {connect}")

        if path is not None:
            # Validate all files exist before proceeding
            for p in path:
                if not os.path.exists(p):
                    typer.echo(f"✗ Error: File not found: {p}", err=True)
                    raise typer.Exit(1)

            typer.echo("Uploading files to remote server...")
            make_default = True
            room_names = get_room_names(path)
            first_room = room_names[0]

            # Open browser before uploading so user can watch progress
            if browser:
                browser_url = f"{connect}/rooms/{first_room}"
                typer.echo(f"Opening browser at {browser_url}")
                webbrowser.open(browser_url)

            for p, room_name in zip(path, room_names):
                typer.echo(f"  Uploading file {p} to room {room_name}")
                read_file(
                    file=p,
                    room=room_name,
                    server_url=connect,
                    start=start,
                    stop=stop,
                    step=step,
                    make_default=make_default,
                )
                make_default = False
            typer.echo("✓ Files uploaded successfully")

        else:
            typer.echo(f"Connected to remote server at {connect}")
            typer.echo("No files to upload.")

            # Open browser to root if no files
            if browser:
                typer.echo(f"Opening browser at {connect}")
                webbrowser.open(connect)

        raise typer.Exit(0)

    # Check for existing server
    # If port is specified, only check that port
    # If port is None, auto-discover (default port first, then smallest)
    server_info = find_running_server(port)

    if server_info is not None:
        typer.echo(
            f"✓ Found existing server (PID: {server_info.pid}, "
            f"Port: {server_info.port}, Version: {server_info.version})"
        )
        typer.echo(f"  Server URL: http://localhost:{server_info.port}")

        # Check version compatibility
        if server_info.version != __version__:
            typer.echo(
                f"⚠ Warning: Server version ({server_info.version}) "
                f"differs from CLI version ({__version__})"
            )
            typer.echo(
                "  Consider running 'zndraw --shutdown' and starting a new server."
            )

        # If files are provided, load them into the existing server
        if path is not None:
            # Validate all files exist before proceeding
            for p in path:
                if not os.path.exists(p):
                    typer.echo(f"✗ Error: File not found: {p}", err=True)
                    raise typer.Exit(1)

            typer.echo("Uploading files to existing server...")
            make_default = True
            server_url = f"http://localhost:{server_info.port}"
            room_names = get_room_names(path)
            first_room = room_names[0]

            # Open browser before uploading so user can watch progress
            if browser:
                browser_url = f"http://localhost:{server_info.port}/rooms/{first_room}"
                typer.echo(f"Opening browser at {browser_url}")
                webbrowser.open(browser_url)

            for p, room_name in zip(path, room_names):
                typer.echo(f"  Uploading file {p} to room {room_name}")
                read_file(
                    file=p,
                    room=room_name,
                    server_url=server_url,
                    start=start,
                    stop=stop,
                    step=step,
                    make_default=make_default,
                )
                make_default = False
            typer.echo("✓ Files uploaded successfully")

        else:
            # Open browser to root if no files
            if browser:
                browser_url = f"http://localhost:{server_info.port}"
                typer.echo(f"Opening browser at {browser_url}")
                webbrowser.open(browser_url)

        typer.echo(f"\nServer is running at http://localhost:{server_info.port}")
        raise typer.Exit(0)

    # No running server found - start a new one
    # Validate files exist before starting server
    if path is not None:
        for p in path:
            if not os.path.exists(p):
                typer.echo(f"✗ Error: File not found: {p}", err=True)
                raise typer.Exit(1)

    # Determine the port to use
    effective_port = port if port is not None else DEFAULT_PORT
    typer.echo(f"Starting new server on port {effective_port}...")

    # Compute room names upfront (if files provided)
    room_names = get_room_names(path) if path else []
    if verbose:
        typer.echo(f"Rooms: {room_names}" if path else "No files loaded on startup.")

    # Build configuration from CLI args (overrides env vars via pydantic-settings)
    config = _build_config(
        port=effective_port,
        host=host,
        redis_url=redis_url,
        storage_type=storage_type,
        storage_path=storage_path,
        storage_url=storage_url,
        storage_database=storage_database,
        media_path=media_path,
        simgen=simgen,
        file_browser=file_browser,
        file_browser_root=file_browser_root,
        celery=celery,
        verbose=verbose,
    )
    set_config(config)

    flask_app = create_app(config=config)

    # Track the first room for browser opening
    first_room = room_names[0] if room_names else None

    # Daemonize if requested (must happen before starting workers and writing PID file)
    if detached:
        # Browser opening is not compatible with detached mode
        if browser:
            typer.echo("Note: Browser will not be opened in detached mode")
            browser = False
        daemonize()

    # Start celery worker after daemonization so its logs go to the log file
    if celery:
        worker = run_celery_worker(config)

    # Generate secure shutdown token for CLI-based shutdown
    shutdown_token = secrets.token_urlsafe(32)

    # Store shutdown token in Flask app config for validation
    flask_app.config["SHUTDOWN_TOKEN"] = shutdown_token

    # Write server info to PID file
    new_server_info = ServerInfo(
        pid=os.getpid(),
        port=config.server_port,
        version=__version__,
        shutdown_token=shutdown_token,
    )
    write_server_info(new_server_info)
    typer.echo(
        f"✓ Server started (PID: {new_server_info.pid}, Port: {config.server_port})"
    )
    typer.echo(f"  Server URL: {config.server_url}")

    # Queue file loading tasks after worker is started
    if path is not None:
        make_default = True
        for p, room_name in zip(path, room_names):
            if verbose:
                typer.echo(f"Loading file {p} into room {room_name}.")
            read_file.delay(
                file=p,
                room=room_name,
                server_url=config.server_url,
                start=start,
                stop=stop,
                step=step,
                make_default=make_default,
            )
            make_default = False
    else:
        # No file provided - just set room name for browser URL
        # Room will be auto-created with "empty" template when frontend joins
        workspace_room = f"workspace-{uuid.uuid4().hex[:8]}"
        first_room = workspace_room
        if verbose:
            typer.echo(f"Starting empty workspace: {workspace_room}")

    # Open browser if requested
    if browser:
        if first_room:
            browser_url = f"http://localhost:{config.server_port}/rooms/{first_room}"
        else:
            browser_url = f"http://localhost:{config.server_port}"
        if verbose:
            typer.echo(f"Opening browser at {browser_url}")
        webbrowser.open(browser_url)

    try:
        socketio.run(
            flask_app, debug=debug, host=config.server_host, port=config.server_port
        )
    finally:
        flask_app.extensions["redis"].flushall()
        if celery:
            # Use SIGKILL directly instead of SIGTERM. Celery with eventlet pool
            # has a bug where the SIGTERM handler tries to print a shutdown message
            # using eventlet's patched os.write(), which raises RuntimeError
            # "do not call blocking functions from the mainloop".
            worker.kill()
            worker.wait()
            typer.echo("Celery worker closed.")
        # Clean up PID file when server stops
        remove_server_info(config.server_port)
        typer.echo("Server stopped.")
