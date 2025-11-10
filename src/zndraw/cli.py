import os
import re
import shutil
import subprocess
import sys
import webbrowser

import typer

from zndraw import __version__
from zndraw.app.tasks import read_file
from zndraw.server import create_app, socketio
from zndraw.server_manager import (
    ServerInfo,
    get_server_status,
    is_process_running,
    read_server_info,
    remove_server_info,
    shutdown_server,
    write_server_info,
)
from zndraw.start_celery import run_celery_worker

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


def path_to_room(path: str) -> str:
    """Convert a file path to a valid room name by replacing non-alphanumeric characters."""
    room = re.sub(r"[^a-zA-Z0-9\-]", "_", path)
    return room


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
    port: int = 5000,
    debug: bool = False,
    verbose: bool = False,
    celery: bool = True,
    storage_path: str = typer.Option(
        "./zndraw-data",
        envvar="ZNDRAW_STORAGE_PATH",
        help="Path to storage directory for trajectory data (LMDB files per room)",
    ),
    redis_url: str | None = typer.Option(
        None,
        help="Redis server URL (e.g., `redis://localhost:6379`). If not provided, an in-memory storage will be used.",
    ),
    host: str = typer.Option(
        "localhost",
        envvar="ZNDRAW_SERVER_HOST",
        help="Server hostname or IP address for the SERVER_URL (e.g., 'example.com' or '192.168.1.1')",
    ),
    force_new_server: bool = typer.Option(
        False,
        "--force-new-server",
        help="Always start a new server instance, ignoring any existing server.",
    ),
    connect: str | None = typer.Option(
        None,
        "--connect",
        help="Connect to a specified remote server (e.g., https://zndraw.myserver.com), bypassing local server checks.",
    ),
    status: bool = typer.Option(
        False,
        "--status",
        help="Check if a local server is running and display its status.",
    ),
    shutdown: bool = typer.Option(
        False,
        "--shutdown",
        help="Stop the local server if it is running.",
    ),
    browser: bool = typer.Option(
        True,
        "--browser/--no-browser",
        help="Automatically open the web browser.",
    ),
    file_browser: bool = typer.Option(
        False,
        "--file-browser/--no-file-browser",
        help="Enable local filesystem browser endpoint.",
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
    simgen: bool = typer.Option(
        False,
        "--simgen/--no-simgen",
        envvar="ZNDRAW_SIMGEN_ENABLED",
        help="Enable SiMGen molecular generation features.",
    ),
):
    """
    Start or connect to a ZnDraw server.

    By default, this command will check if a local server is already running and connect to it.
    If no server is running, it will start a new one.

    Use --status to check server status or --shutdown to stop the server.
    """
    # Validate flag combinations
    if detached and (status or shutdown or connect):
        typer.echo(
            "Error: --detached cannot be used with --status, --shutdown, or --connect",
            err=True,
        )
        raise typer.Exit(1)
    if verbose:
        import logging

        logging.basicConfig(level=logging.DEBUG)

    # Handle --status flag
    if status:
        is_running, server_info, status_message = get_server_status()

        if is_running and server_info is not None:
            typer.echo(f"✓ {status_message}")
            typer.echo(f"  Server URL: http://localhost:{server_info.port}")
            raise typer.Exit(0)
        else:
            if server_info is None:
                typer.echo("✗ No local ZnDraw server is running")
            else:
                typer.echo(f"✗ {status_message}")
                typer.echo("  (You may want to run 'zndraw --shutdown' to clean up)")
            raise typer.Exit(1)
        return

    # Handle --shutdown flag
    if shutdown:
        server_info = read_server_info()

        if server_info is None:
            typer.echo("No server.pid file found. No server to shut down.")
            raise typer.Exit(0)

        if not is_process_running(server_info.pid):
            typer.echo(f"Process {server_info.pid} is not running.")
            typer.echo("Cleaning up server.pid file...")
            remove_server_info()
            raise typer.Exit(0)

        typer.echo(f"Shutting down server (PID: {server_info.pid})...")
        if shutdown_server(server_info):
            if not is_process_running(server_info.pid):
                typer.echo("✓ Server shut down successfully")
                remove_server_info()
                raise typer.Exit(0)
            else:
                typer.echo("⚠ Server process may still be running")
                raise typer.Exit(1)
        else:
            typer.echo("✗ Failed to shut down server")
            raise typer.Exit(1)
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
            first_room = path_to_room(path[0])

            # Open browser before uploading so user can watch progress
            if browser:
                browser_url = f"{connect}/rooms/{first_room}"
                typer.echo(f"Opening browser at {browser_url}")
                webbrowser.open(browser_url)

            for p in path:
                room = path_to_room(p)
                typer.echo(f"  Uploading file {p} to room {room}")
                read_file(
                    file=p,
                    room=room,
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

    # Check for existing server unless force_new_server is set
    if not force_new_server:
        is_running, server_info, status_message = get_server_status()

        if is_running and server_info is not None:
            typer.echo(f"✓ Found existing server: {status_message}")
            typer.echo(f"  Server URL: http://localhost:{server_info.port}")

            # Check version compatibility
            if server_info.version != __version__:
                typer.echo(
                    f"⚠ Warning: Server version ({server_info.version}) "
                    f"differs from CLI version ({__version__})"
                )
                typer.echo(
                    "  Consider running 'zndraw shutdown' and starting a new server."
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
                first_room = path_to_room(path[0])

                # Open browser before uploading so user can watch progress
                if browser:
                    browser_url = f"http://localhost:{server_info.port}/rooms/{first_room}"
                    typer.echo(f"Opening browser at {browser_url}")
                    webbrowser.open(browser_url)

                for p in path:
                    room = path_to_room(p)
                    typer.echo(f"  Uploading file {p} to room {room}")
                    read_file(
                        file=p,
                        room=room,
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
        elif server_info is not None:
            # Server info exists but server is not running - clean up
            typer.echo(f"Found stale server info: {status_message}")
            typer.echo("Cleaning up and starting a new server...")
            remove_server_info()

    # Validate files exist before starting server
    if path is not None:
        for p in path:
            if not os.path.exists(p):
                typer.echo(f"✗ Error: File not found: {p}", err=True)
                raise typer.Exit(1)

    # Start a new server
    if force_new_server:
        typer.echo("Starting a new server instance (--force-new-server)...")
    else:
        typer.echo("No existing server found. Starting a new server...")

    # if one file, promote to template and default
    # if multiple files, print a list of rooms.
    typer.echo(
        f"Rooms: {[path_to_room(p) for p in path]}"
        if path
        else "No files loaded on startup."
    )

    # Load configuration from environment and apply CLI overrides
    from zndraw.config import get_config
    config = get_config()

    # Override config with CLI arguments
    if redis_url is not None:
        config.redis_url = redis_url
    config.storage_path = storage_path
    config.server_host = host
    config.server_port = port
    config.simgen_enabled = simgen
    config.file_browser_enabled = file_browser
    if file_browser_root is not None:
        config.file_browser_root = file_browser_root
    else:
        config.file_browser_root = os.getcwd()
    config.celery_enabled = celery

    # Revalidate after overrides
    config._validate()

    # Check for stale LMDB data on startup when using in-memory mode
    # This can happen if the server was killed without clean shutdown
    if config.redis_url is None and os.path.exists(config.storage_path):
        typer.echo(f"⚠️  Warning: Found existing LMDB data at {config.storage_path}")
        typer.echo("   This data may be stale when using in-memory storage (no Redis).")
        if typer.confirm("   Delete and start fresh?", default=True):
            shutil.rmtree(config.storage_path)
            typer.echo(f"✓ Cleaned storage directory: {config.storage_path}")
        else:
            typer.echo("⚠️  Continuing with existing data (may cause state inconsistencies)")

    flask_app = create_app(config=config)

    # Track the first room for browser opening
    first_room = None
    if path is not None:
        first_room = path_to_room(path[0])

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

    # Write server info to PID file
    server_info = ServerInfo(pid=os.getpid(), port=port, version=__version__)
    write_server_info(server_info)
    typer.echo(f"✓ Server started (PID: {server_info.pid}, Port: {port})")
    typer.echo(f"  Server URL: {config.server_url}")

    # Queue file loading tasks after worker is started
    if path is not None:
        make_default = True
        for idx, p in enumerate(path):
            room = path_to_room(p)
            typer.echo(f"Loading file {p} into room {room}.")
            read_file.delay(
                file=p,
                room=room,
                server_url=config.server_url,
                start=start,
                stop=stop,
                step=step,
                make_default=make_default,
            )
            make_default = False

    # Open browser if requested
    if browser:
        if first_room:
            browser_url = (
                f"http://localhost:{port}/rooms/{first_room}"
            )
            typer.echo(f"Opening browser at {browser_url}")
        else:
            browser_url = f"http://localhost:{port}"
            typer.echo(f"Opening browser at {browser_url}")
        webbrowser.open(browser_url)

    try:
        socketio.run(flask_app, debug=debug, host="0.0.0.0", port=config.server_port)
    finally:
        # Clean up LMDB storage for in-memory mode
        # In-memory mode has ephemeral Redis state, so LMDB data becomes stale
        if config.redis_url is None and os.path.exists(config.storage_path):
            typer.echo(f"🧹 Cleaning LMDB storage: {config.storage_path}")
            shutil.rmtree(config.storage_path)

        flask_app.extensions["redis"].flushall()
        if celery:
            # Properly terminate celery worker and all its child processes
            worker.terminate()
            try:
                worker.wait(timeout=5)  # Wait up to 5 seconds for graceful shutdown
            except subprocess.TimeoutExpired:
                typer.echo(
                    "Celery worker did not shut down gracefully, forcing termination..."
                )
                worker.kill()  # Force kill if it doesn't shut down
                worker.wait()
            typer.echo("Celery worker closed.")
        # Clean up PID file when server stops
        remove_server_info()
        typer.echo("Server stopped.")
