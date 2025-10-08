import os
import re
import shutil
import time
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
    celery: bool = True,
    storage_path: str = "./zndraw-data.zarr",
    redis_url: str = "redis://localhost:6379",
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
        "--browser",
        help="Automatically open the web browser.",
    ),
):
    """
    Start or connect to a ZnDraw server.

    By default, this command will check if a local server is already running and connect to it.
    If no server is running, it will start a new one.

    Use --status to check server status or --shutdown to stop the server.
    """
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
            typer.echo("Uploading files to remote server...")
            make_default = True
            for idx, p in enumerate(path):
                room = path_to_room(p)
                if idx == 0 and browser:
                    browser_url = f"{connect}/rooms/{room}?waitForCreation=true"
                    typer.echo(f"Opening browser at {browser_url}")
                    webbrowser.open(browser_url)
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
                typer.echo("Uploading files to existing server...")
                make_default = True
                server_url = f"http://localhost:{server_info.port}"
                for idx, p in enumerate(path):
                    room = path_to_room(p)
                    # Open browser after uploading files, pointing to first room
                    if browser and idx == 0:
                        browser_url = f"http://localhost:{server_info.port}/rooms/{room}?waitForCreation=true"
                        typer.echo(f"Opening browser at {browser_url}")
                        webbrowser.open(browser_url)
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

    flask_app = create_app(storage_path=storage_path, redis_url=redis_url)
    server_url = f"http://{host}:{port}"
    flask_app.config["SERVER_URL"] = server_url

    # Track the first room for browser opening
    first_room = None
    if path is not None:
        make_default = True
        for idx, p in enumerate(path):
            room = path_to_room(p)
            if idx == 0:
                first_room = room
            typer.echo(f"Loading file {p} into room {room}.")
            read_file.delay(
                file=p,
                room=room,
                server_url=server_url,
                start=start,
                stop=stop,
                step=step,
                make_default=make_default,
            )
            make_default = False

    if celery:
        worker = run_celery_worker()

    # Write server info to PID file
    server_info = ServerInfo(pid=os.getpid(), port=port, version=__version__)
    write_server_info(server_info)
    typer.echo(f"✓ Server started (PID: {server_info.pid}, Port: {port})")
    typer.echo(f"  Server URL: http://localhost:{port}")

    # Open browser if requested
    if browser:
        # If we have a file being loaded, open to that room with waitForCreation flag
        if first_room:
            browser_url = f"http://localhost:{port}/rooms/{first_room}?waitForCreation=true"
            typer.echo(f"Opening browser at {browser_url}")
        else:
            browser_url = f"http://localhost:{port}"
            typer.echo(f"Opening browser at {browser_url}")
        webbrowser.open(browser_url)

    try:
        socketio.run(flask_app, debug=debug, host="0.0.0.0", port=port)
    finally:
        shutil.rmtree("data", ignore_errors=True)
        flask_app.extensions["redis"].flushall()
        if celery:
            worker.terminate()
            worker.wait()
            typer.echo("Celery worker closed.")
        # Clean up PID file when server stops
        remove_server_info()
        typer.echo("Server stopped.")
