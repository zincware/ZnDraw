import contextlib
import io
import traceback

import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import ase
from fastmcp import FastMCP
from zndraw.server_manager import get_server_status, ServerInfo
from zndraw import ZnDraw
import subprocess


# Initialize the MCP Server
mcp = FastMCP("ZnDraw Bridge")


# --- 2. Resources (Reading State) ---
@mcp.resource("docs://geometries")
def get_geometry_docs() -> str:
    """
    How to add / modify geometries in ZnDraw
    """
    from zndraw.geometries import geometries

    return f"""
ZnDraw Geometries are manages via `vis.geometries`!
Available geometries are: {", ".join(geometries.keys())}.
More information about each geometry can be found at `docs://geometry/<geometry_name>`.

You can set geometries using
vis.geometries['GeometryName'] = GeometryClass(...)
"""


@mcp.resource("docs://geometry/{geometry_name}")
def get_geometry_doc(geometry_name: str) -> str:
    """
    Get documentation for a specific geometry.
    """
    from zndraw.geometries import geometries

    geometry_cls = geometries.get(geometry_name)
    if geometry_cls is None:
        return f"Geometry '{geometry_name}' not found."

    return geometry_cls.model_json_schema().__str__()


@mcp.resource("zndraw://state")
def get_current_state() -> tuple[bool, ServerInfo | None, str]:
    """
    Returns the current status of the ZnDraw server.
    """
    return get_server_status()


# --- 3. Tools (Actions) ---


@mcp.tool()
def upload(path: str) -> tuple[bool, str]:
    """Upload a local file to the ZnDraw server."""
    is_running, _, status_message = get_server_status()
    if not is_running:
        return (
            False,
            f"Cannot upload: ZnDraw server is not running. Status: {status_message}",
        )

    result = subprocess.run(["zndraw", path], capture_output=True)
    if result.returncode != 0:
        return False, f"Upload failed: {result.stderr.decode()}"
    return True, f"Upload successful: {result.stdout.decode()}"


@mcp.tool()
def run_zndraw_script(
    python_code: str,
    url: str | None = None,
    room: str = "default",
    user: str | None = None,
    password: str | None = None,
) -> tuple[bool, str]:
    """
    Executes a custom Python script to analyze or manipulate the ZnDraw scene.

    CONTEXT:
    - 'vis' is PRE-LOADED (The ZnDraw connection).
    - 'np' (numpy) is PRE-LOADED.
    - 'ase' is PRE-LOADED.
    - 'go' (plotly.graph_objects) is PRE-LOADED.
    - 'px' (plotly.express) is PRE-LOADED.
    - 'pd' (pandas) is PRE-LOADED.

    INSTRUCTIONS:
    - Use 'print()' to output results to the chat.
    - To modify the scene, call methods on 'vis'.
    - Do not wrap code in markdown blocks (```), just send raw code.
    - vis is `MutableSequence[ase.Atoms]` with extra methods for ZnDraw.
    - update vis.step = len(vis) - 1 to show the progress.

    Further Documentation:
    - use docs://geometries for geometry docs.
    """
    is_running, _, status_message = get_server_status()
    if not is_running:
        return (
            False,
            f"Cannot upload: ZnDraw server is not running. Status: {status_message}",
        )

    vis = ZnDraw(url=url, room=room, user=user, password=password)

    # 1. Prepare the execution sandbox
    sandbox = {
        "vis": vis,
        "np": np,
        "ase": ase,
        "go": go,
        "px": px,
        "pd": pd,
    }

    # 2. Capture Stdout
    output_buffer = io.StringIO()

    try:
        with contextlib.redirect_stdout(output_buffer):
            # 3. Execute the code
            exec(python_code, {}, sandbox)

        result = output_buffer.getvalue()
        return True, f"Output:\n{result}"

    except Exception:
        # 4. Return full traceback so the Agent can debug its own code
        return False, f"Runtime Error:\n{traceback.format_exc()}"


# --- 4. Main Entry Point ---
if __name__ == "__main__":
    print("Starting ZnDraw MCP Server...")
    mcp.run()
