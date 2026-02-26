from __future__ import annotations

from io import StringIO
from typing import Annotated

import typer

from zndraw.cli_agent.connection import get_connection, get_current_step
from zndraw.cli_agent.output import json_print, text_print
from zndraw.schemas import FrameBulkResponse, StatusResponse, StepResponse

frames_app = typer.Typer(name="frames", help="Frame operations")


@frames_app.command("count")
def count(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
) -> None:
    """Get total number of frames in the room."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    resp = conn.get(f"/v1/rooms/{room}/step")
    json_print(StepResponse.model_validate(resp.json()))


@frames_app.command("get")
def get(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    keys: Annotated[
        str | None, typer.Option(help="Comma-separated keys to include")
    ] = None,
    format: Annotated[str, typer.Option(help="Output format: json or xyz")] = "json",
) -> None:
    """Get a single frame by index."""
    import ase.io
    import asebytes
    import msgpack

    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    if index is None:
        index = get_current_step(conn, room)
    params = {}
    if keys:
        params["keys"] = keys

    resp = conn.get(f"/v1/rooms/{room}/frames/{index}", params=params)

    frames_raw = msgpack.unpackb(resp.content, raw=True)
    frame_raw = frames_raw[0] if isinstance(frames_raw, list) else frames_raw
    atoms = asebytes.decode(frame_raw)

    if format == "xyz":
        buf = StringIO()
        ase.io.write(buf, atoms, format="extxyz")
        text_print(buf.getvalue())
    else:
        json_print(
            {
                "symbols": atoms.get_chemical_symbols(),
                "positions": atoms.positions.tolist(),
                "cell": atoms.cell.array.tolist(),
                "pbc": atoms.pbc.tolist(),
                "info": atoms.info,
            }
        )


@frames_app.command("list")
def list_frames(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    start: Annotated[int, typer.Option(help="Start index")] = 0,
    stop: Annotated[int | None, typer.Option(help="Stop index (exclusive)")] = None,
    keys: Annotated[
        str | None, typer.Option(help="Comma-separated keys to include")
    ] = None,
    format: Annotated[str, typer.Option(help="Output format: json or xyz")] = "json",
) -> None:
    """List frames in the room."""
    import ase.io
    import asebytes
    import msgpack

    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    params: dict = {"start": start}
    if stop is not None:
        params["stop"] = stop
    if keys:
        params["keys"] = keys

    resp = conn.get(f"/v1/rooms/{room}/frames", params=params)

    frames_raw = msgpack.unpackb(resp.content, raw=True)
    atoms_list = [asebytes.decode(f) for f in frames_raw]

    if format == "xyz":
        buf = StringIO()
        ase.io.write(buf, atoms_list, format="extxyz")
        text_print(buf.getvalue())
    else:
        result = [
            {
                "symbols": atoms.get_chemical_symbols(),
                "positions": atoms.positions.tolist(),
                "cell": atoms.cell.array.tolist(),
                "pbc": atoms.pbc.tolist(),
                "info": atoms.info,
            }
            for atoms in atoms_list
        ]
        json_print(result)


@frames_app.command("extend")
def extend(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    file: Annotated[str, typer.Option(help="Path to trajectory file")],
) -> None:
    """Extend the room trajectory with frames from a file."""
    import pathlib

    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    path = pathlib.Path(file)
    with path.open("rb") as f:
        resp = conn.post(
            f"/v1/rooms/{room}/trajectory",
            files={"file": (path.name, f, "application/octet-stream")},
        )
    json_print(FrameBulkResponse.model_validate(resp.json()))


@frames_app.command("export")
def export(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    format: Annotated[str, typer.Option(help="Output format (default: xyz)")] = "xyz",
    indices: Annotated[
        str | None, typer.Option(help="Index range, e.g. '0:10'")
    ] = None,
    selection: Annotated[bool, typer.Option(help="Export only selected atoms")] = False,
) -> None:
    """Export trajectory from room to stdout."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])

    if format == "json":
        import asebytes
        import msgpack

        params: dict = {}
        if indices:
            params["indices"] = indices
        if selection:
            params["selection"] = "true"

        resp = conn.get(f"/v1/rooms/{room}/frames", params=params)
        frames_raw = msgpack.unpackb(resp.content, raw=True)
        result = [
            {
                "symbols": atoms.get_chemical_symbols(),
                "positions": atoms.positions.tolist(),
                "cell": atoms.cell.array.tolist(),
                "pbc": atoms.pbc.tolist(),
                "info": atoms.info,
            }
            for atoms in (asebytes.decode(f) for f in frames_raw)
        ]
        json_print(result)
    else:
        params = {"format": format}
        if indices:
            params["indices"] = indices
        if selection:
            params["selection"] = "true"

        resp = conn.get(f"/v1/rooms/{room}/trajectory", params=params)
        text_print(resp.text)


@frames_app.command("delete")
def delete(
    ctx: typer.Context,
    room: Annotated[str, typer.Argument(help="Room ID")],
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
) -> None:
    """Delete a frame by index."""
    conn = get_connection(ctx.obj["url"], ctx.obj["token"])
    if index is None:
        index = get_current_step(conn, room)
    response = conn.delete(f"/v1/rooms/{room}/frames/{index}")
    json_print(StatusResponse.model_validate(response.json()))
