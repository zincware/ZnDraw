from __future__ import annotations

from io import StringIO
from typing import Annotated

import typer

from zndraw.cli_agent.connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from zndraw.cli_agent.output import json_print, text_print
from zndraw.schemas import FrameBulkResponse, StatusResponse, StepResponse

frames_app = typer.Typer(name="frames", help="Frame operations")


@frames_app.command("count")
def count(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get total number of frames in the room."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        data = vis.api.get_step()
        json_print(StepResponse.model_validate(data))
        vis.disconnect()


@frames_app.command("get")
def get(
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    keys: Annotated[
        str | None, typer.Option(help="Comma-separated keys to include")
    ] = None,
    format: Annotated[str, typer.Option(help="Output format: json or xyz")] = "json",
) -> None:
    """Get a single frame by index."""
    import ase.io
    import asebytes
    import msgpack
    from asebytes import atoms_to_dict

    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        if index is None:
            index = vis.step
        params = {}
        if keys:
            params["keys"] = keys

        # Binary endpoint — use raw HTTP
        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/frames/{index}",
            params=params,
            headers=vis.api.get_headers(),
        )
        vis.api.raise_for_status(resp)

        frames_raw = msgpack.unpackb(resp.content, raw=True)
        frame_raw = frames_raw[0] if isinstance(frames_raw, list) else frames_raw
        atoms = asebytes.decode(frame_raw)

        if format == "xyz":
            buf = StringIO()
            ase.io.write(buf, atoms, format="extxyz")
            text_print(buf.getvalue())
        else:
            json_print(atoms_to_dict(atoms))
        vis.disconnect()


@frames_app.command("list")
def list_frames(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
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
    from asebytes import atoms_to_dict

    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        params: dict = {"start": start}
        if stop is not None:
            params["stop"] = stop
        if keys:
            params["keys"] = keys

        resp = vis.api.http.get(
            f"/v1/rooms/{vis.room}/frames",
            params=params,
            headers=vis.api.get_headers(),
        )
        vis.api.raise_for_status(resp)

        frames_raw = msgpack.unpackb(resp.content, raw=True)
        atoms_list = [asebytes.decode(f) for f in frames_raw]

        if format == "xyz":
            buf = StringIO()
            ase.io.write(buf, atoms_list, format="extxyz")
            text_print(buf.getvalue())
        else:
            json_print([atoms_to_dict(atoms) for atoms in atoms_list])
        vis.disconnect()


@frames_app.command("extend")
def extend(
    file: Annotated[str | None, typer.Argument(help="Path to trajectory file")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Extend the room trajectory with frames from a file."""
    import pathlib

    with cli_error_handler():
        room = resolve_room(room)
        if file is None:
            raise typer.BadParameter("FILE is required")
        vis = get_zndraw(url, token, room)
        path = pathlib.Path(file)
        with path.open("rb") as f:
            resp = vis.api.http.post(
                f"/v1/rooms/{vis.room}/trajectory",
                files={"file": (path.name, f, "application/octet-stream")},
                headers=vis.api.get_headers(),
            )
        vis.api.raise_for_status(resp)
        json_print(FrameBulkResponse.model_validate(resp.json()))
        vis.disconnect()


def _expand_indices(raw: str) -> str:
    """Expand colon range syntax to comma-separated indices.

    ``"0:10"`` → ``"0,1,...,9"``; ``"0:10:2"`` → ``"0,2,4,6,8"``.
    Plain comma-separated values pass through unchanged.
    """
    if ":" in raw:
        parts = raw.split(":")
        if not all(p.lstrip("-").isdigit() for p in parts if p):
            raise typer.BadParameter(f"Invalid index range: {raw!r}")
        return ",".join(str(i) for i in range(*[int(p) for p in parts]))
    return raw


@frames_app.command("export")
def export(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    format: Annotated[str, typer.Option(help="Output format (default: xyz)")] = "xyz",
    indices: Annotated[
        str | None, typer.Option(help="Index range, e.g. '0:10'")
    ] = None,
    selection: Annotated[bool, typer.Option(help="Export only selected atoms")] = False,
) -> None:
    """Export trajectory from room to stdout."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)

        # Expand colon range to comma-separated indices
        expanded_indices = _expand_indices(indices) if indices else None

        # Resolve --selection to actual atom indices
        selection_param: str | None = None
        if selection:
            sel = vis.selections.get("particles", ())
            if sel:
                selection_param = ",".join(str(i) for i in sel)

        if format == "json":
            import asebytes
            import msgpack
            from asebytes import atoms_to_dict

            params: dict = {}
            if expanded_indices:
                params["indices"] = expanded_indices
            if selection_param:
                params["selection"] = selection_param

            resp = vis.api.http.get(
                f"/v1/rooms/{vis.room}/frames",
                params=params,
                headers=vis.api.get_headers(),
            )
            vis.api.raise_for_status(resp)
            frames_raw = msgpack.unpackb(resp.content, raw=True)
            json_print([atoms_to_dict(asebytes.decode(f)) for f in frames_raw])
        else:
            params = {"format": format}
            if expanded_indices:
                params["indices"] = expanded_indices
            if selection_param:
                params["selection"] = selection_param

            resp = vis.api.http.get(
                f"/v1/rooms/{vis.room}/trajectory",
                params=params,
                headers=vis.api.get_headers(),
            )
            vis.api.raise_for_status(resp)
            text_print(resp.text)
        vis.disconnect()


@frames_app.command("delete")
def delete(
    index: Annotated[
        int | None, typer.Argument(help="Frame index (default: current step)")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a frame by index."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        if index is None:
            index = vis.step
        vis.api.delete_frame(index)
        json_print(StatusResponse(status="ok"))
        vis.disconnect()
