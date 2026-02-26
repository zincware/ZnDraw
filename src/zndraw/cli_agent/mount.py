from __future__ import annotations

import sys
from typing import Annotated

import typer

from zndraw.cli_agent.connection import resolve_url
from zndraw.cli_agent.output import json_print


def mount_cmd(
    ctx: typer.Context,
    file: Annotated[str, typer.Argument(help="Path to trajectory file")],
    room: Annotated[
        str | None,
        typer.Option(help="Room ID to mount into (default: new room). Must be empty."),
    ] = None,
) -> None:
    """Mount a trajectory file into a room (lazy frame serving).

    Opens the file lazily via asebytes and serves frames on demand
    through the provider system. Blocks until interrupted (Ctrl+C).
    """
    import asebytes

    from zndraw.client import ZnDraw

    url = resolve_url(ctx.obj["url"])
    db = asebytes.ASEIO(file)

    vis = ZnDraw(url=url, room=room)
    vis.mount(db)

    json_print(
        {
            "room_id": vis.room,
            "url": f"{url}/room/{vis.room}",
            "frame_count": len(db),
        }
    )
    sys.stdout.flush()

    try:
        vis.wait()
    except KeyboardInterrupt:
        pass
    finally:
        vis.disconnect()
