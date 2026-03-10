from __future__ import annotations

import sys
from typing import Annotated

import typer

from zndraw.cli_agent.connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    resolve_token,
    resolve_url,
)
from zndraw.cli_agent.output import json_print


def mount_cmd(
    file: Annotated[str, typer.Argument(help="Path to trajectory file")],
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Mount a trajectory file into a room (lazy frame serving).

    Opens the file lazily via asebytes and serves frames on demand
    through the provider system. Blocks until interrupted (Ctrl+C).
    """
    import asebytes

    from zndraw.client import ZnDraw

    resolved_url = resolve_url(url)
    resolved_token = resolve_token(resolved_url, token)
    db = asebytes.ASEIO(file)

    vis = ZnDraw(url=resolved_url, room=room, token=resolved_token)
    vis.mount(db)

    json_print(
        {
            "room_id": vis.room,
            "url": f"{resolved_url}/room/{vis.room}",
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
