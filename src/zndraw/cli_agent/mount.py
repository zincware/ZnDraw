from __future__ import annotations

import sys
from typing import Annotated

import typer

from zndraw.cli_agent.connection import (
    PasswordOpt,
    RoomOpt,
    TokenOpt,
    UrlOpt,
    UserOpt,
    get_connection,
)
from zndraw.cli_agent.output import json_print


def mount_cmd(
    file: Annotated[str, typer.Argument(help="Path to trajectory file")],
    url: UrlOpt = None,
    token: TokenOpt = None,
    user: UserOpt = None,
    password: PasswordOpt = None,
    room: RoomOpt = None,
) -> None:
    """Mount a trajectory file into a room (lazy frame serving).

    Opens the file lazily via open_frames and serves frames on demand
    through the provider system. Blocks until interrupted (Ctrl+C).
    """
    from collections.abc import Iterator

    from zndraw.client import ZnDraw
    from zndraw.io import open_frames

    conn = get_connection(url, token, user, password)
    source = open_frames(file)
    # mount() requires random-access (len + __getitem__);
    # materialise streaming iterators into a list.
    db = list(source) if isinstance(source, Iterator) else source

    vis = ZnDraw(url=conn.base_url, room=room, token=conn.token)
    vis.mount(db)

    json_print(
        {
            "room_id": vis.room,
            "url": f"{conn.base_url}/room/{vis.room}",
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
        conn.close()
