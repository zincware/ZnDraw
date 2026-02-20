"""ZnDraw Agent CLI — structured JSON interface for LLM agents."""

from __future__ import annotations

from typing import Annotated

import typer

from .auth import auth_app
from .bookmarks import bookmarks_app
from .chat import chat_app
from .extensions import extensions_app
from .figures import figures_app
from .frames import frames_app
from .geometries import geometries_app
from .jobs import jobs_app
from .mount import mount_cmd
from .rooms import rooms_app
from .screenshots import screenshots_app
from .selection import selection_app
from .selection_groups import selection_groups_app
from .sessions import sessions_app
from .step import step_app

app = typer.Typer(
    name="zndraw-cli",
    help="ZnDraw Agent CLI — structured JSON interface for LLM agents",
    no_args_is_help=True,
    context_settings={"allow_extra_args": True, "allow_interspersed_args": False},
)

# Sub-apps (resource groups)
app.add_typer(auth_app, name="auth")
app.add_typer(rooms_app, name="rooms")
app.add_typer(frames_app, name="frames")
app.add_typer(step_app, name="step")
app.add_typer(selection_app, name="selection")
app.add_typer(selection_groups_app, name="selection-groups")
app.add_typer(extensions_app, name="extensions")
app.add_typer(jobs_app, name="jobs")
app.add_typer(chat_app, name="chat")
app.add_typer(geometries_app, name="geometries")
app.add_typer(bookmarks_app, name="bookmarks")
app.add_typer(figures_app, name="figures")
app.add_typer(screenshots_app, name="screenshots")
app.add_typer(sessions_app, name="sessions")

# Standalone commands
app.command("mount")(mount_cmd)


@app.callback()
def callback(
    ctx: typer.Context,
    url: Annotated[
        str | None,
        typer.Option(envvar="ZNDRAW_URL", help="ZnDraw server URL"),
    ] = None,
    token: Annotated[
        str | None,
        typer.Option(envvar="ZNDRAW_TOKEN", help="Auth token"),
    ] = None,
) -> None:
    """Global options for server connection."""
    ctx.ensure_object(dict)
    ctx.obj["url"] = url
    ctx.obj["token"] = token
