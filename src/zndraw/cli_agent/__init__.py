"""ZnDraw Agent CLI — structured JSON interface for LLM agents."""

from __future__ import annotations

import typer

from .admin import admin_app
from .auth import auth_app
from .bookmarks import bookmarks_app
from .chat import chat_app
from .extensions import extensions_app
from .figures import figures_app
from .frames import frames_app
from .geometries import geometries_app
from .jobs import jobs_app
from .mount import mount_cmd
from .presets import presets_app
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
)

# Sub-apps (resource groups)
app.add_typer(admin_app, name="admin")
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
app.add_typer(presets_app, name="preset")

# Standalone commands
app.command("mount")(mount_cmd)
