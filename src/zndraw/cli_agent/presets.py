"""CLI commands for visual preset operations."""

from __future__ import annotations

import fnmatch
import pathlib
from typing import Annotated

import typer

from zndraw.geometries.base import BaseGeometry
from zndraw.schemas import Preset, PresetRule

from .connection import (
    RoomOpt,
    TokenOpt,
    UrlOpt,
    cli_error_handler,
    get_zndraw,
    resolve_room,
)
from .output import json_print

presets_app = typer.Typer(name="preset", help="Visual preset operations")


@presets_app.command("list")
def list_presets(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """List all presets in a room (name + description summary)."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        json_print(
            [
                {"name": name, "description": vis.presets[name].description}
                for name in vis.presets
            ]
        )
        vis.disconnect()


@presets_app.command("get")
def get_preset(
    name: Annotated[str | None, typer.Argument(help="Preset name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Get a preset by name."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Preset name is required")
        vis = get_zndraw(url, token, room)
        json_print(vis.presets[name])
        vis.disconnect()


@presets_app.command("load")
def load_preset(
    path: Annotated[
        pathlib.Path | None, typer.Argument(help="Path to preset JSON file")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Load a preset from a JSON file into the room."""
    with cli_error_handler():
        room = resolve_room(room)
        if path is None:
            raise typer.BadParameter("Path to preset JSON file is required")
        vis = get_zndraw(url, token, room)
        preset = vis.presets.load(path)
        json_print(preset)
        vis.disconnect()


@presets_app.command("apply")
def apply_preset(
    name: Annotated[str | None, typer.Argument(help="Preset name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Apply a preset to the room's geometries."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Preset name is required")
        vis = get_zndraw(url, token, room)
        result = vis.presets.apply(name)
        json_print(result)
        vis.disconnect()


@presets_app.command("save")
def save_preset(
    name: Annotated[str | None, typer.Argument(help="Preset name to create")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
    patterns: Annotated[
        list[str] | None,
        typer.Option("--pattern", "-p", help="Geometry key patterns to include"),
    ] = None,
) -> None:
    """Save current geometry state as a new preset.

    Uses the user-provided patterns as rule patterns (for portability across rooms).
    Each pattern that matches at least one geometry becomes a rule containing
    the config from the first matched geometry of each type.

    Example:
        zndraw-cli preset save matt --room my-room -p "particles*" -p "fog" -p "*light*"
    """
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Preset name is required")
        vis = get_zndraw(url, token, room)

        patterns = patterns or ["*"]
        rules: list[PresetRule] = []
        for pattern in patterns:
            type_groups: dict[str, list[BaseGeometry]] = {}
            for key in vis.geometries:
                if fnmatch.fnmatch(key, pattern):
                    geom = vis.geometries[key]
                    gtype = type(geom).__name__
                    type_groups.setdefault(gtype, []).append(geom)

            for gtype, geoms in type_groups.items():
                rules.append(
                    PresetRule(
                        pattern=pattern,
                        geometry_type=gtype,
                        config=geoms[0].model_dump(),
                    )
                )

        preset = Preset(
            name=name,
            description=f"Saved from room {room}",
            rules=rules,
        )
        vis.presets[name] = preset
        json_print(preset)
        vis.disconnect()


@presets_app.command("reset")
def reset_preset(
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Reset all geometries to factory defaults."""
    with cli_error_handler():
        room = resolve_room(room)
        vis = get_zndraw(url, token, room)
        result = vis.presets.apply("@default")
        json_print(result)
        vis.disconnect()


@presets_app.command("export")
def export_preset(
    name: Annotated[str | None, typer.Argument(help="Preset name")] = None,
    output: Annotated[
        pathlib.Path | None, typer.Argument(help="Output JSON path")
    ] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Export a preset from the room to a JSON file for sharing."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Preset name is required")
        if output is None:
            raise typer.BadParameter("Output path is required")
        vis = get_zndraw(url, token, room)
        vis.presets.export(name, output)
        typer.echo(f"Exported preset '{name}' to {output}")
        vis.disconnect()


@presets_app.command("delete")
def delete_preset(
    name: Annotated[str | None, typer.Argument(help="Preset name")] = None,
    url: UrlOpt = None,
    token: TokenOpt = None,
    room: RoomOpt = None,
) -> None:
    """Delete a preset from the room."""
    with cli_error_handler():
        room = resolve_room(room)
        if name is None:
            raise typer.BadParameter("Preset name is required")
        vis = get_zndraw(url, token, room)
        del vis.presets[name]
        typer.echo(f"Deleted preset '{name}'")
        vis.disconnect()
