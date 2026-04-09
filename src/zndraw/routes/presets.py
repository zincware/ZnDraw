"""Presets REST API endpoints for room visual presets."""

import fnmatch
import json
from datetime import UTC, datetime
from functools import lru_cache

from fastapi import APIRouter
from sqlmodel import select

from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    SioDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    InvalidPresetRule,
    NotAuthenticated,
    PresetAlreadyExists,
    PresetNotFound,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.geometries import geometries as geometry_models
from zndraw.models import RoomGeometry, RoomPreset
from zndraw.presets import list_bundled_presets
from zndraw.schemas import (
    Preset,
    PresetApplyResult,
    PresetRule,
    PresetsListResponse,
    StatusResponse,
    deep_merge,
)
from zndraw.socket_events import GeometryInvalidate

router = APIRouter(prefix="/v1/rooms/{room_id}/presets", tags=["presets"])


@lru_cache(maxsize=1)
def _bundled_presets() -> dict[str, Preset]:
    """Load bundled presets from package data (cached, immutable)."""
    result: dict[str, Preset] = {}
    for path in list_bundled_presets():
        data = json.loads(path.read_text())
        preset = Preset.model_validate(data)
        result[preset.name] = preset
    return result


def _validate_rules(rules: list[PresetRule]) -> None:
    """Validate preset rules against known geometry models."""
    for rule in rules:
        if rule.geometry_type is None:
            continue
        model_cls = geometry_models.get(rule.geometry_type)
        if model_cls is None:
            raise InvalidPresetRule.exception(
                f"Unknown geometry type: '{rule.geometry_type}'"
            )
        valid_fields = set(model_cls.model_fields.keys())
        invalid_keys = set(rule.config.keys()) - valid_fields
        if invalid_keys:
            raise InvalidPresetRule.exception(
                f"Invalid config keys for {rule.geometry_type}: {invalid_keys}"
            )


def _row_to_preset(row: RoomPreset) -> Preset:
    """Convert a RoomPreset DB row to a Preset response."""
    rules = [PresetRule.model_validate(r) for r in json.loads(row.rules)]
    return Preset(
        name=row.name,
        description=row.description,
        rules=rules,
        created_at=row.created_at,
        updated_at=row.updated_at,
    )


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_presets(
    session: SessionDep,
    _current_user: CurrentUserDep,
    room_id: str,
) -> PresetsListResponse:
    """List all presets for a room.

    Merges bundled presets with room-level DB presets.
    DB presets override bundled ones with the same name.
    """
    await verify_room(session, room_id)
    result = await session.exec(select(RoomPreset).where(RoomPreset.room_id == room_id))
    rows = result.all()

    # Start with bundled, let DB rows override
    merged: dict[str, Preset] = dict(_bundled_presets())
    for row in rows:
        merged[row.name] = _row_to_preset(row)

    return PresetsListResponse(items=list(merged.values()))


@router.get(
    "/{name}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, PresetNotFound),
)
async def get_preset(
    session: SessionDep,
    _current_user: CurrentUserDep,
    room_id: str,
    name: str,
) -> Preset:
    """Get a single preset by name.

    DB preset takes priority; falls back to bundled preset.
    """
    await verify_room(session, room_id)
    row = await session.get(RoomPreset, (room_id, name))
    if row is not None:
        return _row_to_preset(row)
    bundled = _bundled_presets().get(name)
    if bundled is not None:
        return bundled
    raise PresetNotFound.exception(f"Preset '{name}' not found")


@router.post(
    "",
    status_code=201,
    responses=problem_responses(
        NotAuthenticated,
        RoomNotFound,
        RoomLocked,
        PresetAlreadyExists,
        InvalidPresetRule,
    ),
)
async def create_preset(
    session: SessionDep,
    _room: WritableRoomDep,
    room_id: str,
    request: Preset,
) -> Preset:
    """Create a new preset."""
    _validate_rules(request.rules)

    existing = await session.get(RoomPreset, (room_id, request.name))
    if existing is not None:
        raise PresetAlreadyExists.exception(f"Preset '{request.name}' already exists")

    now = datetime.now(UTC)
    row = RoomPreset(
        room_id=room_id,
        name=request.name,
        description=request.description,
        rules=json.dumps([r.model_dump() for r in request.rules]),
        created_at=now,
        updated_at=now,
    )
    session.add(row)
    await session.commit()
    await session.refresh(row)

    return _row_to_preset(row)


@router.put(
    "/{name}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, InvalidPresetRule
    ),
)
async def upsert_preset(
    session: SessionDep,
    _room: WritableRoomDep,
    room_id: str,
    name: str,
    request: Preset,
) -> Preset:
    """Create or update a preset (idempotent)."""
    _validate_rules(request.rules)

    row = await session.get(RoomPreset, (room_id, name))
    now = datetime.now(UTC)

    if row is None:
        row = RoomPreset(
            room_id=room_id,
            name=name,
            description=request.description,
            rules=json.dumps([r.model_dump() for r in request.rules]),
            created_at=now,
            updated_at=now,
        )
        session.add(row)
    else:
        row.description = request.description
        row.rules = json.dumps([r.model_dump() for r in request.rules])
        row.updated_at = now

    await session.commit()
    await session.refresh(row)

    return _row_to_preset(row)


@router.delete(
    "/{name}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, PresetNotFound
    ),
)
async def delete_preset(
    session: SessionDep,
    _room: WritableRoomDep,
    room_id: str,
    name: str,
) -> StatusResponse:
    """Delete a preset."""
    row = await session.get(RoomPreset, (room_id, name))
    if row is None:
        raise PresetNotFound.exception(f"Preset '{name}' not found")
    await session.delete(row)
    await session.commit()
    return StatusResponse()


@router.post(
    "/{name}/apply",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, PresetNotFound
    ),
)
async def apply_preset(
    session: SessionDep,
    sio: SioDep,
    _room: WritableRoomDep,
    room_id: str,
    name: str,
) -> PresetApplyResult:
    """Apply a preset to all matching geometries in the room.

    Resolves from DB first, then bundled presets.
    The special ``@default`` name resets all geometries to factory defaults.
    """
    if name == "@default":
        from zndraw.routes.rooms import _initialize_default_geometries

        stmt = select(RoomGeometry).where(RoomGeometry.room_id == room_id)
        existing = (await session.exec(stmt)).all()
        existing_keys = [g.key for g in existing]
        for g in existing:
            await session.delete(g)
        await session.flush()
        _initialize_default_geometries(session, room_id)
        await session.commit()

        new_stmt = select(RoomGeometry).where(RoomGeometry.room_id == room_id)
        new_keys = [g.key for g in (await session.exec(new_stmt)).all()]
        all_keys = sorted(set(existing_keys) | set(new_keys))
        for key in all_keys:
            await sio.emit(
                GeometryInvalidate(room_id=room_id, operation="set", key=key),
                room=room_channel(room_id),
            )
        return PresetApplyResult(geometries_updated=all_keys)

    row = await session.get(RoomPreset, (room_id, name))
    if row is not None:
        rules = [PresetRule.model_validate(r) for r in json.loads(row.rules)]
    else:
        bundled = _bundled_presets().get(name)
        if bundled is None:
            raise PresetNotFound.exception(f"Preset '{name}' not found")
        rules = bundled.rules

    stmt = select(RoomGeometry).where(RoomGeometry.room_id == room_id)
    geometries = (await session.exec(stmt)).all()

    updated_keys: list[str] = []
    for geom in geometries:
        matching_configs: list[dict] = []
        for rule in rules:
            if not fnmatch.fnmatch(geom.key, rule.pattern):
                continue
            if rule.geometry_type is not None and rule.geometry_type != geom.type:
                continue
            matching_configs.append(rule.config)

        if not matching_configs:
            continue

        current_config = json.loads(geom.config)
        for config in matching_configs:
            current_config = deep_merge(current_config, config)

        # Validate merged config against geometry model
        model_cls = geometry_models.get(geom.type)
        if model_cls is None:
            continue
        validated = model_cls(**current_config)
        geom.config = validated.model_dump_json()
        updated_keys.append(geom.key)

    await session.commit()

    for key in updated_keys:
        await sio.emit(
            GeometryInvalidate(room_id=room_id, operation="set", key=key),
            room=room_channel(room_id),
        )

    return PresetApplyResult(geometries_updated=updated_keys)
