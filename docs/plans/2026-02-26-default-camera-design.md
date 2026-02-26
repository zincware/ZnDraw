# Default Camera Feature Design

## Problem

New browser sessions always start with the hardcoded `Camera()` defaults.
Users want to set a room-wide default camera so that new sessions inherit
a specific view (position, target, FOV, etc.).

## Design

### Storage

Add a nullable `default_camera: str | None` column to the `Room` SQLModel.
The value is a geometry key (e.g., `"template-cam"`) referencing a `Camera`
in `RoomGeometry`. This follows the "persistent room state lives in SQL"
convention.

### REST API

Two new endpoints on the existing geometries router:

- `GET /v1/rooms/{room_id}/default-camera` → `{"default_camera": "key" | null}`
- `PUT /v1/rooms/{room_id}/default-camera` → body `{"default_camera": "key" | null}`
  - Validates geometry exists in `RoomGeometry` and is type `Camera`
  - Returns 404 if geometry not found, 400 if not a Camera
  - `null` unsets the default

### Server-side camera clone on join

In `room_join()` (socketio.py), when creating a session camera for frontend
clients: if `room.default_camera` is set, load the geometry from SQL and use
its properties (`position`, `target`, `up`, `fov`, `near`, `far`, `zoom`,
`camera_type`) as initial values instead of `Camera()` defaults. No extra
round-trip from the frontend.

### Geometry deletion cleanup

In `delete_geometry()`: if the deleted key matches `room.default_camera`,
set the column to `None` and commit.

### Python client

`vis.default_camera` property on `ZnDraw`:
- Getter: `GET /v1/rooms/{room_id}/default-camera`
- Setter: validates key exists in `vis.geometries` and is a `Camera`,
  then `PUT`. Setting `None` unsets.

### Frontend

- `useDefaultCamera()` hook (React Query GET + mutation PUT)
- Star icon column in `GeometryGrid.tsx` for Camera-type rows
  (excluding session cameras `cam:*`). Click toggles default on/off.

## Files

| File | Change |
|------|--------|
| `src/zndraw/models.py` | Add `default_camera` column to `Room` |
| `src/zndraw/routes/geometries.py` | GET/PUT endpoints + cleanup on delete |
| `src/zndraw/socketio.py` | Clone default camera in `room_join` |
| `src/zndraw/socket_events.py` | Add `default_camera` to `RoomJoinResponse` |
| `src/zndraw/client.py` | `default_camera` property |
| `frontend/src/hooks/useDefaultCamera.ts` | New hook |
| `frontend/src/components/geometry/GeometryGrid.tsx` | Star icon column |
| `frontend/src/myapi/client.ts` | `getDefaultCamera`, `setDefaultCamera` |
| `tests/test_default_camera.py` | Comprehensive tests |

## Not included

- **#858 (job timeout)**: Belongs in `distributed-taskiq`, not zndraw.
- **#863 (arrays.colors dropdown)**: Already handled by `DynamicEnumRenderer`
  fetching live from `useAvailableProperties`. CLI race handled by
  `copy_from="@none"`.
