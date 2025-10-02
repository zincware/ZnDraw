# Socket Events

## connect
Handles a new client connection.

## disconnect
Handles a client disconnection, cleaning up locks and room presence.

## join_room
Joins a client to a specific room.
```json
{
  "room": "room_id",
  "userId": "user_id",
  "clientId": "client_id"
}
```

## selection:set
Sets the selection of atoms.
```json
{
  "indices": [0, 1, 2]
}
```

## frame_selection:set
Sets the selection of frames.
```json
{
  "indices": [0, 1, 2]
}
```

## bookmarks:set
Sets the bookmarks for frames.
```json
{
  "bookmarks": {
    "0": "first frame",
    "10": "tenth frame"
  }
}
```

## set_frame_atomic
Sets the current frame for all clients. Can only be used when no presenter is active.
```json
{
  "frame": 10
}
```

## set_frame_continuous
Continuously updates the frame. Requires the sender to be the presenter.
```json
{
  "frame": 10
}
```

## request_presenter_token
Requests to become the presenter.

## release_presenter_token
Releases the presenter lock.

## lock:acquire
Acquires a lock on a target.
```json
{
  "target": "trajectory:meta"
}
```

## lock:release
Releases a lock on a target.
```json
{
  "target": "trajectory:meta"
}
```

## upload:prepare
Prepares for a file upload, returning a short-lived token.
```json
{
  "action": "append"
}
```

## frame:delete
Deletes a frame.
```json
{
  "frame_id": 10
}
```

---

# Emitted Events

## room_users_update
Sent when a user joins or leaves a room.
```json
{
  "sid1": "user1",
  "sid2": "user2"
}
```

## presenter_update
Sent when the presenter changes.
```json
{
  "presenterSid": "sid"
}
```

## selection:update
Sent when the atom selection is updated.
```json
{
  "indices": [0, 1, 2]
}
```

## frame_selection:update
Sent when the frame selection is updated.
```json
{
  "indices": [0, 1, 2]
}
```

## bookmarks:update
Sent when the bookmarks are updated.
```json
{
  "bookmarks": {
    "0": "first frame"
  }
}
```

## frame_update
Sent when the current frame is updated.
```json
{
  "frame": 10
}
```

## len_frames
Sent when the number of frames changes.
```json
{
  "success": true,
  "count": 20
}
```

## INVALIDATE_SCHEMA
Sent when a client-side extension is registered or a worker disconnects.
```json
{
  "roomId": "room_id",
  "category": "modifiers"
}
```
