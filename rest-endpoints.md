# POST /api/disconnect/<string:client_sid>
Disconnects a client from the room.

# POST /internal/emit
Internal endpoint to emit Socket.IO events.
```json
{
  "event": "event_name",
  "sid": "client_sid",
  "data": {}
}
```

# POST /api/frames/<string:room_id>
Serves multiple frames' data from the room's Zarr store. Can use indices or a slice.
```json
{
  "indices": [0, 1, 2],
  "keys": ["positions", "energy"]
}
```
or
```json
{
  "start": 0,
  "stop": 10,
  "step": 2,
  "keys": ["positions", "energy"]
}
```

# POST /api/rooms/<string:room_id>/frames
Appends, replaces, extends, or inserts a new frame. Requires a bearer token. The data is msgpack encoded.

# GET /api/rooms
List all active rooms.

# GET /api/rooms/<string:room_id>
Get details for a specific room.

# GET /api/templates
List all available room templates.

# GET /api/templates/default
Get the default template.

# PUT /api/templates/default
Set the default template.
```json
{
  "template_id": "template_name"
}
```

# POST /api/rooms/<string:room_id>/promote
Promote a room to a template.
```json
{
  "name": "template_name",
  "description": "template_description"
}
```

# GET /api/exit
Gracefully shut down the server.

# GET /api/rooms/<string:room_id>/schema/<string:category>
Get the schema for a specific category (selections, modifiers, settings) in a room.

# POST /api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>
Logs a user extension action and creates a job.
```json
{
  "userId": "user_id",
  "data": {}
}
```

# GET /api/rooms/<string:room_id>/extension-data/<string:category>/<string:extension>?userId=<user_id>
Get the data for a specific extension for a user.

# GET /api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/workers
Get worker details for a specific extension.

# GET /api/rooms/<string:room_id>/jobs
List active jobs for a room.

# GET /api/rooms/<string:room_id>/jobs/<string:job_id>
Get details for a specific job.

# DELETE /api/rooms/<string:room_id>/jobs/<string:job_id>
Delete a job.

# POST /api/rooms/<string:room_id>/jobs/<string:job_id>/complete
Mark a job as completed.
```json
{
  "result": {},
  "workerId": "worker_id"
}
```

# POST /api/rooms/<string:room_id>/jobs/<string:job_id>/fail
Mark a job as failed.
```json
{
  "error": "error_message",
  "workerId": "worker_id"
}
```

# POST /api/rooms/<string:room_id>/extensions/register
Register a client-side extension.
```json
{
  "name": "extension_name",
  "category": "extension_category",
  "schema": {},
  "clientId": "client_id"
}
```

# GET /api/workers/<string:worker_id>
Get the current state of a worker.

# POST /api/rooms/<string:room_id>/jobs/next
Poll for the next available job for a worker.
```json
{
  "workerId": "worker_id"
}
```

# POST /api/room/<string:room_id>/join
Join a room, creating it if it doesn't exist.
```json
{
  "template": "template_id"
}
```
