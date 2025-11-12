# Remove
- "extension:register" socket handler

# Workflow

1. python-client registers a new extension

```
Client ->|REST| Server: POST /api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/submit
# request includes the sessionId to uniquly identify the worker

# If the extension is not yet registered OR the schame is in agreement with the exiting one
Server --> Client: 202 Accepted + job_id
# if the extension is already registered in that room or in the global registry and has a differet schema
Server --> Client: XXX Rejected, reason=...

Server puts worker into "idle" state, with timestamp
```

2. server emits the new schema to all connected clients in the room if per room or to all if public
```
Server --> |socket| --> client: new schema
Server --> |socket| --> client: incr available workers for given category/extension
```


3. client submits a job to the Server via REST

```
Client ->|REST| Server: POST /api/rooms/<string:room_id>/extensions/<string:category>/<string:extension>/submit
Server --> Client: 202 Accepted + job_id
```

4. server assigns the job to a worker
```
Server --> |socket| --> Worker: assign job id, for workers which are not in "processing"
(this uses socket call and waits for ack!)
And server pushes the task into "submited".
Server puts Worker into "assigned" state, with timestamp
```

5. worker picks up the job and starts processing
```
Worker --> |REST| server: PUT /api/<job:id>/status {status: processing}
Server --> Worker: 200 OK
Server puts Worker into "processing" state, with timestamp
```
6. worker updates progress
```
Worker --> |REST| server: PUT /api/<job:id>/status {status: finished} OR {status: failed}
Server --> Worker: 200 OK
Server puts Worker into "idle" state, with timestamp
```


Additionally, the user can cancel a job:
```
User →|REST| Server: DELETE /api/jobs/{jobId}
Server --> Client: 200 OK
If the task was assigned to a worker, the server notifies the worker via socket:
Server --> |socket| --> Worker: cancel job id
```

# heartbeat and disconnects
The server has a celery has an eventlet (celery in production) background cronjob that checks every 90 s if all available workers have sent a hearbeat in the last 90 s. If not, their state will change to "offline". If they were in the processing state / had a job assigned, that job will be set to failed.

Failed jobs are NOT retried!

On socket connect, we store the sessionId / socketID mapping.
On a disconnect of a client, we find the sessionId, unregister the worker and set its state to offline.
If the worker had a job assigned or was processing, that job is set to failed.

A worker runs a background thread that sends every 90 s a heartbeat to the server via REST:

```
Worker --> |REST| server: PUT /api/heartbeat (with sessionId in payload)
Server --> Worker: 200 OK
```

HEARTBEAT_INTERVAL = 90  # Worker sends every 90s
HEARTBEAT_CHECK_INTERVAL = 90  # Server checks every 90s  
HEARTBEAT_TIMEOUT = 2  # Miss 2 heartbeats = 60s = dead

# Task status
A task can have the following status:
- submited
- assigned
- processing
- finished
- failed

# Task assignment
Tasks are assigned in a round-robin fashion to all available workers that are idle.
Once a worker is set to "idle" e.g. by registering or finishing a task, the server checks if there are pending tasks in the queue that can be assigned to that worker and assigns them.

# ==================== WORKERS ====================
# Workers are independent entities that register with rooms

POST   /api/workers/register
PUT    /api/workers/{workerId}/heartbeat
POST   /api/workers/{workerId}/shutdown
GET    /api/workers/{workerId}
GET    /api/workers/{workerId}/jobs
DELETE /api/workers/{workerId}  # Force unregister

# ==================== ROOMS ====================
# Rooms own jobs and extensions

# Job submission (user action)
POST   /api/rooms/{roomId}/extensions/{category}/{extension}/submit

# Room-scoped queries
GET    /api/rooms/{roomId}/jobs
GET    /api/rooms/{roomId}/extensions
GET    /api/rooms/{roomId}/workers  # Workers registered in this room

# ==================== JOBS ====================
# Jobs are identified by globally unique ID

GET    /api/jobs/{jobId}
PUT    /api/jobs/{jobId}/status
DELETE /api/jobs/{jobId}  # Cancel job

# Global job queries (admin)
GET    /api/jobs?status=failed&roomId=room-123

# ==================== USER ====================
GET    /api/users/{userId}/jobs



1. Server assigns job to worker
   Server -->|socket| Worker: emit('job:assigned', {jobId: '123'})
   Server: job.status = 'assigned', worker.state = 'assigned'
   OR if celery
   run_celery_task.delay(jobId)

2. Worker receives notification and fetches details
   Worker -->|REST| Server: GET /api/jobs/123
   Server --> Worker: 200 OK + {jobId, roomId, ...}

3. Worker validates and starts processing
   Worker -->|REST| Server: PUT /api/jobs/123/status {status: 'processing'}
   Server --> Worker: 200 OK
   Server: job.status = 'processing', worker.state = 'processing'
