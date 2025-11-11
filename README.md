## Quick Start

### Local Development
1. Run the server: `uv run src/server.py`

### Docker Deployment (Recommended)

**Production:**
```bash
cd docker
docker compose up -d
```
Access at `http://localhost` (port 80)

**Development:**
```bash
cd docker/zndraw
docker compose up -d
```
Access at `http://localhost:5000`

For detailed Docker setup instructions, see the [Docker Deployment](#docker-deployment) section below.

# Goal
Split communication into data channels and control channels.
- Data channels: use HTTP PUT/POST to upload/download data. Use query strings to select what should be loaded, e.g. `?data=positions,species` or just`?data=energy`.
- Control channels: use Socket.IO and Redis for state management and locks

Optional: Control logic, the server knows the previous and the next frame and can check what needs to be updated, e.g. only positions or also species, box, ...

Data persistence should be abstracted away via a `DataProvider` interface.
- store e.g. as npy files for fast access once accessed once. (Possibility to hash??)

Data Edits
- not only store the frame id but also the version `{ "frame": 10, "version": 123 }` would also allow undo operations, (also conflict detection -> considering the edit was made on version 1 but the current version is 3, the edit cannot be applied)
- use / broadcast a lock mechanism (redis lock?) e.g. if the data is to be modified we need a lock, data can only be updated if the lock is held by the client that wants to update it (server-side check). Per-frame lock or even more granular? Timeouts: what if a client crashes mid-edit, socketio disconnect event should trigger? Renewals: if an edit takes long, client must refresh lock before expiry. Server authority: server must check that only the lock-holder can commit.
- with `vis.lock` ... ? and `vis.extend` will check if the lock has been aquired, otherwise aquire it.
- To add / remove and to edit, there exist two entries, one for trajectory indices and another one for metadata per frame. To insert at index N, you must first move the last frame (M) to M+1, then M-1 to M, and so on, until you move frame N to N+1.
```
# 1. Trajectory-Level: A single key holding the list of all frames.
"trajectory:indices"  (Sorted Set)
  -> [ 0, 1, 3, 4, ... ]

# 2. Frame-Level: A separate key for each frame's metadata.
"metadata:frame:0"    (Hash)
  -> { "version": 5, "hash": "..." }

"metadata:frame:1"    (Hash)
  -> { "version": 2, "hash": "..." }

"metadata:frame:3"    (Hash)
  -> { "version": 8, "hash": "..." }
```

Rooms
- use `from flask_socketio import join_room, leave_room` and `
```javascript
socket.on('connect', () => {
  socket.emit('join_room', { room: room_id });
});
```
and use keys like `room:project-alpha:lock:frame:10` and `/data/<room>/...`

For playback, we use the `"Presenter Token" Hybrid Model` approach.
1. Acquire Token: When a user starts scrubbing, their client first asks the server for the "presenter token."
socket.emit('request_presenter_token')

2. Server Grants Token: The server uses SETNX in Redis to grant the token to the first client who asks. This token has a short expiry (e.g., 5 seconds).
r.set("room:project-alpha:presenter_lock", client_id, nx=True, ex=5)

3. Stream Updates: If the client gets the token, it then starts sending the stream of set_room_frame messages just like in your model. The server only accepts these messages from the client who currently holds the token.

4. Renew Token: As long as the user is actively scrubbing, the client sends a renew_presenter_token message every few seconds to extend the expiry.

5. Release/Expire: When the user stops scrubbing, their client explicitly sends a release_presenter_token message. If they close the tab or crash, the token automatically expires in Redis after 5 seconds, allowing someone else to take control.

Conflict detection
- Server only accepts update if client’s base version == current version, If mismatch → reject with "conflict", client must rebase on latest version.

Data hashing
- Detect corruption / incomplete transfer.
- Enable deduplication (same frame content doesn’t need to be re-saved).
- store alongside in db `{ "frame": 10, "version": 123, "hash": "sha256:abcd..." }`

Split into room and client-ids?

Use https://msgpack.org/

How to handle `default` room
Store data in `zarr` format. Use padding for variable length data.

When creating a new room, one should be able to specify the base data.


## Summary
✅ Things that are already great
	•	Split channels
	•	Data over HTTP (upload/download, efficient for MB payloads).
	•	Control over Socket.IO + Redis (events, synchronization, locks).
	•	This avoids mixing heavy + light traffic.
	•	DataProvider abstraction
	•	Hides persistence details (HDF5, .npy, object store, …).
	•	Makes it possible to swap implementations without touching business logic.
	•	You can even wrap in-memory cache behind the interface.
	•	Versioning ({frame, version})
	•	Enables undo/redo and conflict detection (if a client edits an outdated version).
	•	Very good idea for collaborative edits.
	•	Locks (Redis or server-side)
	•	Prevents race conditions (two clients trying to edit same frame).
	•	Redis locks scale well if you have multiple server processes.

⸻

🔹 Improvements & Missing Pieces

1. Lock granularity
	•	You probably want per-frame locks, not global locks.
	•	i.e. lock:frame:10 instead of lock:all.
	•	This way two clients can edit different frames simultaneously.
	•	Use Redis SETNX with expiry for robust distributed locks:

r.set(f"lock:frame:{frame_id}", client_id, nx=True, ex=10)



2. Lock lifecycle
	•	You’ll need a strategy for:
	•	Timeouts: what if a client crashes mid-edit?
	•	Renewals: if an edit takes long, client must refresh lock before expiry.
	•	Server authority: server must check that only the lock-holder can commit.

👉 The Redis community has patterns here (Redlock algorithm).

⸻

3. Conflict detection
	•	Versioning is great, but you should enforce:
	•	Server only accepts update if client’s base version == current version.
	•	If not, reject or trigger conflict resolution.

client -> { "frame": 10, "base_version": 123, "new_data": ... }
server -> compares base_version with stored version (123 vs 124)

	•	If mismatch → reject with "conflict", client must rebase on latest version.

⸻

4. Data hashing
	•	Adding hashes has two benefits:
	•	Detect corruption / incomplete transfer.
	•	Enable deduplication (same frame content doesn’t need to be re-saved).
	•	You can store hash in Redis along with version:

{ "frame": 10, "version": 123, "hash": "sha256:abcd..." }



⸻

5. Persistence backend strategy
	•	.npy is great for fast local random access.
	•	For large-scale / multi-machine setup:
	•	Object storage (S3, MinIO) with pre-signed URLs for uploads/downloads.
	•	Zarr/TileDB for chunked array storage.
	•	Your DataProvider abstraction should make it easy to swap between these.

⸻

6. Partial updates (optional optimization)
	•	If edits are often small (e.g. just moving atoms), sending the entire frame blob is overkill.
	•	You could allow patches (delta updates).
	•	Control channel:

{ "frame": 10, "base_version": 123, "delta": { "positions": [...] } }


	•	Server applies delta → increments version → broadcasts update.

Not critical to start with, but worth keeping in mind.

⸻

7. Broadcast strategy
	•	When a frame is updated:
	•	Don’t broadcast the heavy payload.
	•	Only broadcast metadata { "frame": 10, "version": 124 }.
	•	Clients then fetch updated frame data via HTTP if they care about it.

This keeps the control channel lightweight.

⸻

🔹 Refined workflow (with your ideas included)
	1.	Client requests lock
	•	socket.io.emit("lock_request", {frame: 10})
	•	Server acquires Redis lock → grants if free.
	2.	Client uploads data
	•	HTTP POST /upload/frame/10?lock_token=xyz with new blob.
	3.	Server checks & commits
	•	Verifies lock holder.
	•	Verifies base_version matches current version.
	•	Writes new data (DataProvider.save).
	•	Increments version, stores hash.
	•	Releases lock.
	4.	Server broadcasts update
	•	socket.io.emit("frame_updated", {frame: 10, version: 124}).
	5.	Other clients fetch new data
	•	On receiving "frame_updated", request updated frame from data channel.

⸻

✅ This gives you:
	•	Fast uploads/downloads (HTTP).
	•	Strong synchronization (Redis locks + versions).
	•	Robust persistence (swappable backends).
	•	Undo/redo potential (versions + hashes).
	•	No race conditions (lock + version checks).


# Docker Deployment

ZnDraw provides two Docker setups:

1. **Production** (`/docker/`) - Nginx + load balancing + multiple replicas
2. **Development** (`/docker/zndraw/`) - Single instance for local testing

## Prerequisites

- Docker (version 20.10 or later)
- Docker Compose (version 2.0 or later)

## Production Deployment

**Location:** `/docker/docker-compose.yaml`

Full production setup with:
- Nginx reverse proxy + static file serving
- 3 × ZnDraw app replicas (load balanced)
- 2 × Celery workers (parallel task processing)
- Redis with persistence
- Sticky sessions for WebSocket support

```bash
cd docker

# Build and start all services
docker compose up -d

# View logs
docker compose logs -f

# Stop services
docker compose down
```

Access: `http://localhost` (port 80)

**📖 See `/docker/README.md` for full production documentation**

## Development/Testing

**Location:** `/docker/zndraw/docker-compose.yaml`

Simple single-instance setup for development:

```bash
cd docker/zndraw

# Build and start services
docker compose up -d

# View logs
docker compose logs -f

# Stop services
docker compose down
```

Access: `http://localhost:5000`

**📖 See `/docker/zndraw/README.md` for development documentation**

## Configuration

### Environment Variables

Edit the appropriate `docker-compose.yaml` file (`/docker/` or `/docker/zndraw/`) to configure environment variables.

#### Core Server Configuration

| Variable | Description | Default | CLI Argument | Used By |
|----------|-------------|---------|--------------|---------|
| `ZNDRAW_REDIS_URL` | Redis connection URL for coordination | None (in-memory mode) | `--redis-url` | Main app, Celery workers |
| `ZNDRAW_STORAGE_PATH` | Base path for LMDB storage files | `./zndraw-data` | `--storage-path` | Main app, Celery workers |
| `ZNDRAW_SERVER_HOST` | Server bind host | `localhost` | `--host` | Main app |
| `ZNDRAW_SERVER_URL` | Full server URL for callbacks | `http://localhost:5000` | Derived from host/port | Celery workers |
| `ZNDRAW_LOG_LEVEL` | Logging level | `WARNING` | N/A | Main app, Celery workers |

#### Security

| Variable | Description | Default | CLI Argument | Used By |
|----------|-------------|---------|--------------|---------|
| `FLASK_SECRET_KEY` | Flask session secret key | `dev-secret-key-change-in-production` | N/A | Main app |
| `ZNDRAW_ADMIN_USERNAME` | Admin username (optional) | Not set | N/A | Main app |
| `ZNDRAW_ADMIN_PASSWORD` | Admin password (optional) | Not set | N/A | Main app |

**⚠️ Important:** Change `FLASK_SECRET_KEY` to a random secret in production!

#### Upload & Storage Limits

| Variable | Description | Default | CLI Argument | Used By |
|----------|-------------|---------|--------------|---------|
| `ZNDRAW_UPLOAD_TEMP` | Temporary upload directory | `/tmp/zndraw_uploads` | N/A | Main app, Celery workers |
| `ZNDRAW_MAX_UPLOAD_MB` | Maximum upload size in MB | `500` | N/A | Main app |

#### Optional Features

| Variable | Description | Default | CLI Argument | Used By |
|----------|-------------|---------|--------------|---------|
| `ZNDRAW_SIMGEN_ENABLED` | Enable SiMGen molecular generation | `false` | `--simgen/--no-simgen` | Main app |
| `ZNDRAW_FILE_BROWSER_ENABLED` | Enable file browser feature | `false` | `--file-browser` | Main app |
| `ZNDRAW_FILE_BROWSER_ROOT` | File browser root directory | Current directory | `--file-browser-root` | Main app |

#### System (Internal)

| Variable | Description | Default | Used By |
|----------|-------------|---------|---------|
| `OBJC_DISABLE_INITIALIZE_FORK_SAFETY` | Fix for Celery on Apple Silicon | `YES` | Celery workers (macOS only) |


- have the state of the queue available in the chat or a table or something