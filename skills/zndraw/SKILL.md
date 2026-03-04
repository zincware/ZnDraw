---
name: zndraw
description: Use when interacting with a ZnDraw molecular visualization server — reading frames, navigating trajectories, selecting atoms, bookmarking, creating plots, running extensions, or performing computational chemistry analysis on atomic structures
---

# ZnDraw Agent Skill

## Overview

ZnDraw is an interactive visualization platform for atomistic simulations. Two interaction tools:

1. **`uv run zndraw-cli`** — Typer CLI, structured stdout. For CRUD on rooms, frames, selections, bookmarks, figures, extensions, ... .
2. **`uv run python -c "..."`** — For computation requiring numpy/ASE/plotly. **Must use the `ZnDraw` Python client** to connect to the server and modify data. Use `uv run --with <dependencies>` to include extra packages.

## Discover First, Never Guess

The CLI is self-documenting. **Always discover before acting:**

```bash
uv run zndraw-cli --help                                          # all resource groups
uv run zndraw-cli frames --help                                   # verbs + options for frames
uv run zndraw-cli extensions list --room ROOM                     # available extensions (fully qualified names)
uv run zndraw-cli extensions describe --room ROOM NAME            # parameter schema
uv run zndraw-cli frames get --room ROOM 0                        # inspect one frame to learn the data shape
```

Extension names are fully qualified (e.g. `@internal:modifiers:Delete`). Never hardcode names — they vary per server. Run `extensions list` then `extensions describe` to learn parameters before `extensions run`.

## Connection

Both `uv run zndraw-cli` and the Python client require a **running ZnDraw server**. Start one with `uv run zndraw` if none is running. Verify with `uv run zndraw-cli rooms list` (exit code 0 = server reachable).

```bash
# Auto-discovers running server + creates guest token
uv run zndraw-cli rooms list

# Explicit (--url, --token, --room are per-subcommand options)
uv run zndraw-cli rooms list --url http://localhost:8000 --token TOKEN

# Env vars: ZNDRAW_URL, ZNDRAW_TOKEN, ZNDRAW_ROOM
# When working on a single room, export ZNDRAW_ROOM to avoid repeating it
export ZNDRAW_ROOM=my-room-id
uv run zndraw-cli frames count
```

## Authentication

The CLI auto-discovers tokens in this order:
1. `--token` flag or `ZNDRAW_TOKEN` env var
2. Stored token in `~/.zndraw/tokens.json` (validated on use, removed on 401)
3. Guest token (anonymous fallback)

To use your browser identity from the CLI (required for screenshots and session control):

```bash
# Regular user: browser-based login (opens browser for approval)
uv run zndraw-cli auth login

# Login to a specific server (--url is a per-subcommand option)
uv run zndraw-cli auth login --url http://myserver:8000

# Check current identity and admin status
uv run zndraw-cli auth status

# Remove stored token
uv run zndraw-cli auth logout
```

For admin users who need to act as another user:

```bash
# Check if you're admin
uv run zndraw-cli auth status    # look for "is_superuser": true

# List all users
uv run zndraw-cli admin users list

# Login as another user (mints a token)
uv run zndraw-cli admin users login USER_ID
```

After `auth login` or `admin users login`, the token is stored in `~/.zndraw/tokens.json` and automatically used for all subsequent CLI commands — no `--token` needed.

## CLI Quick Reference

Use `<command> --help` for full options. Key patterns:

| Resource | Verbs |
|----------|-------|
| `auth` | `login`, `status`, `logout` |
| `admin` | `users list`, `users login` |
| `rooms` | `list`, `create`, `info`, `open`, `lock`, `unlock`, `set-default` |
| `frames` | `count`, `get`, `list`, `extend FILE`, `export`, `delete` |
| `step` | `get`, `set` |
| `selection` | `get`, `set`, `clear` |
| `selection-groups` | `list`, `get`, `set`, `delete` |
| `bookmarks` | `list`, `set`, `delete` |
| `figures` | `list`, `get`, `set` (`--file` or `--data`), `delete` |
| `chat` | `list`, `send` |
| `extensions` | `list`, `describe`, `run` |
| `jobs` | `list`, `status` |
| `geometries` | `list`, `get`, `set`, `delete`, `types`, `describe` |
| `preset` | `list`, `get`, `load`, `apply`, `save`, `export`, `delete` |
| `screenshots` | `list`, `request`, `get` |
| `sessions` | `list`, `get-camera`, `set-camera` |
| Standalone | `mount FILE` |

**Frames tips:**
- `--format xyz` for extxyz output, default is JSON
- `--keys arrays.positions,info.energy` for server-side field filtering (works on both `get` and `list`). Keys use dotted prefixes: `arrays.*` for per-atom data (positions, numbers, colors, radii), `info.*` for per-frame metadata (energy, connectivity). Structural fields (`cell`, `pbc`) are always included.
- `--start N --stop M` for ranges (don't fetch 10k frames at once)
- `export --indices "0:10"` supports Python range syntax (`start:stop` or `start:stop:step`)
- `export --selection` exports only currently selected atoms
- `extend` goes through ASE file I/O which may reclassify `info` fields (e.g. `energy` → `atoms.calc`). Use Python client `vis.append()` for full fidelity.

**Chat tips:**
- `send` has full markdown support
- Use single quotes for messages with special characters: `'Hello!'`
- Send relevant updates to the chat in addition to printing them.

**Geometry tips:**
- `set KEY --data '{"resolution": 32}'` without `--type` performs a partial update (PATCH, deep merge)
- `set KEY --type Sphere --data '{...}'` with `--type` performs a full replace (PUT)

## Extensions: Discover → Describe → Run → Poll

```bash
# 1. Discover what's available
uv run zndraw-cli extensions list --room ROOM

# 2. Get parameter schema for the extension you need
uv run zndraw-cli extensions describe --room ROOM "@internal:modifiers:AddFromSMILES"

# 3. Run with --key value flags
#    Returns {"job_id": "abc123", "status": "pending"}
uv run zndraw-cli extensions run --room ROOM "@internal:modifiers:AddFromSMILES" --smiles "CCO"

# 4. Poll for completion (job_id is global, no room needed)
uv run zndraw-cli jobs status JOB_ID
```

Three extension categories exist: **modifiers** (edit atoms), **selections** (pick atoms), **analysis** (create plots). Discover which are available via `extensions list`.

**Note:** Custom extensions registered via `vis.register_job()` are scoped to the room they were registered in. Built-in extensions (`@internal:...`) are available in all rooms.

**Known limitation:** Custom extensions registered via `register_job()` are currently not discoverable via `extensions list` / `extensions describe`. Only `@internal:*` extensions appear. Use string dispatch (`vis.run("@internal:...")`) for built-in extensions.

## Visual Presets

Named visual styles (materials, lighting, fog) that can be applied to room geometries. Bundled presets (`matt`, `flat`, `glossy`, `pathtracing`) are always available in every room without needing to be created. Creating a room-level preset with the same name overrides the bundled one; deleting the override makes the bundled version reappear.

```bash
# List presets in a room
uv run zndraw-cli preset list --room ROOM

# Apply a bundled or custom preset
uv run zndraw-cli preset apply --room ROOM matt

# Save current geometry state as a new preset
uv run zndraw-cli preset save --room ROOM my-style -p "particles*" -p "fog" -p "*light*"

# Export/import presets as JSON for sharing
uv run zndraw-cli preset export --room ROOM my-style ./my-style.json
uv run zndraw-cli preset load --room ROOM ./my-style.json
```

Presets use fnmatch patterns (`particles*`, `*light*`, `fog`) and optional `geometry_type` filters to target specific geometries. Rules are deep-merged — only specified config keys are overridden.

## Python Scripting with ZnDraw Client

When the CLI alone isn't enough (filtering, math, custom analysis, Plotly figures), use `uv run python -c "..."` with the **`ZnDraw` Python client**. It connects to the server and provides a `MutableSequence[ase.Atoms]` interface.

```python
from zndraw import ZnDraw
vis = ZnDraw(room="ROOM")
# After `zndraw-cli auth login`, stored token is auto-discovered — no --token needed
# Or with explicit url and token: ZnDraw(url=..., room=..., token="JWT_TOKEN")

# Frames — MutableSequence[ase.Atoms]
len(vis)              # frame count
vis[0]                # frame as ase.Atoms
vis[vis.step]         # current frame
vis.step = 42         # jump to frame
vis.append(atoms)     # append frame
del vis[0]            # delete frame
vis[:]                # all frames as list[ase.Atoms]

# Selections, bookmarks, figures
vis.selection = [0, 1, 2]       # set selected atom indices
vis.bookmarks[42] = "label"     # MutableMapping[int, str]
vis.figures["energy"] = fig     # MutableMapping[str, plotly.Figure]

# Chat
vis.chat.send("message")        # send chat message
vis.chat[0]                      # read message (Sequence)

# Room state
vis.locked = True                # lock/unlock room

# Run extensions — string name + kwargs (recommended)
task_id = vis.run("@internal:modifiers:Delete")
task_id = vis.run("@internal:modifiers:AddFromSMILES", smiles="CCO")

# Discovery
list(vis.extensions)                                    # all extension names
vis.extensions["@internal:modifiers:Delete"]["schema"]   # parameter schema

# Poll task
vis.tasks[task_id]               # TaskResponse
vis.tasks("running")             # filtered view

# Sessions — room-scoped Mapping of active browser sessions
vis.sessions                     # Mapping[str, Session] (all users in room)
list(vis.sessions)               # list of session SIDs

# Screenshots (requires own browser session in the room)
# Screenshots can only be captured for sessions belonging to the authenticated
# user. The Python client creates a new guest user by default — use `token=`
# or the stored CLI token (after `auth login`) to share identity with your
# browser session.
sids = list(vis.sessions)           # list active session SIDs
session = vis.sessions[sids[0]]     # get a session by SID
img = session.screenshot()          # capture screenshot (own sessions only)
img.data                            # PNG bytes
img.save("frame.png")               # save to file

# Visual Presets — MutableMapping[str, Preset]
list(vis.presets)                 # list preset names
vis.presets.apply("matt")        # apply preset to room geometries
vis.presets.load(Path("f.json")) # load preset from JSON file
vis.presets.export("pub", Path("out.json"))  # export to file
vis.presets["custom"] = Preset(name="custom", rules=[...])
del vis.presets["custom"]

# Efficient partial reads (server-side filtering)
vis.get(slice(None), keys=["info.energy"])

# Classmethods (no room needed)
ZnDraw.list_rooms() # or with url="http://localhost:8000"
token = ZnDraw.login(username="...", password="...") # or with url="http://localhost:8000",
```

### Bookmark frames matching a condition

```bash
uv run python -c "
from zndraw import ZnDraw
vis = ZnDraw(room='ROOM')
for i, atoms in enumerate(vis):
    if atoms.info.get('energy', 0) > -40:
        vis.bookmarks[i] = f'E={atoms.info[\"energy\"]:.2f}'
"
```

### Plot energy over time

```bash
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
vis = ZnDraw(room='ROOM')
data = vis.get(slice(None), keys=['info.energy'])
energies = [f['info.energy'] for f in data]
fig = go.Figure(go.Scatter(x=list(range(len(energies))), y=energies, mode='lines+markers'))
fig.update_layout(title='Energy', xaxis_title='Frame', yaxis_title='Energy (eV)')
vis.figures['energy'] = fig
"
```

### Select atoms by element

```bash
uv run python -c "
from zndraw import ZnDraw
vis = ZnDraw(room='ROOM')
atoms = vis[vis.step]
vis.selection = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'H']
"
```

### Compute RDF with ASE

```bash
uv run python -c "
from zndraw import ZnDraw
import plotly.graph_objects as go
from ase.geometry.analysis import Analysis
vis = ZnDraw(room='ROOM')
ana = Analysis(vis[:])
rdf = ana.get_rdf(rmax=6.0, nbins=100, return_dists=True)
fig = go.Figure(go.Scatter(x=rdf[1][0].tolist(), y=rdf[0][0].tolist(), mode='lines'))
fig.update_layout(title='RDF', xaxis_title='r (A)', yaxis_title='g(r)')
vis.figures['rdf'] = fig
"
```

## Mounting Large Files
For large files that support sliced access (`.h5`, `.lmdb`, `.zarr`,`.db`), use `mount` — frames are served lazily on demand.

```bash
# Mount into a new room (auto-generated ID)
uv run zndraw-cli mount trajectory.h5

# Mount into a new specific room
uv run zndraw-cli mount trajectory.h5 --room my-room-id

# Output: {"room_id": "...", "url": "http://...", "frame_count": 50000}
```
This command keeps serving frames until Ctrl+C is pressed.


## Camera and Session Control

Each browser session has an active camera. Use session commands to manage cameras:

```bash
# List sessions to find SIDs
uv run zndraw-cli sessions list --room ROOM

# Get a session's current camera
uv run zndraw-cli sessions get-camera --room ROOM SESSION_ID

# Switch a session's camera (geometry key must be a Camera type)
uv run zndraw-cli sessions set-camera --room ROOM SESSION_ID my-camera
```

Camera geometries can be partially updated like any geometry:
```bash
# Adjust just the FOV without resetting position
uv run zndraw-cli geometries set --room ROOM my-camera --data '{"fov": 45}'
```

## Task Mapping

| User request | Approach |
|-------------|----------|
| "jump to frame 42" | `uv run zndraw-cli step set --room ROOM 42` |
| "how many frames?" | `uv run zndraw-cli frames count --room ROOM` |
| "select all hydrogens" | Python: `vis.selection = [i for i, s in ...]` |
| "bookmark high-energy frames" | Python: iterate `vis`, set `vis.bookmarks[i]` |
| "plot energy over time" | Python: `vis.get(...)` + Plotly → `vis.figures[...]` |
| "compute RDF / analysis" | Python: `vis[:]` + ASE Analysis → `vis.figures[...]` |
| "delete selected atoms" | `extensions describe` the Delete modifier → `extensions run` |
| "add a molecule" | `extensions describe` AddFromSMILES → `extensions run` with params |
| "show me aspirin / molecule" | `extensions list` → find AddFromSMILES → `extensions describe` → `extensions run --smiles "..."` |
| "visualize this large H5 file" | `uv run zndraw-cli mount trajectory.h5` (lazy) |
| "open / show the room" | `uv run zndraw-cli rooms open --room ROOM` |
| "lock the room" | `uv run zndraw-cli rooms lock --room ROOM` or Python: `vis.locked = True` |
| "send a chat message" | `uv run zndraw-cli chat send --room ROOM "msg"` or Python: `vis.chat.send("msg")` |
| "take a screenshot" | `auth login` first, then `screenshots request --room ROOM` (CLI)|
| "list browser sessions" | `uv run zndraw-cli sessions list --room ROOM` (shows all users' sessions) |
| "who am I / admin check" | `uv run zndraw-cli auth status` |
| "apply matt style" | `uv run zndraw-cli preset apply --room ROOM matt` or Python: `vis.presets.apply("matt")` |
| "save current look as preset" | `uv run zndraw-cli preset save --room ROOM my-style -p "particles*" -p "fog"` |
| "share a preset" | `uv run zndraw-cli preset export --room ROOM my-style ./my-style.json` |
| "what can I do?" | `extensions list --room ROOM` + `--help` on each resource group |

## Common Mistakes

- **Running bare `zndraw-cli`** — always use `uv run zndraw-cli`; the binary is not on PATH
- **Guessing extension names/params** — always `list` then `describe` first
- **Fetching everything** — use `vis.get(..., keys=[...])` or CLI `--keys` to limit data
- **Not polling jobs** — `extensions run` returns a `job_id`; check `jobs status JOB_ID` for completion
- **Wrong room** — `rooms list` to discover, don't guess UUIDs
- **Piping CLI to Python for mutations** — use the `ZnDraw` client instead; it connects directly and handles serialization
