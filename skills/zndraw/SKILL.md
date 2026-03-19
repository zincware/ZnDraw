---
name: zndraw
description: Use when interacting with a ZnDraw molecular visualization server — reading frames, navigating trajectories, selecting atoms, bookmarking, creating plots, running extensions, or performing computational chemistry analysis on atomic structures
---

# ZnDraw Agent Skill

## Overview

ZnDraw is an interactive visualization platform for atomistic simulations. Two interaction tools:

1. **`uv run zndraw-cli`** — Typer CLI, structured stdout. For CRUD on rooms, frames, selections, bookmarks, figures, extensions, ... .
2. **`uv run python -c "..."`** — **Last resort** for computation requiring numpy/ASE/plotly that no extension can handle. **Must use the `ZnDraw` Python client** to connect to the server and modify data. Use `uv run --with <dependencies>` to include extra packages.

## Extensions First — MANDATORY

**For ANY action that modifies atoms, creates structures, selects atoms, or runs analysis, you MUST check extensions before writing Python code or using the CLI directly.**

The server has extensions (modifiers, selections, analysis) that handle most tasks. Extensions you don't expect exist — e.g. PackBox for building simulation boxes, various selection methods, analysis tools. **You cannot know what's available without checking.**

**Required workflow for every action task:**

```
1. uv run zndraw-cli extensions list          # ALWAYS do this first
2. Scan the list — does an extension match the task?
   YES → extensions describe NAME → extensions run NAME
   NO  → Use CLI commands or Python client as fallback
```

**Do NOT skip step 1.** Even if you think you know how to do it in Python, an extension may do it better, faster, and with proper server-side integration.

**Red flags — STOP and check extensions first:**
- You're about to write `uv run python -c "..."` for anything other than reading data or making plots
- You're about to use ASE to build/modify structures programmatically
- You think "I'll just write a quick script to..."
- The task involves creating, deleting, duplicating, or transforming atoms

The CLI is self-documenting. Use `--help` on any command to discover options:

```bash
uv run zndraw-cli --help                                          # all resource groups
uv run zndraw-cli frames --help                                   # verbs + options for frames
uv run zndraw-cli frames get 0                                    # inspect one frame to learn the data shape
```

Extension names are fully qualified (e.g. `@internal:modifiers:Delete`). Always quote them in shell commands (`@` and `:` may trigger shell expansion). Never hardcode names — they vary per server.

## Connection & Session Setup

Both `uv run zndraw-cli` and the Python client require a **running ZnDraw server**. Start one with `uv run zndraw` if none is running. Verify with `uv run zndraw-cli rooms list` (exit code 0 = server reachable).

**At the start of every session, set the room once.** All commands with `--room` read from `ZNDRAW_ROOM` automatically — never pass `--room` manually after this.

```bash
# Discover rooms
uv run zndraw-cli rooms list

# Set room for all subsequent commands — do this ONCE
export ZNDRAW_ROOM=<room-id>

# All commands now use this room automatically
uv run zndraw-cli frames count          # no --room needed
uv run zndraw-cli extensions list       # no --room needed
uv run zndraw-cli step get              # no --room needed
```

For explicit server connection: `uv run zndraw-cli rooms list --url http://localhost:8000 --token TOKEN`

Env vars: `ZNDRAW_URL`, `ZNDRAW_TOKEN`, `ZNDRAW_ROOM`

## Authentication

The CLI auto-discovers tokens in this order:
1. `--token` flag / `ZNDRAW_TOKEN` env var, OR `--user` / `ZNDRAW_USER` + `--password` / `ZNDRAW_PASSWORD` flags (mutually exclusive)
2. Stored token in `~/.zndraw/tokens.json` (validated on use, removed on 401)
3. Guest token (anonymous fallback)

```bash
uv run zndraw-cli auth status           # check current identity
uv run zndraw-cli auth login            # browser-based login (see below)
uv run zndraw-cli auth logout           # remove stored token
```

For admin users who need to act as another user:

```bash
uv run zndraw-cli auth status           # look for "is_superuser": true
uv run zndraw-cli admin users list
uv run zndraw-cli admin users login USER_ID
```

After `auth login` or `admin users login`, the token is stored in `~/.zndraw/tokens.json` and automatically used for all subsequent CLI commands — no `--token` needed.

## Screenshots, GIFs & Browser Sessions

Screenshots and GIF capture require **two prerequisites**:
1. **You must be logged in as the user** (not as a guest)
2. **The user must have the room open in their browser**

### Auth login is INTERACTIVE — you CANNOT approve it yourself

`auth login` opens a browser window where **the human user must click "Approve"**. This is a device-code flow — the CLI prints a code and waits up to 5 minutes for browser approval.

**You MUST NOT** use Playwright, Selenium, or any browser automation to approve the login — no programmatic workarounds.
- Attempt to create sessions or capture screenshots without proper auth
**Required workflow for screenshots/GIFs:**

```bash
# 1. Check if already logged in
uv run zndraw-cli auth status
# Look for "token_source": "stored" — means you're logged in

# 2. If NOT logged in: start login and ASK THE USER to approve
uv run zndraw-cli auth login
# ⚠️ STOP HERE — tell the user: "Please approve the login in your browser window"
# WAIT for the user to confirm they approved it

# 3. Verify login succeeded
uv run zndraw-cli auth status

# 4. Ensure room is open in user's browser, then capture
uv run zndraw-cli screenshots request           # single screenshot
uv run zndraw-cli gif capture --orbit -o out.gif # orbit GIF
```

**If `auth status` shows `token_source: "stored"`, skip steps 2-3** — you're already logged in.

### GIF capture tips

- **Hide helper geometries before capture** — curves, helper lines, etc. are visible in GIFs. Toggle them off first:
  ```bash
  uv run zndraw-cli geometries toggle my-curve    # hide before capture
  uv run zndraw-cli gif capture --curve my-curve --curve-step 0.02 -o out.gif  # smaller step = more frames
  uv run zndraw-cli geometries toggle my-curve    # show again after
  ```
- **Use presets for visual styles** — don't manually create PathTracing/lighting geometries. Use `presets apply pathtracing` instead.
- **Increase `--delay` for path tracing** — default `0.02s` is too fast for path tracing to converge. Use `--delay 0.5` or higher:
  ```bash
  uv run zndraw-cli presets apply pathtracing
  uv run zndraw-cli gif capture --orbit --delay 0.5 -o pathtraced.gif
  ```
- **`--center` and `--target`** — `--center x,y,z` sets the orbit center, `--target x,y,z` sets where the camera looks. Default: target = center.

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
| `geometries` | `list`, `get`, `set`, `delete`, `toggle`, `set-prop`, `types`, `describe` |
| `preset` | `list`, `get`, `load`, `apply`, `save`, `export`, `delete`, `reset` |
| `screenshots` | `list`, `request`, `get` |
| `gif` | `capture` |
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
- Escape sequences `\n` and `\t` are processed — `send 'line1\nline2'` produces a multi-line message
- Send relevant updates to the chat in addition to printing them.

**Geometry tips:**
- `set KEY --data '{"resolution": 32}'` without `--type` performs a partial update (PATCH, deep merge)
- `set KEY --type Sphere --data '{...}'` with `--type` performs a full replace (PUT)
- `toggle KEY` flips a geometry's `active` state (on/off)
- `set-prop KEY material.opacity 0.5` sets a single property using dot notation (builds nested PATCH)

### Isosurface data contract

`Isosurface` renders a 3D surface from volumetric data stored in a frame's `atoms.info` dict.

- `cube_key`: dotted path to the info key (e.g. `info.orbital_homo`)
- The referenced dict must contain:
  - `grid`: `np.ndarray` of shape `(Nx, Ny, Nz)` — scalar field values
  - `origin`: `np.ndarray` of shape `(3,)` — world-space origin (Angstrom)
  - `cell`: `np.ndarray` of shape `(3, 3)` — axis vectors spanning the grid (Angstrom)
- Parameters: `isovalue` (float, default `0.02`), `resolution` (0–1, default `1.0`), `opacity` (0–1, default `0.6`), `color` (hex, default `#2244CC`)

```python
from zndraw import ZnDraw
from zndraw.geometries import Isosurface
vis = ZnDraw(room="ROOM")
vis.geometries["homo"] = Isosurface(cube_key="info.orbital_homo", isovalue=0.02)
```

Use positive and negative isovalues to show both lobes of an orbital.

## Extensions: Discover → Describe → Run

```bash
# 1. Discover what's available
uv run zndraw-cli extensions list

# 2. Get parameter schema — shows exact --flag names and types
uv run zndraw-cli extensions describe "@internal:modifiers:AddFromSMILES"
# Output tells you: --smiles (str), --optimize (bool), etc.

# 3. Run with --key value flags matching the schema (NOT --params JSON)
uv run zndraw-cli extensions run "@internal:modifiers:AddFromSMILES" --smiles "CCO"

# 3b. Run and block until completed (preferred)
uv run zndraw-cli extensions run "@internal:selections:All" --wait --timeout 30

# 4. Manual polling fallback (if not using --wait)
uv run zndraw-cli jobs status JOB_ID
```

**Important:** Extension parameters are passed as `--flag value` pairs (from `describe` output), NOT as `--params '{"key": "value"}'` JSON.

Three extension categories exist: **modifiers** (edit atoms), **selections** (pick atoms), **analysis** (create plots). Discover which are available via `extensions list`.

**Note:** Custom extensions registered via `vis.register_job()` are scoped to the room they were registered in.

## Visual Presets

Named visual styles (materials, lighting, fog) that can be applied to room geometries. Bundled presets (`matt`, `flat`, `glossy`, `pathtracing`) are always available in every room without needing to be created. Creating a room-level preset with the same name overrides the bundled one; deleting the override makes the bundled version reappear.

```bash
uv run zndraw-cli presets list                                     # list presets
uv run zndraw-cli presets apply matt                               # apply bundled preset
uv run zndraw-cli presets reset                                    # reset to factory defaults
uv run zndraw-cli presets save my-style -p "particles*" -p "fog"   # save current state
uv run zndraw-cli presets export my-style ./my-style.json          # export for sharing
uv run zndraw-cli presets load ./my-style.json                     # import from file
```

Presets use fnmatch patterns (`particles*`, `*light*`, `fog`) and optional `geometry_type` filters to target specific geometries. Rules are deep-merged — only specified config keys are overridden.

## Python Scripting with ZnDraw Client

**Only use Python when no extension covers the task** — typically for custom data reads, math, filtering, or Plotly figures. If you haven't run `extensions list` yet, do that first.

Use `uv run python -c "..."` with the **`ZnDraw` Python client**. It connects to the server and provides a `MutableSequence[ase.Atoms]` interface.

```python
from zndraw import ZnDraw, Extension, Category
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

# Run extensions — returns TaskHandle
task = vis.run("@internal:modifiers:Delete")
task = vis.run("@internal:modifiers:AddFromSMILES", smiles="CCO")
task.wait(timeout=30)            # block until completed/failed
task.status                      # "pending" | "running" | "completed" | "failed"
task.id                          # task ID string

# Discovery
list(vis.extensions)                                    # all extension names
vis.extensions["@internal:modifiers:Delete"]["schema"]   # parameter schema

# Task handles
vis.tasks[task_id]               # TaskHandle (with .wait(), .status, .id)
vis.tasks("running")             # filtered view

# Sessions — room-scoped Mapping of active browser sessions
vis.sessions                     # Mapping[str, Session] (all users in room)
list(vis.sessions)               # list of session SIDs

# Screenshots (requires own browser session — see "Screenshots, GIFs" section)
sids = list(vis.sessions)           # list active session SIDs
session = vis.sessions[sids[0]]     # get a session by SID
img = session.screenshot()          # capture screenshot (own sessions only)
img.data                            # PNG bytes
img.save("frame.png")               # save to file

# Visual Presets — MutableMapping[str, Preset]
list(vis.presets)                 # list preset names
vis.presets.apply("matt")        # apply preset to room geometries
vis.presets.apply("@default")    # reset all geometries to factory defaults
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

### Select atoms by element (Python fallback)

Only use this if `extensions list` shows no suitable selection extension (e.g. no ByElement/BySymbol):

```bash
uv run python -c "
from zndraw import ZnDraw
vis = ZnDraw(room='ROOM')
atoms = vis[vis.step]
vis.selection = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'H']
"
```

### Compute RDF (Python fallback)

Only use if `extensions list` shows no RDF analysis extension:

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
For large files that support sliced access (`.h5`, `.lmdb`, `.zarr`, `.db`), use `mount` — frames are served lazily on demand.

```bash
uv run zndraw-cli mount trajectory.h5                 # new room (auto-generated ID)
uv run zndraw-cli mount trajectory.h5 --room my-room  # specific room
# Output: {"room_id": "...", "url": "http://...", "frame_count": 50000}
```
This command keeps serving frames until Ctrl+C is pressed.


## Camera and Session Control

Each browser session has an active camera. Use session commands to manage cameras:

```bash
uv run zndraw-cli sessions list                          # find SIDs
uv run zndraw-cli sessions get-camera SESSION_ID         # current camera
uv run zndraw-cli sessions set-camera SESSION_ID my-cam  # switch camera
```

Camera geometries can be partially updated like any geometry:
```bash
uv run zndraw-cli geometries set my-camera --data '{"fov": 45}'
```

## Task Mapping

| User request | Approach |
|-------------|----------|
| "build a box / pack molecules" | `extensions list` → find PackBox or similar → `describe` → `run` |
| "delete selected atoms" | `extensions list` → find Delete modifier → `describe` → `run` |
| "add a molecule" | `extensions list` → find AddFromSMILES → `describe` → `run` with params |
| "show me aspirin / molecule" | `extensions list` → find AddFromSMILES → `describe` → `run --smiles "..."` |
| "select all hydrogens" | `extensions list` → check for element selection extension first; fallback: Python `vis.selection = [...]` |
| "duplicate / replicate atoms" | `extensions list` → find Replicate or similar → `describe` → `run` |
| "compute RDF / analysis" | `extensions list` → check for analysis extension first; fallback: Python ASE + Plotly |
| "jump to frame 42" | `uv run zndraw-cli step set 42` |
| "how many frames?" | `uv run zndraw-cli frames count` |
| "bookmark high-energy frames" | Python: iterate `vis`, set `vis.bookmarks[i]` |
| "plot energy over time" | Python: `vis.get(...)` + Plotly → `vis.figures[...]` |
| "visualize this large H5 file" | `uv run zndraw-cli mount trajectory.h5` (lazy) |
| "open / show the room" | `uv run zndraw-cli rooms open` |
| "lock the room" | `uv run zndraw-cli rooms lock` or Python: `vis.locked = True` |
| "send a chat message" | `uv run zndraw-cli chat send "msg"` or Python: `vis.chat.send("msg")` |
| "take a screenshot" | Check `auth status` → `auth login` if needed (ask user to approve in browser) → `screenshots request` |
| "create a GIF / orbit animation" | Check `auth status` → `auth login` if needed (ask user) → ensure room open in browser → `gif capture --orbit -o out.gif` |
| "create a path-traced GIF" | `presets apply pathtracing` → `gif capture --orbit --delay 0.5 -o out.gif` |
| "GIF along a camera path" | Curve must exist as a geometry (created via UI or `geometries set`) → `geometries toggle CURVE` to hide → `gif capture --curve CURVE ...` → toggle CURVE back |
| "list browser sessions" | `uv run zndraw-cli sessions list` |
| "who am I / admin check" | `uv run zndraw-cli auth status` |
| "apply matt style" | `uv run zndraw-cli presets apply matt` or Python: `vis.presets.apply("matt")` |
| "reset to defaults" | `uv run zndraw-cli presets reset` or Python: `vis.presets.apply("@default")` |
| "save current look as preset" | `uv run zndraw-cli presets save my-style -p "particles*" -p "fog"` |
| "share a preset" | `uv run zndraw-cli presets export my-style ./my-style.json` |
| "toggle fog off" | `uv run zndraw-cli geometries toggle fog` |
| "make particles transparent" | `uv run zndraw-cli geometries set-prop particles material.opacity 0.5` |
| "run extension and wait" | `extensions run EXT --wait` or Python: `vis.run(EXT).wait()` |
| "visualize volumetric data / orbitals" | `geometries set KEY --type Isosurface --data '{"cube_key": "info.KEY", "isovalue": 0.02}'` |
| "what can I do?" | `extensions list` + `--help` on each resource group |

## Common Mistakes

- **Writing Python before checking extensions** — the #1 mistake. An extension likely already does what you're about to script. `extensions list` takes 1 second; writing+debugging a Python script takes minutes. Always check first.
- **Trying to automate `auth login` approval** — `auth login` opens a browser for the **human user** to approve. You CANNOT use Playwright, Selenium, or any automation. Ask the user to approve in their browser and wait for confirmation.
- **Screenshots/GIFs without login** — guest tokens create a different user identity than the browser session. You must `auth login` first, and the user must have the room open in their browser.
- **Repeating `--room` on every command** — set `export ZNDRAW_ROOM=...` once at session start. Never pass `--room` manually after that.
- **Running bare `zndraw-cli`** — always use `uv run zndraw-cli`; the binary is not on PATH
- **Guessing extension names/params** — extensions change per server and per room. Always `list` then `describe` first; never assume you know what's available.
- **Fetching everything** — use `vis.get(..., keys=[...])` or CLI `--keys` to limit data
- **Not waiting for tasks** — use `--wait` flag or Python `task.wait()` instead of manual polling
- **Wrong room** — `rooms list` to discover, don't guess UUIDs
- **Piping CLI to Python for mutations** — use the `ZnDraw` client instead; it connects directly and handles serialization
- **Leaving helper geometries visible in GIFs** — curves, helper lines show up in captures. Toggle them off with `geometries toggle KEY` before `gif capture`.
- **Manual path tracing setup** — use `presets apply pathtracing`, not manual `geometries set --type PathTracing`. Presets are single-command and correct.
- **Default `--delay` with path tracing** — `0.02s` produces noisy frames. Use `--delay 0.5` or higher for path tracing to accumulate samples.
- **`screenshots get` instead of `screenshots request`** — `get ID` downloads an existing screenshot by its database ID. `request` captures a new one.
