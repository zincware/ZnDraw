---
name: zndraw
description: Use when interacting with a ZnDraw molecular visualization server — reading frames, navigating trajectories, selecting atoms, bookmarking, creating plots, running extensions, or performing computational chemistry analysis on atomic structures
---

# ZnDraw Agent Skill

## Overview

ZnDraw is an interactive visualization platform for atomistic simulations. Two interaction tools:

1. **`zndraw-cli`** — Typer CLI, structured stdout. For CRUD on rooms, frames, selections, bookmarks, figures, extensions, ... .
2. **`uv run python -c "..."`** — For computation requiring numpy/ASE/plotly. **Must use the `ZnDraw` Python client** to connect to the server and modify data. Use `uv run --with <dependencies>` to include extra packages.

## Discover First, Never Guess

The CLI is self-documenting. **Always discover before acting:**

```bash
zndraw-cli --help                            # all resource groups
zndraw-cli frames --help                     # verbs + options for frames
zndraw-cli extensions list ROOM              # available extensions (fully qualified names)
zndraw-cli extensions describe ROOM NAME     # parameter schema
zndraw-cli frames get ROOM 0                 # inspect one frame to learn the data shape
```

Extension names are fully qualified (e.g. `@internal:modifier:Delete`). Never hardcode names — they vary per server. Run `extensions list` then `extensions describe` to learn parameters before `extensions run`.

## Connection

Both `zndraw-cli` and the Python client require a **running ZnDraw server**. Start one with `uv run zndraw` if none is running. Verify with `zndraw-cli rooms list` (exit code 0 = server reachable).

```bash
# Auto-discovers running server + creates guest token
zndraw-cli rooms list

# Explicit
zndraw-cli --url http://localhost:1234 --token TOKEN rooms list

# Env vars: ZNDRAW_URL, ZNDRAW_TOKEN
# When working on a single room, export ZNDRAW_ROOM to avoid repeating it
export ZNDRAW_ROOM=my-room-id
```

## CLI Quick Reference

Use `<command> --help` for full options. Key patterns:

| Resource | Verbs |
|----------|-------|
| `rooms` | `list`, `create`, `info` |
| `frames` | `count`, `get`, `list`, `extend`, `export`, `delete` |
| `step` | `get`, `set` |
| `selection` | `get`, `set`, `clear` |
| `selection-groups` | `list`, `get`, `set`, `delete` |
| `bookmarks` | `list`, `set`, `delete` |
| `figures` | `list`, `get`, `set`, `delete` |
| `chat` | `list`, `send` |
| `extensions` | `list`, `describe`, `run` |
| `jobs` | `list`, `status` |
| `geometries` | `list`, `get`, `set`, `delete` |
| `screenshots` | `list`, `request`, `get` |
| `sessions` | `list`, `settings` |
| Standalone | `mount FILE` |

**Frames tips:**
- `--format xyz` for extxyz output, default is JSON
- `--keys positions,info.energy` for server-side field filtering
- `--start N --stop M` for ranges (don't fetch 10k frames at once)

** Chat tips:**
- `send` has full markdown support
- send relevant updates to the chat in addtition to printing them.
## Extensions: Discover → Describe → Run → Poll

```bash
# 1. Discover what's available
zndraw-cli extensions list ROOM

# 2. Get parameter schema for the extension you need
zndraw-cli extensions describe ROOM "@internal:modifier:AddFromSMILES"

# 3. Run with --key value flags
#    Returns {"job_id": "abc123", "status": "pending"}
zndraw-cli extensions run ROOM "@internal:modifier:AddFromSMILES" --smiles "CCO"

# 4. Poll for completion (job_id is global, no room needed)
zndraw-cli jobs status JOB_ID
```

Three extension categories exist: **modifiers** (edit atoms), **selections** (pick atoms), **analysis** (create plots). Discover which are available via `extensions list`.

## Python Scripting with ZnDraw Client

When the CLI alone isn't enough (filtering, math, custom analysis, Plotly figures), use `uv run python -c "..."` with the **`ZnDraw` Python client**. It connects to the server and provides a `MutableSequence[ase.Atoms]` interface.

```python
from zndraw import ZnDraw
vis = ZnDraw(url="http://localhost:1234", room="ROOM")

len(vis)              # frame count
vis[0]                # frame as ase.Atoms
vis[vis.step]         # current frame
vis.step = 42         # jump to frame
vis.append(atoms)     # append frame
del vis[0]            # delete frame
vis[:]                # all frames as list[ase.Atoms]

vis.selection = [0, 1, 2]       # set selected atom indices
vis.bookmarks[42] = "label"     # MutableMapping[int, str]
vis.figures["energy"] = fig     # MutableMapping[str, plotly.Figure]
vis.log("message")              # send chat message

# Efficient partial reads (server-side filtering)
vis.get(slice(None), keys=["info.energy"])
```

### Bookmark frames matching a condition

```bash
uv run python -c "
from zndraw import ZnDraw
vis = ZnDraw(url='http://localhost:1234', room='ROOM')
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
vis = ZnDraw(url='http://localhost:1234', room='ROOM')
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
vis = ZnDraw(url='http://localhost:1234', room='ROOM')
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
vis = ZnDraw(url='http://localhost:1234', room='ROOM')
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
zndraw-cli mount trajectory.h5

# Mount into a new specific room
zndraw-cli mount trajectory.h5 --room my-room-id

# Output: {"room_id": "...", "url": "http://...", "frame_count": 50000}
```
This command keeps serving frames until Ctrl+C is pressed.


## Task Mapping

| User request | Approach |
|-------------|----------|
| "jump to frame 42" | `zndraw-cli step set ROOM 42` |
| "how many frames?" | `zndraw-cli frames count ROOM` |
| "select all hydrogens" | Python: `vis.selection = [i for i, s in ...]` |
| "bookmark high-energy frames" | Python: iterate `vis`, set `vis.bookmarks[i]` |
| "plot energy over time" | Python: `vis.get(...)` + Plotly → `vis.figures[...]` |
| "compute RDF / analysis" | Python: `vis[:]` + ASE Analysis → `vis.figures[...]` |
| "delete selected atoms" | `extensions describe` the Delete modifier → `extensions run` |
| "add a molecule" | `extensions describe` AddFromSMILES → `extensions run` with params |
| "visualize this large H5 file" | `zndraw-cli mount trajectory.h5` (lazy) |
| "what can I do?" | `extensions list ROOM` + `--help` on each resource group |

## Common Mistakes

- **Guessing extension names/params** — always `list` then `describe` first
- **Fetching everything** — use `vis.get(..., keys=[...])` or CLI `--keys` to limit data
- **Not polling jobs** — `extensions run` returns a `job_id`; check `jobs status JOB_ID` for completion
- **Wrong room** — `rooms list` to discover, don't guess UUIDs
- **Piping CLI to Python for mutations** — use the `ZnDraw` client instead; it connects directly and handles serialization
