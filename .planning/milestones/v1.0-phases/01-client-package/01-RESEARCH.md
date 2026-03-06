# Phase 1: Client Package - Research

**Researched:** 2026-03-05
**Domain:** Python module decomposition / package refactoring
**Confidence:** HIGH

## Summary

Phase 1 converts the monolithic `src/zndraw/client.py` (2298 lines, ~81 KB) into a `src/zndraw/client/` package with six focused modules. The file has clear section boundaries (marked by `# ===` comment headers) that map directly to the target modules: serialization helpers, exceptions, lock, API manager, socket manager, and the main ZnDraw class.

The primary risk is breaking existing imports. Analysis of `src/` and `tests/` reveals 30+ import sites that reference `from zndraw.client import ...`. Python's `__init__.py` re-export mechanism makes this safe: the package `__init__.py` re-exports all previously-public names, so `from zndraw.client import ZnDraw` continues to resolve without downstream changes. The one exception is `RoomLockedError` and `ZnDrawError`, which per user decision move to `src/zndraw/exceptions.py` (shared), while `NotConnectedError` stays in `client/exceptions.py`.

This is a pure structural refactor with no functional changes. The existing 499 integration tests serve as the regression gate. New unit tests target only the extracted serialization helpers and exception hierarchy.

**Primary recommendation:** Use Python's package `__init__.py` as a compatibility shim -- re-export every name that was previously importable from `zndraw.client`, then update internal imports to point directly at the new submodules.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- `ZnDrawError` and `RoomLockedError` move to `src/zndraw/exceptions.py` (shared between server and client code; `RoomLockedError` is already referenced by server-side `raise_for_client()`)
- `NotConnectedError` goes in `client/exceptions.py` (client-only, imports `ZnDrawError` from `zndraw.exceptions`)
- RFC 9457 `ProblemType` classes stay in `exceptions.py` -- no mixing concerns, just colocating the shared base exception
- `_estimate_frame_size` goes in `serialization.py` (operates on frame data structures, same domain as other serialization helpers)
- Other module assignments per REQUIREMENTS: `serialization.py`, `lock.py`, `api.py`, `socket.py`, `core.py`
- Behavioral tests for serialization helpers: round-trip tests with real ASE Atoms objects created via `molify.smiles2atoms`
- Smoke tests for exception hierarchy (correct inheritance, instantiation)
- No mocked unit tests for APIManager/SocketManager -- existing 499 integration tests cover that
- Test location: `tests/test_client/` directory (mirrors `tests/test_cli_agent/` pattern)

### Claude's Discretion
- `client/__init__.py` re-export surface (what's importable from `zndraw.client` directly)
- Import migration strategy for internal consumers (clean break, update imports in-place)
- Circular import resolution between client submodules
- Exact module placement for any remaining borderline helpers

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| CLNT-01 | `client.py` converted to `src/zndraw/client/` package with `__init__.py` | Package structure pattern verified in codebase (`storage/`, `extensions/`, `cli_agent/`); `__init__.py` re-export pattern documented below |
| CLNT-02 | Serialization helpers extracted to `client/serialization.py` | Lines 79-163 identified: `atoms_to_json_dict`, `json_dict_to_atoms`, `raw_frame_to_atoms`, `_estimate_frame_size`, `_TARGET_CHUNK_BYTES`, `_MAX_CHUNK_FRAMES`; import dependencies documented |
| CLNT-03 | Exception classes extracted to `client/exceptions.py` | User decision: `ZnDrawError` + `RoomLockedError` go to shared `src/zndraw/exceptions.py`; only `NotConnectedError` in `client/exceptions.py` |
| CLNT-04 | `ZnDrawLock` extracted to `client/lock.py` | Lines 188-246; depends on `APIManager` (type hint only at module level); imports `EDIT_LOCK_REFRESH` lazily from `zndraw.schemas` |
| CLNT-05 | `APIManager` extracted to `client/api.py` | Lines 254-1136; self-contained dataclass; references `ZnDrawError`, `RoomLockedError`, `PROBLEM_TYPES`, `ProblemDetail` |
| CLNT-06 | `SocketManager` extracted to `client/socket.py` | Lines 1143-1263; references `ZnDraw` (TYPE_CHECKING), `SyncClientWrapper`, socket events |
| CLNT-07 | `ZnDraw` main class extracted to `client/core.py` | Lines 1270-2298; references all other client modules; largest single class |
| CLNT-08 | `from zndraw import ZnDraw` continues to work | Chain: `zndraw/__init__.py` -> `zndraw.client` -> `client/__init__.py` re-exports `ZnDraw` from `client.core` |
| CLNT-09 | All 499 existing tests pass unchanged | `__init__.py` re-exports all names currently imported via `from zndraw.client import X`; full import inventory documented below |
| CLNT-10 | Unit tests added for extracted modules | Serialization round-trip tests with `molify.smiles2atoms`; exception hierarchy smoke tests; test location: `tests/test_client/` |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python | >=3.11 | Package system with `__init__.py` | Native feature, no dependencies |
| pytest | 9.0.2 | Test framework | Already in dev dependencies |
| ase | >=3.27.0 | Atoms objects for serialization tests | Already a project dependency |
| molify | (available) | Create test Atoms fixtures via `smiles2atoms` | User preference per CONTEXT.md |
| asebytes | >=0.3.0a3 | `encode`/`decode` for frame serialization | Already used in `client.py` |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| ruff | (dev) | Format + lint after refactoring | Run `uv run ruff format .` and `uv run ruff check --select I --fix .` after every module creation |
| pyright | (dev) | Type checking | Run `uv run pyright .` to verify no type regressions |

### Alternatives Considered
None -- this is a pure refactoring phase using only existing project tooling.

**Installation:**
No new dependencies needed.

## Architecture Patterns

### Recommended Project Structure
```
src/zndraw/client/
    __init__.py          # Re-exports all public names for backwards compatibility
    serialization.py     # atoms_to_json_dict, json_dict_to_atoms, raw_frame_to_atoms,
                         #   _estimate_frame_size, _TARGET_CHUNK_BYTES, _MAX_CHUNK_FRAMES
    exceptions.py        # NotConnectedError only (ZnDrawError + RoomLockedError -> zndraw.exceptions)
    lock.py              # ZnDrawLock dataclass
    api.py               # APIManager dataclass (all REST API methods)
    socket.py            # SocketManager dataclass (Socket.IO connection + handlers)
    core.py              # ZnDraw(MutableSequence[ase.Atoms]) main class

tests/test_client/
    __init__.py          # Empty (needed for pytest discovery in subdirectory)
    test_serialization.py  # Round-trip tests for serialization helpers
    test_exceptions.py     # Exception hierarchy smoke tests
```

### Pattern 1: Package `__init__.py` Re-Export (Compatibility Shim)

**What:** The `client/__init__.py` re-exports every name that was previously importable from `zndraw.client`, preserving all existing import paths.

**When to use:** Always, for this refactoring. This is how `storage/__init__.py` and `extensions/__init__.py` already work in this codebase.

**Example:**
```python
# src/zndraw/client/__init__.py
"""ZnDraw Python client package."""

from zndraw.client.api import APIManager
from zndraw.client.core import ZnDraw
from zndraw.client.exceptions import NotConnectedError
from zndraw.client.lock import ZnDrawLock
from zndraw.client.serialization import (
    _estimate_frame_size,
    _TARGET_CHUNK_BYTES,
    _MAX_CHUNK_FRAMES,
    atoms_to_json_dict,
    json_dict_to_atoms,
    raw_frame_to_atoms,
)
from zndraw.client.socket import SocketManager

# Re-export shared exceptions so `from zndraw.client import ZnDrawError` still works
from zndraw.exceptions import ZnDrawError, RoomLockedError  # noqa: F401 (re-export)

# Re-export accessors that were previously importable via `from zndraw.client import Sessions`
from zndraw.accessors import Sessions  # noqa: F401 (re-export)

__all__ = [
    "APIManager",
    "NotConnectedError",
    "RoomLockedError",
    "Sessions",
    "SocketManager",
    "ZnDraw",
    "ZnDrawError",
    "ZnDrawLock",
    "_TARGET_CHUNK_BYTES",
    "_MAX_CHUNK_FRAMES",
    "_estimate_frame_size",
    "atoms_to_json_dict",
    "json_dict_to_atoms",
    "raw_frame_to_atoms",
]
```

### Pattern 2: Shared Exception Placement

**What:** `ZnDrawError` and `RoomLockedError` move to `src/zndraw/exceptions.py` (the shared exceptions module). This eliminates the lazy import `from zndraw.client import RoomLockedError` that currently exists in `RoomLocked.raise_for_client()`.

**When to use:** Per user decision. These exceptions are referenced by both client code and server-side exception handlers.

**Example:**
```python
# In src/zndraw/exceptions.py, add at top (after existing imports):
class ZnDrawError(Exception):
    """Base exception for ZnDraw client errors."""

class RoomLockedError(ZnDrawError):
    """Raised when the room is locked (admin lock or edit lock by another user)."""

# Then RoomLocked.raise_for_client() changes from:
#   from zndraw.client import RoomLockedError
#   raise RoomLockedError(...)
# to simply:
#   raise RoomLockedError(...)
```

### Pattern 3: TYPE_CHECKING for Cross-Module References

**What:** Use `TYPE_CHECKING` blocks to avoid circular imports between client submodules. Several cross-references exist (e.g., `SocketManager` references `ZnDraw`, `ZnDrawLock` references `APIManager`).

**When to use:** When module A needs module B's type for annotations but module B imports from module A at runtime.

**Example:**
```python
# client/socket.py
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from zndraw.client.core import ZnDraw

# At runtime, ZnDraw is only used as a type annotation.
# The dataclass field `zndraw: ZnDraw` works because of `from __future__ import annotations`.
```

### Anti-Patterns to Avoid
- **Do NOT use underscore-prefixed module names** (`_serialization.py`): The early research suggested `_serialization.py` etc., but the codebase convention for packages like `storage/`, `extensions/`, `cli_agent/` uses plain names without underscores. Follow the existing convention.
- **Do NOT re-export names that were never publicly imported**: Only re-export what the import inventory shows is actually used by external code.
- **Do NOT update test imports**: The `__init__.py` re-exports ensure all existing `from zndraw.client import X` paths continue to work. Tests must pass unchanged (CLNT-09).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Import compatibility | Custom import hooks or `sys.modules` tricks | `__init__.py` re-exports | Standard Python mechanism, used elsewhere in this codebase |
| Circular import resolution | Runtime import workarounds | `TYPE_CHECKING` + `from __future__ import annotations` | Clean, standard, already used throughout `client.py` |
| Test Atoms fixtures | Manual `ase.Atoms(symbols=..., positions=...)` | `molify.smiles2atoms("O")` | User preference; produces realistic molecules with proper connectivity |

## Common Pitfalls

### Pitfall 1: Missing Re-Export Breaks External Import
**What goes wrong:** A name that was importable from `from zndraw.client import X` stops working because it wasn't added to the `__init__.py` re-exports.
**Why it happens:** The import inventory is incomplete. Some imports are buried inside test functions (e.g., `from zndraw.client import _estimate_frame_size` inside `TestEstimateFrameSize`).
**How to avoid:** Use the complete import inventory below. Every name must be re-exported.
**Warning signs:** Any test importing from `zndraw.client` fails with `ImportError`.

### Pitfall 2: Circular Import Between Client Submodules
**What goes wrong:** `core.py` imports from `api.py`, `socket.py`, `serialization.py`, `lock.py`; but `socket.py` references `ZnDraw` type. If these are runtime imports, Python raises `ImportError`.
**Why it happens:** The original monolithic file had no circular import issue because everything was in one module. Splitting creates potential cycles.
**How to avoid:** Use `TYPE_CHECKING` + `from __future__ import annotations` consistently. The current `client.py` already uses this pattern for `FrameSource`, `SessionItem`, `ZnDrawTqdm`. Apply the same approach for cross-submodule type references.
**Warning signs:** `ImportError` at import time mentioning partially initialized module.

**Specific circular dependencies to handle:**
- `socket.py` has `zndraw: ZnDraw` field -> use `TYPE_CHECKING` for `ZnDraw`
- `lock.py` has `api: APIManager` field -> safe, `APIManager` is in `api.py` which has no back-reference to `lock.py`
- `core.py` imports all other modules -> safe, it's the leaf; no other submodule imports from `core.py` at runtime

### Pitfall 3: `ZnDrawError` / `RoomLockedError` Import Cycle with `exceptions.py`
**What goes wrong:** Moving `ZnDrawError` and `RoomLockedError` to `src/zndraw/exceptions.py` could create a cycle: `exceptions.py` imports from `client` (current line 336), and `client` imports from `exceptions`.
**Why it happens:** Currently `RoomLocked.raise_for_client()` has `from zndraw.client import RoomLockedError`. After the move, this becomes a local reference, eliminating the cycle.
**How to avoid:** Move the exception classes FIRST, then update `RoomLocked.raise_for_client()` to reference them directly (no lazy import needed). The `api.py` module imports `ZnDrawError` and `RoomLockedError` from `zndraw.exceptions` -- this is already the direction of the dependency (client -> exceptions), not the reverse.

### Pitfall 4: `Sessions` Re-Export
**What goes wrong:** `test_sessions.py` line 12 has `from zndraw.client import Sessions`. The `Sessions` class is defined in `accessors.py`, not `client.py`, but was importable because `client.py` did `from zndraw.accessors import Sessions`. After the split, this import path must still work.
**Why it happens:** The monolithic `client.py` had a top-level import of `Sessions` which made it importable as `zndraw.client.Sessions`.
**How to avoid:** Include `Sessions` in the `__init__.py` re-exports.

### Pitfall 5: `atoms_to_json_dict` Import in Test
**What goes wrong:** `test_routes_frames.py` line 24 has `from zndraw.client import atoms_to_json_dict`. Must be re-exported.
**Why it happens:** Test uses the serialization helper directly.
**How to avoid:** Include in `__init__.py` re-exports. Already listed in the inventory.

### Pitfall 6: Old `client.py` File Left Behind
**What goes wrong:** If the old `client.py` file is not deleted, Python may import from it instead of the new `client/` package (file takes precedence over package in some configurations).
**Why it happens:** Python module resolution: a `.py` file and a directory with the same name cannot coexist.
**How to avoid:** Delete `client.py` as part of creating the `client/` directory. This is actually enforced by Python -- you cannot have both `client.py` and `client/` in the same directory. The deletion is step 1, not a cleanup step.

## Code Examples

### Complete Import Inventory (Verified from Codebase)

Every `from zndraw.client import X` found in `src/` and `tests/`:

**Source code (`src/`):**
| Import | File | Line | Context |
|--------|------|------|---------|
| `ZnDraw` | `__init__.py` | 3 | Top-level re-export |
| `ZnDraw` | `cli.py` | 27 | Top-level import |
| `ZnDraw` | `executor.py` | 58 | Lazy import inside function |
| `ZnDraw` | `cli_agent/mount.py` | 25 | TYPE_CHECKING import |
| `ZnDraw` | `tqdm.py` | 11 | TYPE_CHECKING import |
| `RoomLockedError` | `exceptions.py` | 336 | Lazy import in `raise_for_client()` |
| `RoomLockedError` | `cli_agent/connection.py` | 359 | Lazy import inside function |
| `RoomLockedError, ZnDrawError` | `cli_agent/connection.py` | 368 | Lazy import inside function |
| `APIManager` | `accessors.py` | 31 | TYPE_CHECKING import |
| `APIManager` | `cli_agent/rooms.py` | 41 | Lazy import inside function |

**Tests (`tests/`):**
| Import | File | Line |
|--------|------|------|
| `ZnDraw` | `test_client.py` | 13 (via `from zndraw import ZnDraw`) |
| `_estimate_frame_size` | `test_client.py` | 36, 42, 48 |
| `ZnDrawError` | `test_client.py` | 139, 153 |
| `RoomLockedError, ZnDraw` | `test_geometry_ownership.py` | 12 |
| `RoomLockedError` | `test_edit_lock_integration.py` | 13 |
| `Sessions` | `test_sessions.py` | 12 |
| `atoms_to_json_dict` | `test_routes_frames.py` | 24 |
| `ZnDraw` | `test_default_camera.py` | 17 |
| `ZnDraw` | `test_client_source.py` | 10 |
| `ZnDraw` | `test_cli_agent/test_mount.py` | 11 |
| `ZnDraw` | `test_client_token_discovery.py` | 27, 59, 79 |

### Module Dependency Graph (Runtime Imports Only)

```
serialization.py  -> ase, base64, asebytes (encode, decode), zndraw.enrichment, zndraw.connectivity
exceptions.py     -> zndraw.exceptions (imports ZnDrawError)
lock.py           -> threading, logging, api.py (APIManager), zndraw.schemas (EDIT_LOCK_REFRESH)
api.py            -> httpx, json, msgpack, zndraw.exceptions (PROBLEM_TYPES, ProblemDetail, ZnDrawError*, RoomLockedError*)
socket.py         -> socketio, logging, zndraw_socketio, zndraw.socket_events
core.py           -> api.py, socket.py, lock.py, serialization.py, exceptions.py,
                     zndraw.accessors, zndraw.geometries.camera, zndraw_joblib, uuid, etc.
```

*Note:* `ZnDrawError` and `RoomLockedError` will be in `zndraw.exceptions` after the move, so `api.py` imports them from there.

### Serialization Module Content (Exact Boundaries)

```python
# client/serialization.py - lines 79-163 of current client.py
# Top-level imports needed:
import base64
from typing import Any
import ase
from asebytes import decode, encode
from zndraw.enrichment import add_colors, add_radii

# Functions to extract:
def atoms_to_json_dict(atoms, connectivity_threshold=1000): ...
def json_dict_to_atoms(data): ...
def raw_frame_to_atoms(frame): ...

_TARGET_CHUNK_BYTES = 2_000_000
_MAX_CHUNK_FRAMES = 1000

def _estimate_frame_size(frame): ...
```

### Test Fixtures Pattern

```python
# tests/test_client/test_serialization.py
import pytest
from molify import smiles2atoms

from zndraw.client.serialization import (
    atoms_to_json_dict,
    json_dict_to_atoms,
    _estimate_frame_size,
)


@pytest.mark.parametrize("smiles", ["O", "CCO", "c1ccccc1"])
def test_round_trip(smiles: str) -> None:
    """Atoms survive encode -> decode round trip."""
    atoms = smiles2atoms(smiles)
    json_dict = atoms_to_json_dict(atoms)
    recovered = json_dict_to_atoms(json_dict)
    assert len(recovered) == len(atoms)
    # positions should match within tolerance
    assert recovered.positions == pytest.approx(atoms.positions, abs=1e-6)
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Monolithic `client.py` (2298 lines) | `client/` package with 6 modules | This phase | Each module has single responsibility |
| `RoomLockedError` in client, lazy-imported by `exceptions.py` | `RoomLockedError` in `exceptions.py`, imported by client | This phase | Eliminates circular lazy import |

**Deprecated/outdated:**
- `vis.log()` and `vis.progress_bar()` are already deprecated in `core.py` -- carry them over as-is, do not remove (out of scope per REQUIREMENTS)

## Open Questions

1. **`_decode_raw_frame` staticmethod placement**
   - What we know: `ZnDraw._decode_raw_frame()` is a static method on the ZnDraw class (line 2131). It uses `msgpack` and `msgpack_numpy`. It could go in `serialization.py` or stay in `core.py`.
   - What's unclear: Whether it conceptually belongs with the other serialization helpers or is tightly coupled to the `ZnDraw.get()` method.
   - Recommendation: Keep it in `core.py` as a static method since it's only used by `ZnDraw.get()` and its import footprint (`msgpack_numpy`) differs from the other serialization helpers. This is within Claude's discretion per CONTEXT.md.

2. **`log` logger variable**
   - What we know: `log = logging.getLogger(__name__)` exists at module level in `client.py`. After the split, each module should have its own logger.
   - What's unclear: Whether the logger name `zndraw.client` matters for any log filtering configuration.
   - Recommendation: Each module uses `logging.getLogger(__name__)` which gives names like `zndraw.client.core`, `zndraw.client.api`, etc. This follows Python convention and is strictly more specific than the old single logger.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | pytest 9.0.2 |
| Config file | None (convention-based discovery via `pyproject.toml`) |
| Quick run command | `uv run pytest tests/test_client/ -x -q` |
| Full suite command | `uv run pytest tests/ -x` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| CLNT-01 | Package exists with `__init__.py` | smoke | `uv run python -c "from zndraw.client import ZnDraw"` | N/A (import check) |
| CLNT-02 | Serialization helpers work | unit | `uv run pytest tests/test_client/test_serialization.py -x` | Wave 0 |
| CLNT-03 | Exception hierarchy correct | unit | `uv run pytest tests/test_client/test_exceptions.py -x` | Wave 0 |
| CLNT-04 | `ZnDrawLock` importable | smoke | `uv run python -c "from zndraw.client.lock import ZnDrawLock"` | N/A |
| CLNT-05 | `APIManager` importable | smoke | `uv run python -c "from zndraw.client.api import APIManager"` | N/A |
| CLNT-06 | `SocketManager` importable | smoke | `uv run python -c "from zndraw.client.socket import SocketManager"` | N/A |
| CLNT-07 | `ZnDraw` in `core.py` | smoke | `uv run python -c "from zndraw.client.core import ZnDraw"` | N/A |
| CLNT-08 | `from zndraw import ZnDraw` works | regression | `uv run python -c "from zndraw import ZnDraw; print(ZnDraw)"` | N/A |
| CLNT-09 | All 499 existing tests pass | regression | `uv run pytest tests/ -x` | Existing |
| CLNT-10 | Unit tests for extracted modules | unit | `uv run pytest tests/test_client/ -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `uv run pytest tests/test_client/ -x -q && uv run python -c "from zndraw import ZnDraw"`
- **Per wave merge:** `uv run pytest tests/ -x`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `tests/test_client/__init__.py` -- empty file for pytest discovery
- [ ] `tests/test_client/test_serialization.py` -- covers CLNT-02, CLNT-10
- [ ] `tests/test_client/test_exceptions.py` -- covers CLNT-03, CLNT-10
- No framework install needed -- pytest 9.0.2 already available

## Sources

### Primary (HIGH confidence)
- Direct codebase analysis of `src/zndraw/client.py` (2298 lines, every line read)
- Direct codebase analysis of all files importing from `zndraw.client` (grep of `src/` and `tests/`)
- Existing package patterns: `src/zndraw/storage/__init__.py`, `src/zndraw/extensions/__init__.py`, `src/zndraw/cli_agent/__init__.py`
- `src/zndraw/exceptions.py` -- full read, verified `RoomLocked.raise_for_client()` lazy import on line 336
- `.planning/phases/01-client-package/01-CONTEXT.md` -- user decisions
- `.planning/REQUIREMENTS.md` -- requirement definitions
- `tests/conftest.py` -- test infrastructure patterns

### Secondary (MEDIUM confidence)
- `.planning/research/ARCHITECTURE.md` -- prior project-level research (aligned with findings; diverges on underscore-prefixed module names -- this research follows codebase convention of plain names)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - no new libraries needed, all tools already in project
- Architecture: HIGH - direct codebase analysis, patterns verified from existing packages
- Pitfalls: HIGH - every import path verified via grep, circular dependency analysis done on actual code
- Validation: HIGH - pytest infrastructure already exists, test patterns verified

**Research date:** 2026-03-05
**Valid until:** Indefinite (pure refactoring, no external dependency changes)
