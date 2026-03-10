# Phase 1: Client Package - Context

**Gathered:** 2026-03-05
**Status:** Ready for planning

<domain>
## Phase Boundary

Convert the monolithic `client.py` (2297 lines) into a `src/zndraw/client/` package with single-responsibility modules. Pure structural refactor — no functional changes. All 499 existing tests must pass unchanged.

</domain>

<decisions>
## Implementation Decisions

### Exception placement
- `ZnDrawError` and `RoomLockedError` move to `src/zndraw/exceptions.py` (shared between server and client code; `RoomLockedError` is already referenced by server-side `raise_for_client()`)
- `NotConnectedError` goes in `client/exceptions.py` (client-only, imports `ZnDrawError` from `zndraw.exceptions`)
- RFC 9457 `ProblemType` classes stay in `exceptions.py` — no mixing concerns, just colocating the shared base exception

### Module boundaries
- `_estimate_frame_size` goes in `serialization.py` (operates on frame data structures, same domain as other serialization helpers)
- Other module assignments per REQUIREMENTS: `serialization.py`, `lock.py`, `api.py`, `socket.py`, `core.py`

### Unit test scope
- Behavioral tests for serialization helpers: round-trip tests with real ASE Atoms objects created via `molify.smiles2atoms`
- Smoke tests for exception hierarchy (correct inheritance, instantiation)
- No mocked unit tests for APIManager/SocketManager — existing 499 integration tests cover that
- Test location: `tests/test_client/` directory (mirrors `tests/test_cli_agent/` pattern)

### Claude's Discretion
- `client/__init__.py` re-export surface (what's importable from `zndraw.client` directly)
- Import migration strategy for internal consumers (clean break, update imports in-place)
- Circular import resolution between client submodules
- Exact module placement for any remaining borderline helpers

</decisions>

<code_context>
## Existing Code Insights

### Reusable Assets
- `src/zndraw/exceptions.py`: RFC 9457 ProblemType base with `raise_for_client()` pattern — will receive `ZnDrawError` + `RoomLockedError`
- `src/zndraw/__init__.py`: Already re-exports `from zndraw.client import ZnDraw` — update to `from zndraw.client.core import ZnDraw` or `from zndraw.client import ZnDraw`
- `src/zndraw/accessors.py`: Accessor classes (Selections, Geometries, etc.) imported by client.py — not part of this split

### Established Patterns
- Package-as-directory with `__init__.py` barrel: `geometries/`, `extensions/`, `storage/`, `cli_agent/`
- Test subdirectory pattern: `tests/test_cli_agent/` mirrors `src/zndraw/cli_agent/`
- Clean break on internal imports: PROJECT.md confirms no external consumers of internal classes

### Integration Points
- `src/zndraw/__init__.py` line 3: `from zndraw.client import ZnDraw` — must keep working (CLNT-08)
- `src/zndraw/exceptions.py` line 336: `from zndraw.client import RoomLockedError` — will change to local reference
- `src/zndraw/executor.py` line 58: `from zndraw.client import ZnDraw` — update to new path
- `src/zndraw/cli.py` line 27: `from zndraw.client import ZnDraw` — update to new path
- `src/zndraw/cli_agent/connection.py`: imports `ZnDraw`, `RoomLockedError`, `ZnDrawError` — update all

### Current class/function layout in client.py
- Lines ~79-168: Serialization helpers (`atoms_to_json_dict`, `json_dict_to_atoms`, `raw_frame_to_atoms`, `_estimate_frame_size`)
- Lines ~170-186: Exception classes (`ZnDrawError`, `NotConnectedError`, `RoomLockedError`)
- Lines ~188-252: `ZnDrawLock` class
- Lines ~254-1142: `APIManager` class
- Lines ~1144-1269: `SocketManager` class
- Lines ~1271+: `ZnDraw(MutableSequence[ase.Atoms])` class

</code_context>

<specifics>
## Specific Ideas

- Use `molify.smiles2atoms` for creating test Atoms fixtures (user preference over manual Atoms construction)
- Follow existing test subdirectory pattern (`tests/test_client/test_serialization.py`, etc.)

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-client-package*
*Context gathered: 2026-03-05*
