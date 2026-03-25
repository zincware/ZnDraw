# Plan: Restore Deprecated `progress_tracker`

## Goal
Restore the `v0.6.0` `vis.progress_tracker(...)` method as a deprecated shim to ensure backward compatibility for applications like the SiMGen demo. Also, fix outdated documentation that still references the old API instead of the new `ZnDrawTqdm` class.

## Strategy: Standalone Shim vs. `ZnDrawTqdm` Wrapper
We will use a **Standalone Shim** (`_DeprecatedProgressTracker`) rather than trying to subclass/wrap the new `ZnDrawTqdm`.
- **Why?** The old `v0.6.0` API used absolute percentage updates (`tracker.update(progress=50)`). The new `ZnDrawTqdm` inherits from standard `tqdm`, where `.update(n)` expects a *relative increment*. Hacking `tqdm` to accept absolute percentages via a `progress` kwarg is brittle and confusing. A lightweight standalone shim that directly calls the `self.vis.api.progress_*` REST methods perfectly replicates the `v0.6.0` behavior in just 20 lines of code without polluting the new codebase.

## Steps

### 1. Implement Deprecated Shim in `src/zndraw/client/core.py`
- Create `_DeprecatedProgressTracker(vis, description)` class.
  - Implement `__enter__` to call `vis.api.progress_start(...)`.
  - Implement `__exit__` to call `vis.api.progress_complete(...)`.
  - Implement `update(description=None, progress=None)` which maps the absolute `progress` percentage to the new API format: `vis.api.progress_update(..., n=progress, total=100)`.
- Add the `vis.progress_tracker(description)` method back to the `ZnDraw` class.
- Decorate it with `@typing_extensions.deprecated("Use ZnDrawTqdm instead")`.
- Remove the mistakenly added `vis.progress_bar` shim (since `progress_bar` was never a real API in `v0.6.0`, it was a typo of `progress_tracker`).

### 2. Update Documentation
- **`docs/source/python-api.rst`**:
  - Remove references to `vis.progress_tracker()`.
  - Add standard usage examples for `ZnDrawTqdm` (e.g., `for item in ZnDrawTqdm(items, vis=vis, description="Processing"): ...`).
- **`README.md`**:
  - Verify if progress bars are mentioned and update to `ZnDrawTqdm` if necessary.

### 3. Write Tests
- Add a test in `tests/test_client_api.py` (or similar) to ensure `vis.progress_tracker()`:
  1. Emits a `DeprecationWarning`.
  2. Successfully connects and updates progress using the absolute percentage format without crashing.

### 4. Verification
- Run `uv run pytest tests/` to ensure no existing logic is broken.
- Ensure type-checking passes (`uv run pyright .`).
