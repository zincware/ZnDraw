# Global Extension Bugs — Design Spec

**Date:** 2026-04-01
**Scope:** Four bugs found during global extension investigation

## Bug 1: `list_extensions()` defaults to `@internal`, hiding global extensions

### Problem

`APIManager.list_extensions()` in `src/zndraw/client/api.py:854` uses
`room or "@internal"` as fallback. Since callers (`Extensions` accessor,
CLI `extensions list`) never pass `room`, the method always queries
`GET /v1/joblib/rooms/@internal/jobs` — which only returns `@internal`
extensions, never `@global` or room-scoped ones.

The server-side `_room_job_filter(room_id)` in `src/zndraw_joblib/router.py:230`
correctly returns `@global + @internal + room-scoped` jobs when queried with a
real room ID. The frontend already does this correctly.

### Symptoms

- `uv run zndraw-cli extensions list` never shows `@global` extensions
- `list(vis.extensions)` only returns `@internal` extensions
- `"@global:modifiers:X" in vis.extensions` returns `False` (Mapping protocol violation)
- Existing tests pass because they bypass `list_extensions()` and query the HTTP API directly

### Fix

Change `src/zndraw/client/api.py:854`:

```python
# Before
job_room = room or "@internal"

# After
job_room = room or self.room_id
```

No changes needed in `accessors.py` or `cli_agent/extensions.py` — they call
`list_extensions()` without `room`, which now correctly falls back to `self.room_id`.

### Test

Add test in `tests/zndraw/worker/` that registers a global extension, then verifies
`vis.api.list_extensions()` (no room arg) includes it — the exact code path the
existing tests bypass.

---

## Bug 2: No global exception handler — bare 500s with no diagnostics

### Problem

`src/zndraw/app.py` registers exception handlers for `ProblemError`,
`JoblibProblemError`, and `RequestValidationError` only. Any unhandled exception
falls through to Starlette's default handler, which returns plain text
`Internal Server Error` with no RFC 9457 body, no traceback in the response,
and potentially no server-side logging.

This made the remote 500 error (Bug 3) impossible to diagnose.

### Fix

Add `InternalServerError(ProblemType)` to `src/zndraw/exceptions.py`:

```python
class InternalServerError(ProblemType):
    """An unexpected error occurred while processing the request.

    This error is returned when an unhandled exception reaches the
    top-level exception handler. The ``detail`` field contains the
    exception message; the full traceback is logged server-side.
    """

    title: ClassVar[str] = "Internal Server Error"
    status: ClassVar[int] = 500
```

Register in `PROBLEM_TYPES`.

Add `import logging` and `logger = logging.getLogger(__name__)` to `src/zndraw/app.py`.

Add catch-all in `src/zndraw/app.py`:

```python
@app.exception_handler(Exception)
async def _unhandled_exception_handler(
    request: Request, exc: Exception
) -> JSONResponse:
    """Catch-all for unhandled exceptions — log and return RFC 9457."""
    logger.error(
        "Unhandled %s on %s %s",
        type(exc).__name__, request.method, request.url.path,
        exc_info=True,
    )
    return await problem_exception_handler(
        request, InternalServerError.exception(detail=str(exc))
    )
```

Logs full traceback server-side, delegates to the existing canonical
`problem_exception_handler` for RFC 9457 response.

Refactor `_validation_exception_handler` to delegate instead of building its
own `JSONResponse`:

```python
@app.exception_handler(RequestValidationError)
async def _validation_exception_handler(
    request: Request, exc: RequestValidationError
) -> JSONResponse:
    """Convert FastAPI validation errors to RFC 9457 problem detail."""
    detail = "; ".join(
        f"{'.'.join(str(x) for x in e['loc'])}: {e['msg']}" for e in exc.errors()
    )
    return await problem_exception_handler(
        request, UnprocessableContent.exception(detail=detail)
    )
```

### Test

Assert that an unhandled exception returns `application/problem+json` with
`status: 500` and `type: "/v1/problems/internal-server-error"`.

---

## Bug 3: Regression test for nested dict+numpy frames

### Problem

Appending frames with nested dicts containing numpy arrays in `atoms.info`
(e.g. `{"grid": np.array(...), "origin": np.array(...), "cell": np.array(...)}`)
causes a 500 on the remote server. Works locally with the same version (0.7.0a7).
Root cause is unknown — Bug 2's fix is needed to get diagnostics.

The `Isosurface` geometry feature depends on this data shape.

### Fix

Add regression test in `tests/zndraw/test_routes_frames.py` that:

1. Creates `ase.Atoms` with `atoms.info["cube_data"] = {"grid": np.randn(10,10,10), ...}`
2. Serializes via `atoms_to_json_dict`
3. POSTs to the frames endpoint
4. Asserts 201
5. GETs it back and verifies the nested dict round-trips

This ensures the code path works locally and catches future regressions.

---

## Bug 4: Misleading `JobManager.register()` docstring

### Problem

`src/zndraw_joblib/client.py:375` says `room` "Defaults to `@global`". This is
technically true for `JobManager.register()` in isolation, but
`ZnDraw.register_job()` resolves `None` to `self.room` before calling it. The
effective default for end users is room-scoped, not global.

### Fix

Update the docstring to clarify:

```
room
    Room scope. Defaults to ``"@global"`` when called directly.
    ``ZnDraw.register_job()`` resolves *None* to the client's room.
```

---

## Summary

| Bug | Files | Change |
|-----|-------|--------|
| 1 | `src/zndraw/client/api.py` | `room or self.room_id` |
| 2 | `src/zndraw/exceptions.py`, `src/zndraw/app.py` | `InternalServerError` + catch-all handler |
| 3 | `tests/zndraw/test_routes_frames.py` | Regression test for nested dict+numpy |
| 4 | `src/zndraw_joblib/client.py` | Docstring clarification |

Plus tests for Bug 1 (list_extensions includes @global) and Bug 2 (500 returns RFC 9457).
