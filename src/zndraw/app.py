from pathlib import Path

import socketio
from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from starlette.staticfiles import StaticFiles
from zndraw_joblib import (
    ProblemException as JoblibProblemException,
    problem_exception_handler as joblib_problem_exception_handler,
    router as joblib_router,
)
from zndraw_joblib.dependencies import (
    verify_writable_room as joblib_verify_writable_room,
)

from zndraw.database import lifespan
from zndraw.dependencies import get_writable_room_id
from zndraw.exceptions import (
    ProblemException,
    UnprocessableContent,
    problem_exception_handler,
)
from zndraw.routes.admin import router as admin_router
from zndraw.routes.auth import router as auth_router
from zndraw.routes.bookmarks import router as bookmarks_router
from zndraw.routes.chat import router as chat_router
from zndraw.routes.edit_lock import router as edit_lock_router
from zndraw.routes.figures import router as figures_router
from zndraw.routes.frames import router as frames_router
from zndraw.routes.geometries import default_camera_router, router as geometries_router
from zndraw.routes.problems import router as problems_router
from zndraw.routes.progress import router as progress_router
from zndraw.routes.rooms import router as rooms_router
from zndraw.routes.screenshots import router as screenshots_router
from zndraw.routes.selection_groups import router as selection_groups_router
from zndraw.routes.server_settings import router as server_settings_router
from zndraw.routes.step import router as step_router
from zndraw.routes.tools import router as tools_router
from zndraw.routes.trajectory import router as trajectory_router
from zndraw.routes.utility import router as utility_router
from zndraw.socketio import tsio

app = FastAPI(title="ZnDraw API", lifespan=lifespan)

# Override joblib's verify_writable_room to enforce room locks
app.dependency_overrides[joblib_verify_writable_room] = get_writable_room_id

# Register exception handlers
app.add_exception_handler(ProblemException, problem_exception_handler)
app.add_exception_handler(JoblibProblemException, joblib_problem_exception_handler)


@app.exception_handler(RequestValidationError)
async def _validation_exception_handler(
    _request: Request, exc: RequestValidationError
) -> JSONResponse:
    """Convert FastAPI validation errors to RFC 9457 problem detail."""
    detail = "; ".join(
        f"{'.'.join(str(x) for x in e['loc'])}: {e['msg']}" for e in exc.errors()
    )
    problem = UnprocessableContent.create(detail=detail)
    return JSONResponse(
        status_code=422,
        content=problem.model_dump(exclude_none=True),
        media_type="application/problem+json",
    )


# Include routers
app.include_router(admin_router)
app.include_router(auth_router)
app.include_router(bookmarks_router)
app.include_router(chat_router)
app.include_router(figures_router)
app.include_router(frames_router)
app.include_router(geometries_router)
app.include_router(default_camera_router)
app.include_router(edit_lock_router)
app.include_router(problems_router)
app.include_router(progress_router)
app.include_router(server_settings_router)
app.include_router(rooms_router)
app.include_router(screenshots_router)
app.include_router(selection_groups_router)
app.include_router(step_router)
app.include_router(tools_router)
app.include_router(trajectory_router)
app.include_router(utility_router)
app.include_router(joblib_router)

# Serve built frontend assets and SPA catch-all.
# Only mount if static/ exists — during development, Vite dev server proxies instead.
_static_dir = Path(__file__).parent / "static"
if _static_dir.is_dir():
    # Serve hashed JS/CSS bundles at /assets/
    app.mount("/assets", StaticFiles(directory=_static_dir / "assets"))

    # SPA catch-all: any path not matched by API routers serves index.html
    @app.get("/{path:path}", include_in_schema=False)
    async def _spa_catch_all(path: str) -> FileResponse:
        return FileResponse(_static_dir / "index.html")


# Mount Socket.IO
socket_app = socketio.ASGIApp(tsio, other_asgi_app=app)
