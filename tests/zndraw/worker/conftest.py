"""Shared fixtures and test extensions for worker E2E tests."""

import threading
from collections.abc import Callable
from pathlib import Path
from typing import Any, ClassVar

import ase
import ase.io
import pytest

from zndraw import ZnDraw
from zndraw_joblib.client import Category, Extension
from zndraw_joblib.schemas import JobSummary, TaskResponse

_WATER = ase.Atoms("H2O", positions=[[0, 0, 0], [0, 0, 1], [1, 0, 0]])
_ASEBYTES_EXTENSIONS = frozenset({".h5", ".h5md", ".lmdb"})


def _write_water(path: Path, frames: list[ase.Atoms]) -> None:
    """Write *frames* to *path*, routing asebytes formats through ``ASEIO``."""
    if path.suffix in _ASEBYTES_EXTENSIONS:
        import asebytes

        with asebytes.ASEIO(str(path)) as db:
            db.extend(frames)
    else:
        ase.io.write(path, frames)


@pytest.fixture
def water_xyz(tmp_path: Path) -> Path:
    """Write a 3-atom H2O xyz file to ``tmp_path / water.xyz`` and return its path.

    Shared by any test that needs a real on-disk structure file (e.g. the
    ``@internal:modifiers:LoadFile`` e2e path).
    """
    path = tmp_path / "water.xyz"
    _write_water(path, [_WATER])
    return path


@pytest.fixture(params=[".xyz", ".h5", ".h5md"])
def water_file(request: pytest.FixtureRequest, tmp_path: Path) -> Path:
    """Parametric fixture: one H2O frame written in every format LoadFile must handle.

    Covers both branches of ``zndraw.io.open_frames`` — ``ase.io.iread`` for
    streaming text formats and ``asebytes.ASEIO`` for random-access formats.
    ``.lmdb`` is omitted: the asebytes LMDB backend returns mmap-backed
    read-only arrays which downstream ``atoms.set_pbc(...)`` can't mutate —
    a pre-existing issue that also breaks the ``zndraw <file>.lmdb`` CLI
    path and is unrelated to #923.
    """
    suffix = request.param
    path = tmp_path / f"water{suffix}"
    _write_water(path, [_WATER])
    return path


@pytest.fixture
def water_trajectory_h5(tmp_path: Path) -> Path:
    """Five-frame H2O trajectory in h5 format, for slice tests against LoadFile."""
    path = tmp_path / "trj.h5"
    _write_water(path, [_WATER for _ in range(5)])
    return path


# =============================================================================
# Test Extension Classes
# =============================================================================


class _Echo(Extension):
    """Modifier that writes a bookmark to prove room correctness."""

    category: ClassVar[Category] = Category.MODIFIER
    value: str = "hello"

    def run(self, vis: Any, **_kwargs: Any) -> None:
        vis.bookmarks[0] = self.value


_Echo.__name__ = "Echo"


class _EchoAnalysis(Extension):
    """Analysis that writes a figure to prove room correctness."""

    category: ClassVar[Category] = Category.ANALYSIS
    value: float = 0.0

    def run(self, vis: Any, **_kwargs: Any) -> None:
        import plotly.graph_objects as go

        fig = go.Figure(layout={"title": str(self.value)})
        vis.figures["echo_analysis"] = fig


_EchoAnalysis.__name__ = "EchoAnalysis"


class _FailingExtension(Extension):
    """Extension that always raises an error."""

    category: ClassVar[Category] = Category.MODIFIER

    def run(self, _vis: Any, **_kwargs: Any) -> None:
        raise RuntimeError("intentional failure")


_FailingExtension.__name__ = "FailingExtension"


# =============================================================================
# Extension class fixtures
# =============================================================================


@pytest.fixture
def Echo() -> type[_Echo]:  # noqa: N802
    """The Echo test extension class."""
    return _Echo


@pytest.fixture
def EchoAnalysis() -> type[_EchoAnalysis]:  # noqa: N802
    """The EchoAnalysis test extension class."""
    return _EchoAnalysis


@pytest.fixture
def FailingExtension() -> type[_FailingExtension]:  # noqa: N802
    """The FailingExtension test extension class."""
    return _FailingExtension


# =============================================================================
# Helper fixtures (return callables)
# =============================================================================


@pytest.fixture
def wait_for_task() -> Callable[..., TaskResponse]:
    """Long-poll GET until task reaches terminal status."""

    def _wait(vis: ZnDraw, task_id: str, *, timeout: int = 30) -> TaskResponse:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/tasks/{task_id}",
            headers={**vis.api.get_headers(), "Prefer": f"wait={timeout}"},
            timeout=timeout + 5,
        )
        resp.raise_for_status()
        return TaskResponse.model_validate(resp.json())

    return _wait


@pytest.fixture
def get_job_list() -> Callable[..., list[JobSummary]]:
    """GET job list for a room (read-only assertion helper)."""

    def _get(vis: ZnDraw, room_id: str = "@global") -> list[JobSummary]:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/rooms/{room_id}/jobs",
            headers=vis.api.get_headers(),
        )
        resp.raise_for_status()
        return [JobSummary.model_validate(j) for j in resp.json()["items"]]

    return _get


@pytest.fixture
def get_task() -> Callable[..., TaskResponse]:
    """GET a single task by ID (read-only assertion helper)."""

    def _get(vis: ZnDraw, task_id: str) -> TaskResponse:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/tasks/{task_id}",
            headers=vis.api.get_headers(),
        )
        resp.raise_for_status()
        return TaskResponse.model_validate(resp.json())

    return _get


@pytest.fixture
def run_worker_loop() -> Callable[..., tuple[threading.Thread, threading.Event]]:
    """Start a real worker listen() loop in a background thread.

    Returns (thread, stop_event). Call stop.set(); thread.join(timeout=5)
    to shut down.
    """

    def _start(worker: ZnDraw) -> tuple[threading.Thread, threading.Event]:
        stop = threading.Event()

        def _loop() -> None:
            try:
                for task in worker.jobs.listen(polling_interval=0.5, stop_event=stop):
                    worker.jobs.start(task)
                    try:
                        worker._execute_task(task)
                    except Exception as e:  # noqa: BLE001
                        worker.jobs.fail(task, str(e))
                    else:
                        worker.jobs.complete(task)
            except Exception:  # noqa: BLE001, S110
                pass  # Worker disconnected — stop gracefully

        thread = threading.Thread(target=_loop, daemon=True)
        thread.start()
        return thread, stop

    return _start
