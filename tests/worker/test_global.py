"""E2E tests for @global job scope.

Workers use ZnDraw.jobs exclusively for register, claim, submit, disconnect,
and status updates. Raw httpx is ONLY for read-only GET assertions.
"""

from typing import Any, ClassVar

import ase
import pytest
from zndraw_joblib.client import Category, ClaimedTask, Extension
from zndraw_joblib.schemas import JobResponse

from zndraw import ZnDraw

# =============================================================================
# Single Worker, Single Job (dev mode)
# =============================================================================


def test_register_global_job(server, Echo, get_job_list):
    """worker.jobs.register(Echo) makes Echo visible in GET job list."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


def test_registered_job_has_correct_schema(server, Echo):
    """The registered job's schema matches Echo.model_json_schema()."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        resp = worker.api.http.get(
            f"{worker.api.base_url}/v1/joblib/rooms/{worker.room}/jobs/@global:modifiers:Echo",
            headers=worker.api.get_headers(),
        )
        resp.raise_for_status()
        job = JobResponse.model_validate(resp.json())
        assert job.schema_ == Echo.model_json_schema()
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


def test_submit_and_auto_complete(server, Echo, wait_for_task):
    """Submitting a task auto-completes via the worker's claim loop."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        submitter.append(ase.Atoms("H"))
        worker.jobs.register(Echo)
        task_id = submitter.jobs.submit(
            Echo(value="test"), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed"
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_manual_control_no_auto_pickup(server, Echo, get_task):
    """With auto_pickup=False, tasks stay pending until manually claimed."""
    worker = ZnDraw(url=server, auto_pickup=False)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        task_id = submitter.jobs.submit(
            Echo(value="manual"), room=submitter.room, job_room="@global"
        )
        task = get_task(submitter, task_id)
        assert task.status == "pending"

        claimed = worker.jobs.claim()
        assert claimed is not None
        assert isinstance(claimed, ClaimedTask)
        worker.jobs.cancel(claimed)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_worker_auto_claims_task(server, Echo, wait_for_task):
    """Auto-serve claim loop picks up submitted tasks and runs them."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        submitter.append(ase.Atoms("H"))
        worker.jobs.register(Echo)
        task_id = submitter.jobs.submit(
            Echo(value="claim_test"), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed"
        assert submitter.bookmarks[0] == "claim_test"
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_full_lifecycle(server, Echo, run_worker_loop, wait_for_task):
    """Register -> submit -> worker loop -> verify bookmark in correct room."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        task_id = submitter.jobs.submit(
            Echo(value="proof"), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert submitter.bookmarks[0] == "proof"

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_claim_empty_queue_returns_none(server, Echo):
    """claim() returns None when no tasks are pending."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        claimed = worker.jobs.claim()
        assert claimed is None
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# Room Correctness Verification
# =============================================================================


def test_vis_connected_to_submitter_room(server, Echo, run_worker_loop, wait_for_task):
    """Worker writes bookmark to submitter's room, not another room."""
    worker = ZnDraw(url=server)
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        room_a.append(ase.Atoms("H"))
        room_b.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        # Submit from room_a only
        task_id = room_a.jobs.submit(
            Echo(value="room_a_proof"), room=room_a.room, job_room="@global"
        )
        result = wait_for_task(room_a, task_id)
        assert result.status == "completed"

        # Bookmark should exist in room_a
        assert room_a.bookmarks[0] == "room_a_proof"

        # Bookmark should NOT exist in room_b
        with pytest.raises(KeyError):
            _ = room_b.bookmarks[0]

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        room_a.disconnect()
        room_b.disconnect()


# =============================================================================
# Single Worker, Multiple Jobs
# =============================================================================


def test_register_two_jobs_same_worker(server, Echo, EchoAnalysis, get_job_list):
    """Registering Echo + EchoAnalysis on the same ZnDraw shows both in job list."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        worker.jobs.register(EchoAnalysis)

        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
        assert "EchoAnalysis" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


def test_worker_claims_from_multiple_jobs(
    server, Echo, EchoAnalysis, run_worker_loop, wait_for_task
):
    """Worker with two jobs claims tasks from both sequentially."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        worker.jobs.register(EchoAnalysis)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        tid1 = submitter.jobs.submit(
            Echo(value="echo_proof"), room=submitter.room, job_room="@global"
        )
        tid2 = submitter.jobs.submit(
            EchoAnalysis(value=42.0), room=submitter.room, job_room="@global"
        )

        r1 = wait_for_task(submitter, tid1)
        r2 = wait_for_task(submitter, tid2)
        assert r1.status == "completed", f"Echo failed: {r1.error}"
        assert r2.status == "completed", f"EchoAnalysis failed: {r2.error}"

        assert submitter.bookmarks[0] == "echo_proof"
        assert "echo_analysis" in submitter.figures

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


# =============================================================================
# Multiple Workers, Single Job
# =============================================================================


def test_two_workers_register_same_job(server, Echo):
    """Two workers registering the same job both appear as workers."""
    w1 = ZnDraw(url=server)
    w2 = ZnDraw(url=server)
    try:
        w1.jobs.register(Echo)
        w2.jobs.register(Echo)

        resp = w1.api.http.get(
            f"{w1.api.base_url}/v1/joblib/rooms/{w1.room}/jobs/@global:modifiers:Echo",
            headers=w1.api.get_headers(),
        )
        resp.raise_for_status()
        job = JobResponse.model_validate(resp.json())
        assert len(job.workers) == 2
    finally:
        w1.jobs.disconnect()
        w1.disconnect()
        w2.jobs.disconnect()
        w2.disconnect()


def test_two_workers_claim_different_tasks(
    server, Echo, run_worker_loop, wait_for_task
):
    """Two workers each claim a different task - no double-claim."""
    w1 = ZnDraw(url=server)
    w2 = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        w1.jobs.register(Echo)
        w2.jobs.register(Echo)
        submitter.append(ase.Atoms("H"))

        t1, stop1 = run_worker_loop(w1, server)
        t2, stop2 = run_worker_loop(w2, server)

        tid1 = submitter.jobs.submit(
            Echo(value="task1"), room=submitter.room, job_room="@global"
        )
        tid2 = submitter.jobs.submit(
            Echo(value="task2"), room=submitter.room, job_room="@global"
        )

        r1 = wait_for_task(submitter, tid1)
        r2 = wait_for_task(submitter, tid2)
        assert r1.status == "completed"
        assert r2.status == "completed"

        assert r1.worker_id is not None
        assert r2.worker_id is not None

        stop1.set()
        stop2.set()
        t1.join(timeout=5)
        t2.join(timeout=5)
    finally:
        w1.jobs.disconnect()
        w1.disconnect()
        w2.jobs.disconnect()
        w2.disconnect()
        submitter.disconnect()


# =============================================================================
# Multiple Workers, Multiple Jobs
# =============================================================================


def test_two_workers_two_jobs(
    server, Echo, EchoAnalysis, run_worker_loop, wait_for_task
):
    """Worker A: Echo, Worker B: EchoAnalysis -> each claims its own."""
    w_echo = ZnDraw(url=server)
    w_analysis = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        w_echo.jobs.register(Echo)
        w_analysis.jobs.register(EchoAnalysis)
        submitter.append(ase.Atoms("H"))

        t1, stop1 = run_worker_loop(w_echo, server)
        t2, stop2 = run_worker_loop(w_analysis, server)

        tid1 = submitter.jobs.submit(
            Echo(value="echo_val"), room=submitter.room, job_room="@global"
        )
        tid2 = submitter.jobs.submit(
            EchoAnalysis(value=99.0), room=submitter.room, job_room="@global"
        )

        r1 = wait_for_task(submitter, tid1)
        r2 = wait_for_task(submitter, tid2)
        assert r1.status == "completed", f"Echo failed: {r1.error}"
        assert r2.status == "completed", f"EchoAnalysis failed: {r2.error}"

        assert submitter.bookmarks[0] == "echo_val"
        assert "echo_analysis" in submitter.figures

        stop1.set()
        stop2.set()
        t1.join(timeout=5)
        t2.join(timeout=5)
    finally:
        w_echo.jobs.disconnect()
        w_echo.disconnect()
        w_analysis.jobs.disconnect()
        w_analysis.disconnect()
        submitter.disconnect()


# =============================================================================
# Schema & Validation
# =============================================================================


def test_schema_conflict_409(server, Echo):
    """Registering a job with same name but different schema raises ValueError (409)."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)

        # Local class with same __name__ but different fields
        class Echo(Extension):  # type: ignore[no-redef]
            category: ClassVar[Category] = Category.MODIFIER
            different_field: int = 0

            def run(self, vis: Any, **kwargs: Any) -> Any:
                pass

        with pytest.raises(ValueError, match="Schema mismatch"):
            worker.jobs.register(Echo)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


def test_failing_extension_auto_marked_failed(server, FailingExtension, wait_for_task):
    """A failing extension is auto-marked as 'failed' by the claim loop."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        submitter.append(ase.Atoms("H"))
        worker.jobs.register(FailingExtension)
        task_id = submitter.jobs.submit(
            FailingExtension(), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "failed"
        assert result.error is not None
        assert "intentional failure" in result.error
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_invalid_state_transition_409(server, Echo, run_worker_loop, wait_for_task):
    """Transitioning a completed task to running returns 409."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        task_id = submitter.jobs.submit(
            Echo(value="x"), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed"

        # Try invalid transition: completed -> running via raw PATCH
        resp = submitter.api.http.patch(
            f"{submitter.api.base_url}/v1/joblib/tasks/{task_id}",
            headers=submitter.api.get_headers(),
            json={"status": "running"},
        )
        assert resp.status_code == 409

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


# =============================================================================
# Error Handling
# =============================================================================


def test_fail_records_error_message(
    server, FailingExtension, run_worker_loop, wait_for_task
):
    """manager.fail() records the error message on the task."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(FailingExtension)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        task_id = submitter.jobs.submit(
            FailingExtension(), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "failed"
        assert result.error is not None
        assert "intentional failure" in result.error

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


# =============================================================================
# Worker Cleanup & Disconnect
# =============================================================================


def test_worker_disconnect_clears_state(server, Echo):
    """Disconnecting a worker clears its worker_id and registry."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        assert worker.jobs.worker_id is not None
        assert len(worker.jobs) > 0

        worker.jobs.disconnect()

        assert worker.jobs.worker_id is None
        assert len(worker.jobs) == 0
    finally:
        worker.disconnect()


def test_worker_disconnect_fails_claimed_tasks(server, Echo, get_task):
    """Disconnecting a worker fails its claimed (non-terminal) tasks."""
    worker = ZnDraw(url=server, auto_pickup=False)
    submitter = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo)
        task_id = submitter.jobs.submit(
            Echo(value="x"), room=submitter.room, job_room="@global"
        )

        # Claim but don't complete
        claimed = worker.jobs.claim()
        assert claimed is not None

        # Disconnect -> task should be failed
        worker.jobs.disconnect()

        task = get_task(submitter, task_id)
        assert task.status == "failed"
    finally:
        worker.disconnect()
        submitter.disconnect()


def test_auto_serve_completes_task_on_submit(server, wait_for_task):
    """Auto-serve picks up and completes a task without manual claim()."""

    # Use a unique extension so no other test's workers are linked to this job
    class AutoOnly(Extension):
        category: ClassVar[Category] = Category.MODIFIER

        def run(self, vis: Any, **kwargs: Any) -> None:
            vis.bookmarks[0] = "auto"

    AutoOnly.__name__ = "AutoOnly"

    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server)
    try:
        submitter.append(ase.Atoms("H"))
        worker.jobs.register(AutoOnly)
        task_id = submitter.jobs.submit(
            AutoOnly(), room=submitter.room, job_room="@global"
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed"
        assert submitter.bookmarks[0] == "auto"
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


def test_worker_connect_disconnect_connect(server, Echo, get_job_list):
    """A worker can register, disconnect, then a new worker registers same job."""
    w1 = ZnDraw(url=server)
    try:
        w1.jobs.register(Echo)
        w1.jobs.disconnect()
    finally:
        w1.disconnect()

    w2 = ZnDraw(url=server)
    try:
        w2.jobs.register(Echo)
        jobs = get_job_list(w2, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        w2.jobs.disconnect()
        w2.disconnect()


# =============================================================================
# Auth Mode (server_auth)
# =============================================================================


def test_guest_cannot_register_global_job(server_auth, Echo):
    """Guest user cannot register @global jobs (requires admin)."""
    guest = ZnDraw(url=server_auth)
    try:
        with pytest.raises(PermissionError):
            guest.jobs.register(Echo)
    finally:
        guest.disconnect()


def test_admin_can_register_global_job(server_auth, Echo, get_job_list):
    """Admin user can register @global jobs."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    try:
        admin.jobs.register(Echo)
        jobs = get_job_list(admin, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        admin.jobs.disconnect()
        admin.disconnect()


def test_guest_can_submit_to_global_job(server_auth, Echo, wait_for_task):
    """Guest can submit tasks to @global jobs registered by admin."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    try:
        admin.jobs.register(Echo)
    except Exception:
        admin.disconnect()
        raise

    guest = ZnDraw(url=server_auth)
    try:
        guest.append(ase.Atoms("H"))
        task_id = guest.jobs.submit(
            Echo(value="guest_submit"), room=guest.room, job_room="@global"
        )
        assert isinstance(task_id, str)
        result = wait_for_task(guest, task_id)
        assert result.status == "completed"
        assert guest.bookmarks[0] == "guest_submit"
    finally:
        admin.jobs.disconnect()
        admin.disconnect()
        guest.disconnect()


def test_admin_full_lifecycle(server_auth, Echo, run_worker_loop, wait_for_task):
    """Admin registers, guest submits, admin worker loop completes with bookmark."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    try:
        admin.jobs.register(Echo)
    except Exception:
        admin.disconnect()
        raise

    guest = ZnDraw(url=server_auth)
    try:
        guest.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(admin, server_auth)

        task_id = guest.jobs.submit(
            Echo(value="admin_proof"), room=guest.room, job_room="@global"
        )
        result = wait_for_task(guest, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert guest.bookmarks[0] == "admin_proof"

        stop.set()
        thread.join(timeout=5)
    finally:
        admin.jobs.disconnect()
        admin.disconnect()
        guest.disconnect()
