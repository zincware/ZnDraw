"""E2E tests for @internal (server-side) extension execution.

@internal extensions are registered at startup and executed by the in-process
TaskIQ worker. The full path is:
  vis.run(extension) -> submit_task -> TaskIQ -> execute_extension
  -> ZnDraw(url, room_id) -> extension.run(vis)
"""

import ase
import pytest
from zndraw_joblib.schemas import JobResponse

from zndraw import ZnDraw
from zndraw.extensions.modifiers import Delete
from zndraw.extensions.selections import All, NoneSelection

# =============================================================================
# Registration (read-only verification)
# =============================================================================


def test_internal_jobs_listed_on_startup(server, get_job_list):
    """Server registers built-in extensions as @internal jobs at startup."""
    vis = ZnDraw(url=server)
    try:
        jobs = get_job_list(vis, room_id="@internal")
        names = {j.name for j in jobs}
        assert "NoneSelection" in names
        assert "All" in names
        assert "Delete" in names
    finally:
        vis.disconnect()


def test_internal_job_has_valid_schema(server):
    """An @internal job's schema contains a title or properties key."""
    vis = ZnDraw(url=server)
    try:
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/jobs/@internal:selections:All",
            headers=vis.api.get_headers(),
        )
        resp.raise_for_status()
        job = JobResponse.model_validate(resp.json())
        assert "properties" in job.schema_ or "title" in job.schema_
    finally:
        vis.disconnect()


# =============================================================================
# Submission
# =============================================================================


def test_submit_returns_task_id(server):
    """vis.run() for an @internal extension returns a non-empty task ID."""
    vis = ZnDraw(url=server)
    try:
        vis.append(ase.Atoms("H"))
        task_id = vis.run(All())
        assert isinstance(task_id, str)
        assert len(task_id) > 0
    finally:
        vis.disconnect()


def test_submit_nonexistent_internal_errors(server, Echo):
    """Submitting a task for a non-existent @internal job raises an error."""
    vis = ZnDraw(url=server)
    try:
        # Echo is NOT registered as @internal (only built-in extensions are)
        with pytest.raises(KeyError, match="not found"):
            vis.jobs.submit(Echo(), room=vis.room, job_room="@internal")
    finally:
        vis.disconnect()


# =============================================================================
# Full E2E (side-effect verification)
# =============================================================================


def test_selection_all_e2e(server, wait_for_task):
    """All() selects every atom in the current frame."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("HHHHH", positions=[[i, 0, 0] for i in range(5)])
        vis.append(atoms)

        task_id = vis.run(All())
        result = wait_for_task(vis, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert vis.selection == tuple(range(5))
    finally:
        vis.disconnect()


def test_selection_none_e2e(server, wait_for_task):
    """NoneSelection() clears the selection."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
        vis.append(atoms)
        vis.selection = [0, 1]

        task_id = vis.run(NoneSelection())
        result = wait_for_task(vis, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert vis.selection == ()
    finally:
        vis.disconnect()


def test_modifier_delete_e2e(server, wait_for_task):
    """Delete() removes selected atoms and appends a new frame."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("HHHHH", positions=[[i, 0, 0] for i in range(5)])
        vis.append(atoms)
        vis.selection = [1, 3]

        task_id = vis.run(Delete())
        result = wait_for_task(vis, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        # Delete appends a new frame with selected atoms removed
        assert len(vis) >= 2
        new_frame = vis[-1]
        assert len(new_frame) == 3
    finally:
        vis.disconnect()


def test_internal_task_with_invalid_payload_fails(server, wait_for_task):
    """Submitting an @internal extension with invalid payload reaches 'failed' status."""
    vis = ZnDraw(url=server)
    try:
        vis.append(ase.Atoms("H"))

        # Random requires 'count' (no default), so empty payload is invalid.
        # Use raw HTTP to bypass client-side model validation.
        resp = vis.api.http.post(
            f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/tasks/@internal:selections:Random",
            headers=vis.api.get_headers(),
            json={"payload": {}},
        )
        assert resp.status_code == 202
        task_id = resp.json()["id"]

        result = wait_for_task(vis, task_id)
        assert result.status == "failed"
        assert result.error is not None
    finally:
        vis.disconnect()


def test_sequential_internal_tasks(server, wait_for_task):
    """Multiple sequential @internal tasks all complete."""
    vis = ZnDraw(url=server)
    try:
        atoms = ase.Atoms("HH", positions=[[0, 0, 0], [1, 0, 0]])
        vis.append(atoms)

        for i in range(3):
            task_id = vis.run(NoneSelection())
            result = wait_for_task(vis, task_id)
            assert result.status == "completed", f"Task {i} failed: {result.error}"
    finally:
        vis.disconnect()
