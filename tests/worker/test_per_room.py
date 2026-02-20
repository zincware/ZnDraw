"""E2E tests for per-room scoped jobs.

Per-room jobs are registered with room=vis.room and are only visible
to clients in that room. Cross-room isolation is a key property.
"""

import uuid

import ase

from zndraw import ZnDraw

# =============================================================================
# Registration & Isolation
# =============================================================================


def test_register_room_scoped_job(server, Echo, get_job_list):
    """Can register a job scoped to a specific room."""
    vis = ZnDraw(url=server)
    try:
        vis.jobs.register(Echo, room=vis.room)
        jobs = get_job_list(vis, room_id=vis.room)
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_room_job_not_visible_in_other_room(server, Echo, get_job_list):
    """A room-scoped job is not visible from a different room."""
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        room_a.jobs.register(Echo, room=room_a.room)

        jobs = get_job_list(room_b, room_id=room_b.room)
        room_scoped = [j for j in jobs if j.full_name.startswith(room_a.room)]
        assert len(room_scoped) == 0
    finally:
        room_a.jobs.disconnect()
        room_a.disconnect()
        room_b.disconnect()


def test_room_job_visible_in_own_room(server, Echo, get_job_list):
    """A room-scoped job is visible from its own room."""
    vis = ZnDraw(url=server)
    try:
        vis.jobs.register(Echo, room=vis.room)
        jobs = get_job_list(vis, room_id=vis.room)
        room_scoped = [j for j in jobs if j.full_name.startswith(vis.room)]
        assert len(room_scoped) == 1
        assert room_scoped[0].name == "Echo"
    finally:
        vis.jobs.disconnect()
        vis.disconnect()


def test_global_job_visible_from_any_room(server, Echo, get_job_list):
    """@global jobs are visible from any room."""
    registrar = ZnDraw(url=server)
    viewer = ZnDraw(url=server)
    try:
        registrar.jobs.register(Echo)  # defaults to @global
        jobs = get_job_list(viewer, room_id=viewer.room)
        global_jobs = [j for j in jobs if j.full_name.startswith("@global")]
        names = {j.name for j in global_jobs}
        assert "Echo" in names
    finally:
        registrar.jobs.disconnect()
        registrar.disconnect()
        viewer.disconnect()


# =============================================================================
# Task Submission & Claiming (single worker)
# =============================================================================


def test_submit_to_room_job_auto_completes(server, Echo, wait_for_task):
    """Submitting to a room-scoped job auto-completes via claim loop."""
    worker = ZnDraw(url=server)
    try:
        worker.append(ase.Atoms("H"))
        worker.jobs.register(Echo, room=worker.room)
        task_id = worker.jobs.submit(
            Echo(value="room_test"),
            room=worker.room,
            job_room=worker.room,
        )
        result = wait_for_task(worker, task_id)
        assert result.status == "completed"
        assert result.room_id == worker.room
        assert worker.bookmarks[0] == "room_test"
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


def test_full_room_lifecycle(server, Echo, run_worker_loop, wait_for_task):
    """Register in room -> submit -> worker loop -> verify bookmark."""
    worker = ZnDraw(url=server)
    submitter = ZnDraw(url=server, room=worker.room)
    try:
        worker.jobs.register(Echo, room=worker.room)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server)

        task_id = submitter.jobs.submit(
            Echo(value="room_proof"),
            room=submitter.room,
            job_room=worker.room,
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert submitter.bookmarks[0] == "room_proof"

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()


# =============================================================================
# Task Submission & Claiming (two workers)
# =============================================================================


def test_two_workers_same_room_job(server, Echo, run_worker_loop, wait_for_task):
    """Two workers in the same room claim different tasks via listen() loops."""
    room_id_shared = uuid.uuid4().hex
    w1 = ZnDraw(url=server, room=room_id_shared)
    w2 = ZnDraw(url=server, room=room_id_shared)
    submitter = ZnDraw(url=server, room=room_id_shared)
    try:
        w1.jobs.register(Echo, room=room_id_shared)
        w2.jobs.register(Echo, room=room_id_shared)
        submitter.append(ase.Atoms("H"))

        t1, stop1 = run_worker_loop(w1, server)
        t2, stop2 = run_worker_loop(w2, server)

        tid1 = submitter.jobs.submit(
            Echo(value="t1"), room=room_id_shared, job_room=room_id_shared
        )
        tid2 = submitter.jobs.submit(
            Echo(value="t2"), room=room_id_shared, job_room=room_id_shared
        )

        r1 = wait_for_task(submitter, tid1)
        r2 = wait_for_task(submitter, tid2)
        assert r1.status == "completed"
        assert r2.status == "completed"

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


def test_two_workers_submit_and_complete(server, Echo, run_worker_loop, wait_for_task):
    """Both workers run via listen() loops, sequential tasks verified."""
    room_id_shared = uuid.uuid4().hex
    w1 = ZnDraw(url=server, room=room_id_shared)
    w2 = ZnDraw(url=server, room=room_id_shared)
    submitter = ZnDraw(url=server, room=room_id_shared)
    try:
        w1.jobs.register(Echo, room=room_id_shared)
        w2.jobs.register(Echo, room=room_id_shared)
        submitter.append(ase.Atoms("H"))

        t1, stop1 = run_worker_loop(w1, server)
        t2, stop2 = run_worker_loop(w2, server)

        tid1 = submitter.jobs.submit(
            Echo(value="first"), room=room_id_shared, job_room=room_id_shared
        )
        r1 = wait_for_task(submitter, tid1)
        assert r1.status == "completed"
        assert submitter.bookmarks[0] == "first"

        tid2 = submitter.jobs.submit(
            Echo(value="second"), room=room_id_shared, job_room=room_id_shared
        )
        r2 = wait_for_task(submitter, tid2)
        assert r2.status == "completed"
        assert submitter.bookmarks[0] == "second"

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
# Worker Connect/Disconnect
# =============================================================================


def test_worker_disconnect_in_room_clears_state(server, Echo):
    """Disconnecting a room-scoped worker clears its worker_id and registry."""
    worker = ZnDraw(url=server)
    try:
        worker.jobs.register(Echo, room=worker.room)
        assert worker.jobs.worker_id is not None
        assert len(worker.jobs) > 0

        worker.jobs.disconnect()

        assert worker.jobs.worker_id is None
        assert len(worker.jobs) == 0
    finally:
        worker.disconnect()


def test_worker_disconnect_fails_claimed_tasks_in_room(server, Echo, get_task):
    """Disconnecting a room-scoped worker fails its claimed tasks."""
    worker = ZnDraw(url=server, auto_pickup=False)
    checker = ZnDraw(url=server, room=worker.room)
    try:
        worker.jobs.register(Echo, room=worker.room)
        task_id = worker.jobs.submit(
            Echo(value="x"), room=worker.room, job_room=worker.room
        )

        claimed = worker.jobs.claim()
        assert claimed is not None

        worker.jobs.disconnect()

        task = get_task(checker, task_id)
        assert task.status == "failed"
    finally:
        worker.disconnect()
        checker.disconnect()


def test_worker_disconnect_reconnect_room(server, Echo, get_job_list):
    """Worker registers in room, disconnects, new worker registers same job."""
    room_id = ZnDraw(url=server).room

    w1 = ZnDraw(url=server, room=room_id)
    try:
        w1.jobs.register(Echo, room=room_id)
        w1.jobs.disconnect()
    finally:
        w1.disconnect()

    w2 = ZnDraw(url=server, room=room_id)
    try:
        w2.jobs.register(Echo, room=room_id)
        jobs = get_job_list(w2, room_id=room_id)
        room_jobs = [j for j in jobs if j.full_name.startswith(room_id)]
        names = {j.name for j in room_jobs}
        assert "Echo" in names
    finally:
        w2.jobs.disconnect()
        w2.disconnect()


# =============================================================================
# Cross-Room Isolation
# =============================================================================


def test_cross_room_claim_isolation(server, Echo):
    """Worker in room_A cannot claim tasks submitted to room_B's job."""
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        room_a.jobs.register(Echo, room=room_a.room)
        room_b.jobs.register(Echo, room=room_b.room)

        # Submit to room_B's job
        room_b.jobs.submit(
            Echo(value="b_task"),
            room=room_b.room,
            job_room=room_b.room,
        )

        # room_A worker should NOT be able to claim room_B's task
        claimed = room_a.jobs.claim()
        assert claimed is None
    finally:
        room_a.jobs.disconnect()
        room_a.disconnect()
        room_b.jobs.disconnect()
        room_b.disconnect()


def test_two_rooms_two_workers(server, Echo, run_worker_loop, wait_for_task):
    """Each room's worker only processes its own tasks, correct side-effects."""
    w_a = ZnDraw(url=server)
    w_b = ZnDraw(url=server)
    sub_a = ZnDraw(url=server, room=w_a.room)
    sub_b = ZnDraw(url=server, room=w_b.room)
    try:
        w_a.jobs.register(Echo, room=w_a.room)
        w_b.jobs.register(Echo, room=w_b.room)
        sub_a.append(ase.Atoms("H"))
        sub_b.append(ase.Atoms("H"))

        ta, stop_a = run_worker_loop(w_a, server)
        tb, stop_b = run_worker_loop(w_b, server)

        tid_a = sub_a.jobs.submit(
            Echo(value="proof_a"), room=w_a.room, job_room=w_a.room
        )
        tid_b = sub_b.jobs.submit(
            Echo(value="proof_b"), room=w_b.room, job_room=w_b.room
        )

        ra = wait_for_task(sub_a, tid_a)
        rb = wait_for_task(sub_b, tid_b)
        assert ra.status == "completed"
        assert rb.status == "completed"

        assert sub_a.bookmarks[0] == "proof_a"
        assert sub_b.bookmarks[0] == "proof_b"

        stop_a.set()
        stop_b.set()
        ta.join(timeout=5)
        tb.join(timeout=5)
    finally:
        w_a.jobs.disconnect()
        w_a.disconnect()
        w_b.jobs.disconnect()
        w_b.disconnect()
        sub_a.disconnect()
        sub_b.disconnect()


# =============================================================================
# Auth Mode — Room-scoped
# =============================================================================


def test_guest_can_register_room_job(server_auth, Echo, get_job_list):
    """Guest user can register a room-scoped job (no admin required)."""
    guest = ZnDraw(url=server_auth)
    try:
        guest.jobs.register(Echo, room=guest.room)
        jobs = get_job_list(guest, room_id=guest.room)
        room_jobs = [j for j in jobs if j.full_name.startswith(guest.room)]
        names = {j.name for j in room_jobs}
        assert "Echo" in names
    finally:
        guest.jobs.disconnect()
        guest.disconnect()


def test_guest_room_full_lifecycle(server_auth, Echo, run_worker_loop, wait_for_task):
    """Guest: register -> submit -> worker loop -> bookmark verified."""
    worker = ZnDraw(url=server_auth)
    submitter = ZnDraw(url=server_auth, room=worker.room)
    try:
        worker.jobs.register(Echo, room=worker.room)
        submitter.append(ase.Atoms("H"))

        thread, stop = run_worker_loop(worker, server_auth)

        task_id = submitter.jobs.submit(
            Echo(value="guest_proof"),
            room=submitter.room,
            job_room=worker.room,
        )
        result = wait_for_task(submitter, task_id)
        assert result.status == "completed", f"Task failed: {result.error}"

        assert submitter.bookmarks[0] == "guest_proof"

        stop.set()
        thread.join(timeout=5)
    finally:
        worker.jobs.disconnect()
        worker.disconnect()
        submitter.disconnect()
