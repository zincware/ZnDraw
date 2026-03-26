"""Tests for run_kwargs support in registration and task claiming."""

from typing import Any, ClassVar

from zndraw_joblib.client import (
    Category,
    ClaimedTask,
    Extension,
    JobManager,
)


class Echo(Extension):
    """Extension that records what run() received."""

    category: ClassVar[Category] = Category.MODIFIER
    value: int = 0

    def run(self, vis: Any, **kwargs: Any) -> None:
        pass


# -- register() accepts run_kwargs -----------------------------------------


def test_register_stores_run_kwargs(api, client):
    """register(cls, run_kwargs={...}) should store kwargs in the registry."""
    manager = JobManager(api)
    model = object()  # non-serializable, worker-local
    manager.register(Echo, run_kwargs={"model": model})

    full_name = "@global:modifiers:Echo"
    assert full_name in manager
    # run_kwargs should be retrievable from the registry entry
    assert manager.get_run_kwargs(full_name) == {"model": model}


def test_register_without_run_kwargs_defaults_empty(api, client):
    """register(cls) without run_kwargs should default to empty dict."""
    manager = JobManager(api)
    manager.register(Echo)

    full_name = "@global:modifiers:Echo"
    assert manager.get_run_kwargs(full_name) == {}


def test_register_decorator_with_run_kwargs(api, client):
    """@manager.register(run_kwargs={...}) decorator form should work."""
    manager = JobManager(api)
    model = object()

    @manager.register(run_kwargs={"model": model})
    class Custom(Extension):
        category: ClassVar[Category] = Category.MODIFIER

        def run(self, vis: Any, **kwargs: Any) -> None:
            pass

    full_name = "@global:modifiers:Custom"
    assert full_name in manager
    assert manager.get_run_kwargs(full_name) == {"model": model}


# -- ClaimedTask carries run_kwargs ----------------------------------------


def test_claimed_task_has_run_kwargs():
    """ClaimedTask should expose run_kwargs attribute."""
    ext = Echo(value=1)
    task = ClaimedTask(
        task_id="abc",
        job_name="@global:modifiers:Echo",
        room_id="room_1",
        extension=ext,
        run_kwargs={"model": "gpt"},
    )
    assert task.run_kwargs == {"model": "gpt"}


def test_claimed_task_run_kwargs_defaults_empty():
    """ClaimedTask without run_kwargs should default to empty dict."""
    ext = Echo(value=1)
    task = ClaimedTask(
        task_id="abc",
        job_name="@global:modifiers:Echo",
        room_id="room_1",
        extension=ext,
    )
    assert task.run_kwargs == {}


# -- claim() populates run_kwargs from registry ----------------------------


def test_claim_populates_run_kwargs(api, client):
    """claim() should attach stored run_kwargs to the ClaimedTask."""
    model = object()
    manager = JobManager(api)
    manager.register(Echo, run_kwargs={"model": model})

    # Submit a task
    resp = client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Echo",
        json={"payload": {"value": 42}},
    )
    assert resp.status_code == 202

    claimed = manager.claim()
    assert claimed is not None
    assert claimed.run_kwargs == {"model": model}
    assert claimed.extension.value == 42


def test_claim_without_run_kwargs_gives_empty(api, client):
    """claim() for a job registered without run_kwargs should give empty dict."""
    manager = JobManager(api)
    manager.register(Echo)

    resp = client.post(
        "/v1/joblib/rooms/room_1/tasks/@global:modifiers:Echo",
        json={"payload": {"value": 7}},
    )
    assert resp.status_code == 202

    claimed = manager.claim()
    assert claimed is not None
    assert claimed.run_kwargs == {}
