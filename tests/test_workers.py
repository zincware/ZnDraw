import pytest
import requests
from conftest import get_jwt_auth_headers

from zndraw import ZnDraw
from zndraw.app.job_manager import JobStatus
from zndraw.extensions import Extension, Category
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections


class ModifierExtension(Extension):
    category = Category.MODIFIER

    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        kwargs["info"].update({"parameter": self.parameter})
        if kwargs["info"].get("raise") is True:
            raise ValueError("Test error")

class RaiseOnParameterExtension(Extension):
    category = Category.MODIFIER

    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        if self.parameter == 1:
            raise ValueError("Parameter cannot be 1")

class SelectionExtension(Extension):
    category = Category.SELECTION

    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        kwargs["info"].update({"parameter": self.parameter})
        if kwargs["info"].get("raise") is True:
            raise ValueError("Test error")


@pytest.mark.parametrize("category", ["modifiers", "selections"])
def test_register_extensions(server, category):
    room = "testroom"
    user = "testuser"
    if category == "modifiers":
        mod = ModifierExtension
        default_keys = set(modifiers.keys())
    elif category == "selections":
        mod = SelectionExtension
        default_keys = set(selections.keys())
    else:
        raise ValueError("Unknown category")
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    # Schema endpoint now returns a list of extension objects
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys

    for default_mod in default_keys:
        response = requests.get(
            f"{server}/api/rooms/{room}/extensions/{category}/{default_mod}/workers"
        )
        assert response.status_code == 200
        response_json = response.json()
        assert response_json == {
            "idleWorkers": [],
            "progressingWorkers": [],
            "queueLength": 0,
            "totalWorkers": 0,
        }

    # connect extension
    vis.register_extension(mod)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys | {mod.__name__}

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idleWorkers": [vis.sid],
        "progressingWorkers": [],
        "queueLength": 0,
        "totalWorkers": 1,
    }

    # disconnect extension via client disconnect
    vis.disconnect()

    # Wait for async cleanup to complete
    import time
    time.sleep(0.5)

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys

    # connect two workers
    vis1 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis2 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis1.register_extension(mod)
    vis2.register_extension(mod)

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert set(response_json["idleWorkers"]) == {vis1.sid, vis2.sid}
    assert response_json["totalWorkers"] == 2

    vis1.disconnect()

    # Wait for async cleanup to complete
    time.sleep(0.5)

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == [vis2.sid]
    assert response_json["totalWorkers"] == 1

    vis2.disconnect()

    # Wait for async cleanup to complete
    time.sleep(0.5)

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == []
    assert response_json["totalWorkers"] == 0


@pytest.mark.parametrize("category", ["modifiers", "selections"])
def test_run_client_extensions(server, category):
    room = "testroom"
    user = "testuser"
    if category == "modifiers":
        mod = ModifierExtension
        default_keys = set(modifiers.keys())
    elif category == "selections":
        mod = SelectionExtension
        default_keys = set(selections.keys())
    else:
        raise ValueError("Unknown category")
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # Client-registered extensions use /private endpoint
    # post request with data = {}
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    # get job status - should be "assigned" since idle worker exists
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__
    assert response_json["status"] == JobStatus.ASSIGNED  # Job assigned to idle worker immediately
    assert response_json["data"] == {"parameter": 42}

    # check all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["id"] == jobId

    # Now we emulate the worker starting the job (new push-based architecture)
    # Step 1: Worker gets job details
    response = requests.get(
        f"{server}/api/jobs/{jobId}",
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId
    assert response_json["data"] == {"parameter": 42}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__

    # Step 2: Worker transitions job to processing
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # get the worker state
    response = requests.get(f"{server}/api/workers/{vis.sid}", headers=get_jwt_auth_headers(server, user))
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": False,
        "currentJob": jobId,
    }

    # check job status again - should be processing since worker picked it up
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.PROCESSING
    assert response_json["worker_id"] == vis.sid

    # now we emulate completing the job by the worker /api/rooms/<string:room_id>/jobs/<string:job_id>/status"
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.COMPLETED

    # get the worker state again
    response = requests.get(f"{server}/api/workers/{vis.sid}", headers=get_jwt_auth_headers(server, user))
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": True,
        "currentJob": None,
    }

    # Check jobs list
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["status"] == JobStatus.COMPLETED

    # let's queue two more job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 44}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId3 = response_json.pop("jobId")
    # jobId3 is the only pending job (jobId2 was assigned and removed from queue)
    assert response_json == {
        "queuePosition": 0,  # No jobs ahead of this one
        "status": "success",
    }

    # Worker picks up next job (new push-based architecture)
    # Step 1: Get job details for jobId2
    response = requests.get(
        f"{server}/api/jobs/{jobId2}",
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId2
    assert response_json["data"] == {"parameter": 43}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__

    # Step 2: Start processing jobId2
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId2}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.PROCESSING
    assert response_json["worker_id"] == vis.sid
    # check job 3 status - should be PENDING since jobId2 is still processing
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId3}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.PENDING
    # check worker state
    response = requests.get(f"{server}/api/workers/{vis.sid}", headers=get_jwt_auth_headers(server, user))
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": False,
        "currentJob": jobId2,
    }

    # complete job 2
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId2}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.COMPLETED
    # check worker state - jobId3 should have been auto-assigned after jobId2 completed
    response = requests.get(f"{server}/api/workers/{vis.sid}", headers=get_jwt_auth_headers(server, user))
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": False,
        "currentJob": jobId3,
    }
    # verify jobId3 is now ASSIGNED
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId3}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.ASSIGNED
    # Worker picks up next job (jobId3 should now be assigned since jobId2 is complete)
    # Step 1: Get job details for jobId3
    response = requests.get(
        f"{server}/api/jobs/{jobId3}",
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId3
    assert response_json["data"] == {"parameter": 44}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__

    # Step 2: Start processing jobId3
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId3}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId3}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # check job list
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 3
    for job in response_json:
        assert job["status"] == JobStatus.COMPLETED


def test_worker_finish_nonstarted_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    # Job is assigned but not processing yet - try to complete it without processing first
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 400
    assert response.json() == {"error": f"Job must be in 'processing' state to complete (current: {JobStatus.ASSIGNED.value})"}

    # Start processing the job first
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200

    # Try to finish with wrong worker id
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": "wrong-worker-id"},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker ID does not match job's worker"}

    # Finish the job with correct worker id
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.COMPLETED
    assert response_json["worker_id"] == vis.sid


def test_worker_fail_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }

    # Job is assigned but not processing yet - try to fail it without processing first
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={
            "status": "failed",
            "workerId": vis.sid,
            "error": "Something went wrong",
        },
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 400
    assert response.json() == {"error": f"Job must be in 'processing' state to fail (current: {JobStatus.ASSIGNED.value})"}

    # Start processing the job first
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200

    # Try to fail with wrong worker id
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={
            "status": "failed",
            "workerId": "wrong-worker-id",
            "error": "Something went wrong",
        },
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker ID does not match job's worker"}

    # Fail the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={
            "status": "failed",
            "workerId": vis.sid,
            "error": "Something went wrong",
        },
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.FAILED
    assert response_json["worker_id"] == vis.sid
    assert response_json["error"] == "Something went wrong"

    # get all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["status"] == JobStatus.FAILED
    assert response_json[0]["error"] == "Something went wrong"
    assert response_json[0]["data"] == {"parameter": 42}
    assert response_json[0]["id"] == jobId


def test_delete_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension

    # First, register the extension to make it available
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # Disconnect the worker so jobs stay in PENDING state (not ASSIGNED)
    vis.socket.disconnect()

    # Submit a job (will be PENDING since no idle workers)
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    # Job is pending (not assigned yet since worker disconnected)
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }

    # delete the job (should succeed since job is PENDING, not ASSIGNED/PROCESSING)
    response = requests.delete(
        f"{server}/api/rooms/{room}/jobs/{jobId}",
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 404
    assert response.json() == {"error": "Job not found"}

    # get all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 0

    # delete a job that does not exist
    response = requests.delete(
        f"{server}/api/rooms/{room}/jobs/non-existing-job-id",
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 404
    assert response.json() == {"error": "Job not found"}


# test_celery_task removed - uses obsolete /api/jobs/next polling endpoint
# Celery workers now use push-based task dispatch via Celery tasks,
# not HTTP polling


@pytest.mark.parametrize("category", ["modifiers", "selections"])
def test_register_extensions_reconnect_with_queue(server, category):
    room = "testroom"
    user = "testuser"
    if category == "modifiers":
        mod = ModifierExtension
        default_keys = set(modifiers.keys())
    elif category == "selections":
        mod = SelectionExtension
        default_keys = set(selections.keys())
    else:
        raise ValueError("Unknown category")
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys | {mod.__name__}

    # submit first job - will be assigned to idle worker
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    response.raise_for_status()
    jobId1 = response.json().pop("jobId")
    # check job status - should be ASSIGNED
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId1}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.ASSIGNED

    # submit second job - worker is busy so this will stay PENDING
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    response.raise_for_status()
    jobId2 = response.json().pop("jobId")
    # check job status - should be PENDING
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.PENDING

    # disconnect client
    requests.post(f"{server}/api/disconnect/{vis.sid}").raise_for_status()

    # check job1 status - should be FAILED since worker had it assigned and disconnected
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId1}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.FAILED

    # check job2 status - should still be PENDING since it was never assigned
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.PENDING

    # Extension should remain in schema because job2 is still pending
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys | {mod.__name__}

    # get available workers for the job
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == []
    assert response_json["totalWorkers"] == 0

    # try submit a job -> should not fail even though there are no workers
    # jobId2 is still pending, so this new job will be position 1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 44}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId3 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 1,  # jobId2 is at position 0
        "status": "success",
    }
    # reconnect client
    assert not vis.socket.connected
    vis.connect()
    assert vis.socket.connected
    # vis.client.sio.sleep(5)  # give some time to process reconnect

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    # Worker should have jobId2 assigned, so not idle
    assert response_json["idleWorkers"] == []
    assert response_json["totalWorkers"] == 1

    # Verify jobId2 got auto-assigned to the reconnected worker
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == JobStatus.ASSIGNED
    assert response_json["worker_id"] == vis.sid

    # disconnect, change auto_pickup_jobs to True and reconnect
    # requests.post(f"{server}/api/disconnect/{vis.sid}").raise_for_status()
    # vis.client.auto_pickup_jobs = True
    # vis.client.connect()
    # assert vis.client.sio.connected
    # vis.client.sio.sleep(10)  # give some time to process reconnect and pickup
    # # check that both jobs are completed
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == JobStatus.COMPLETED
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == JobStatus.COMPLETED


@pytest.mark.parametrize("category", ["modifiers", "selections"])
def test_register_extensions_reconnect_without_queue(server, category):
    room = "testroom"
    user = "testuser"
    if category == "modifiers":
        mod = ModifierExtension
        default_keys = set(modifiers.keys())
    elif category == "selections":
        mod = SelectionExtension
        default_keys = set(selections.keys())
    else:
        raise ValueError("Unknown category")
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names = {ext["name"] for ext in response_json}
    assert returned_names == default_keys | {mod.__name__}
    # disconnect client
    vis.disconnect()
    vis.socket.sio.sleep(1)  # give some time to process disconnect

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    returned_names_after = {ext["name"] for ext in response_json}
    assert (
        returned_names_after == default_keys
    )  # no more jobs in queue, so the modifier should be gone

    # try submit a job -> should fail because extension no longer exists
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/private/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
        headers=get_jwt_auth_headers(server, user),
    )
    assert response.status_code == 404
    response_data = response.json()
    assert "not found" in response_data["error"].lower()
    assert response_data["code"] == "EXTENSION_NOT_FOUND"


def test_submit_task_via_vis_run(server):
    vis = ZnDraw(url=server, room="testroom", user="testuser", auto_pickup_jobs=False)
    vis.register_extension(ModifierExtension)

    vis.run(ModifierExtension(parameter=123))

    response = requests.get(f"{server}/api/rooms/testroom/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["status"] == JobStatus.ASSIGNED
    assert response_json[0]["data"] == {"parameter": 123}

@pytest.mark.parametrize("public", [True, False])
def test_submit_task_twice_via_vis_run(server, public):
    vis1 = ZnDraw(url=server, room="testroom", user="testuser", auto_pickup_jobs=True)
    vis1.register_extension(RaiseOnParameterExtension, public=public)

    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")

    job = vis2.run(RaiseOnParameterExtension(parameter=0), public=public)
    job.wait(timeout=5)
    assert job.status == JobStatus.COMPLETED
    # submit again
    job = vis2.run(RaiseOnParameterExtension(parameter=1), public=public)
    job.wait(timeout=5)
    assert job.status == JobStatus.FAILED
    # # submit again
    job1 = vis2.run(RaiseOnParameterExtension(parameter=0), public=public)
    job2 = vis2.run(RaiseOnParameterExtension(parameter=0), public=public)
    assert job1.status == JobStatus.ASSIGNED
    assert job2.status == JobStatus.PENDING
    job1.wait(timeout=5)
    assert job1.status == JobStatus.COMPLETED
    job2.wait(timeout=5)
    assert job2.status == JobStatus.COMPLETED

@pytest.mark.parametrize("public", [True, False])
def test_submit_task_twice_via_vis_run_two_extensions(server, public):
    w1 = ZnDraw(url=server, room="testroom", user="testuser", auto_pickup_jobs=True)
    w2 = ZnDraw(url=server, room="testroom", user="testuser", auto_pickup_jobs=True)
    w1.register_extension(RaiseOnParameterExtension, public=public)
    w2.register_extension(RaiseOnParameterExtension, public=public)

    vis = ZnDraw(url=server, room="testroom", user="testuser2")
    job1 = vis.run(RaiseOnParameterExtension(parameter=0), public=public)
    job2 = vis.run(RaiseOnParameterExtension(parameter=0), public=public)
    assert job1.status == JobStatus.ASSIGNED
    assert job2.status == JobStatus.ASSIGNED
    job1.wait(timeout=5)
    job2.wait(timeout=5)
    assert job1.status == JobStatus.COMPLETED
    assert job2.status == JobStatus.COMPLETED

    # submit 3
    job1 = vis.run(RaiseOnParameterExtension(parameter=0), public=public)
    job2 = vis.run(RaiseOnParameterExtension(parameter=0), public=public)
    job3 = vis.run(RaiseOnParameterExtension(parameter=0), public=public)
    assert job1.status == JobStatus.ASSIGNED
    assert job2.status == JobStatus.ASSIGNED
    assert job3.status == JobStatus.PENDING 
    job1.wait(timeout=5)
    job2.wait(timeout=5)
    job3.wait(timeout=5)
    assert job1.status == JobStatus.COMPLETED
    assert job2.status == JobStatus.COMPLETED
    assert job3.status == JobStatus.COMPLETED


def test_submit_task_twice_via_vis_register_twice_single_worker(server):
    # Test that a single worker can register the same task public and per room, but only picks up one at a time
    w1 = ZnDraw(url=server, room="testroom", user="testuser", auto_pickup_jobs=True)
    w1.register_extension(RaiseOnParameterExtension, public=True)
    w1.register_extension(RaiseOnParameterExtension, public=False)

    vis = ZnDraw(url=server, room="testroom", user="testuser2")
    job1 = vis.run(RaiseOnParameterExtension(parameter=0), public=True)
    job2 = vis.run(RaiseOnParameterExtension(parameter=0), public=False)
    assert job1.status == JobStatus.ASSIGNED
    assert job2.status == JobStatus.PENDING
    job1.wait(timeout=5)
    assert job1.status == JobStatus.COMPLETED
    job2.wait(timeout=5)
    assert job2.status == JobStatus.COMPLETED