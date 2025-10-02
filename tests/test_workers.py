import json
import time

import pytest
import rdkit2ase
import requests

from zndraw.extensions import Extension, ExtensionType
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections
from zndraw.utils import atoms_to_dict, update_colors_and_radii
from zndraw.zndraw import ZnDraw


class ModifierExtension(Extension):
    category = ExtensionType.MODIFIER

    parameter: int

    def run(self, vis: ZnDraw, **kwargs):
        kwargs["info"].update({"parameter": self.parameter})
        if kwargs["info"].get("raise") is True:
            raise ValueError("Test error")


class SelectionExtension(Extension):
    category = ExtensionType.SELECTION

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
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys

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
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys | {mod.__name__}

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idleWorkers": [vis.client.sid],
        "progressingWorkers": [],
        "queueLength": 0,
        "totalWorkers": 1,
    }

    # disconnect extension via client disconnect
    vis.client.sio.disconnect()

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys

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
    assert set(response_json["idleWorkers"]) == {vis1.client.sid, vis2.client.sid}
    assert response_json["totalWorkers"] == 2

    vis1.client.sio.disconnect()
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == [vis2.client.sid]
    assert response_json["totalWorkers"] == 1

    vis2.client.sio.disconnect()
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

    # /api/rooms/${roomId}/extensions/${category}/${extension}/submit
    # post request with data = {}
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    # get job status
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__
    assert response_json["status"] == "queued"
    assert response_json["data"] == {"parameter": 42}

    # check all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["id"] == jobId

    # /api/jobs/next?worker_id=<worker_id>
    # now we emulate picking up the job by the worker
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId
    assert response_json["data"] == {"parameter": 42}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__

    # get the worker state
    response = requests.get(f"{server}/api/workers/{vis.client.sid}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": False,
        "currentJob": jobId,
    }

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "running"
    assert response_json["worker_id"] == vis.client.sid

    # now we emulate completing the job by the worker /api/rooms/<string:room_id>/jobs/<string:job_id>/status"
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "completed"

    # get the worker state again
    response = requests.get(f"{server}/api/workers/{vis.client.sid}")
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
    assert response_json[0]["status"] == "completed"

    # let's queue two more job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 44}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId3 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 1,
        "status": "success",
    }

    # pick up next job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId2
    assert response_json["data"] == {"parameter": 43}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "running"
    assert response_json["worker_id"] == vis.client.sid
    # check job 3 status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId3}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "queued"
    # check worker state
    response = requests.get(f"{server}/api/workers/{vis.client.sid}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": False,
        "currentJob": jobId2,
    }
    # Try requesting next job while one is running
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker is not idle"}

    # complete job 2
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId2}/status",
        json={"status": "completed", "workerId": vis.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "completed"
    # check worker state
    response = requests.get(f"{server}/api/workers/{vis.client.sid}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": True,
        "currentJob": None,
    }
    # pick up next job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId3
    assert response_json["data"] == {"parameter": 44}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId3}/status",
        json={"status": "completed", "workerId": vis.client.sid},
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
        assert job["status"] == "completed"


def test_run_different_client_different_extensions(server):
    room = "testroom"
    user = "testuser"
    mod1 = ModifierExtension
    mod2 = SelectionExtension
    vis1 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis2 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis1.register_extension(mod1)
    vis2.register_extension(mod2)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod1.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod1.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    # queue job for mod2
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/selections/{mod2.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
    )
    assert response.status_code == 200
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/selections/{mod2.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
    )
    assert response.status_code == 200

    # assert that there are 4 jobs in total
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 4

    # pick up job for vis1
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis1.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json["jobId"]
    assert response_json["category"] == "modifiers"
    assert response_json["extension"] == mod1.__name__
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis1.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # pick up next job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis1.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json["jobId"]
    assert response_json["category"] == "modifiers"
    assert response_json["extension"] == mod1.__name__
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId2}/status",
        json={"status": "completed", "workerId": vis1.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}
    # try pick up another job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis1.client.sid}
    )
    assert response.status_code == 400
    assert response.json() == {"error": "No jobs available"}

    # queue another job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod1.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId3 = response_json.pop("jobId")
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis1.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId3
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId3}/status",
        json={"status": "completed", "workerId": vis1.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # pick up job for vis2
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis2.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json["jobId"]
    assert response_json["category"] == "selections"
    assert response_json["extension"] == mod2.__name__
    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis2.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # list all jobs, there should be one queued job and 2+1+1 completed jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 5
    completed_jobs = [job for job in response_json if job["status"] == "completed"]
    queued_jobs = [job for job in response_json if job["status"] == "queued"]
    assert len(completed_jobs) == 4
    assert len(queued_jobs) == 1


def test_run_different_client_same_extensions(server):
    room = "testroom"
    user = "testuser"
    vis1 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis2 = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis1.register_extension(ModifierExtension)
    vis2.register_extension(ModifierExtension)

    # queue 4 job
    for idx in range(4):
        _ = requests.post(
            f"{server}/api/rooms/{room}/extensions/modifiers/{ModifierExtension.__name__}/submit",
            json={"data": {"parameter": 42 + idx}, "userId": user},
        )

    # pick up job for vis1
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis1.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId1 = response_json["jobId"]
    assert response_json["category"] == "modifiers"
    assert response_json["extension"] == ModifierExtension.__name__
    # pick up job for vis2
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis2.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json["jobId"]
    assert response_json["category"] == "modifiers"
    assert response_json["extension"] == ModifierExtension.__name__

    # check jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 4
    running_jobs = [job for job in response_json if job["status"] == "running"]
    queued_jobs = [job for job in response_json if job["status"] == "queued"]
    assert len(running_jobs) == 2
    assert len(queued_jobs) == 2
    # finish the jobs for both
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId1}/status",
        json={"status": "completed", "workerId": vis1.client.sid},
    )
    assert response.status_code == 200
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId2}/status",
        json={"status": "completed", "workerId": vis2.client.sid},
    )
    assert response.status_code == 200

    # check jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 4
    completed_jobs = [job for job in response_json if job["status"] == "completed"]
    queued_jobs = [job for job in response_json if job["status"] == "queued"]
    assert len(completed_jobs) == 2
    assert len(queued_jobs) == 2


def test_worker_finish_nonstarted_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    # complete job without starting it
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.client.sid},
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Job is not running"}
    # pick up the job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId
    assert response_json["data"] == {"parameter": 42}

    # finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": "wrong-worker-id"},
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker ID does not match job's worker ID"}

    # finish the job with correct worker id
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": vis.client.sid},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "completed"
    assert response_json["worker_id"] == vis.client.sid


def test_worker_fail_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }

    # fail the job without starting it
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "failed", "workerId": vis.client.sid, "error": "Something went wrong"},
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Job is not running"}

    # pick up the job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId
    assert response_json["data"] == {"parameter": 42}

    # fail the job with wrong worker id
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "failed", "workerId": "wrong-worker-id", "error": "Something went wrong"},
    )
    assert response.status_code == 400
    assert response.json() == {"error": "Worker ID does not match job's worker ID"}

    # fail the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "failed", "workerId": vis.client.sid, "error": "Something went wrong"},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "failed"
    assert response_json["worker_id"] == vis.client.sid
    assert response_json["error"] == "Something went wrong"

    # get all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["status"] == "failed"
    assert response_json[0]["error"] == "Something went wrong"
    assert response_json[0]["data"] == {"parameter": 42}
    assert response_json[0]["id"] == jobId


def test_delete_job(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user, auto_pickup_jobs=False)
    vis.register_extension(mod)

    # queue job for mod1
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }

    # delete the job
    response = requests.delete(f"{server}/api/rooms/{room}/jobs/{jobId}")
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
    response = requests.delete(f"{server}/api/rooms/{room}/jobs/non-existing-job-id")
    assert response.status_code == 404
    assert response.json() == {"error": "Job not found"}


def test_worker_pickup_task(server):
    room = "testroom"
    user = "testuser"
    mod = ModifierExtension
    vis = ZnDraw(url=server, room=room, user=user)
    shared_dict = {}

    vis.register_extension(mod, run_kwargs={"info": shared_dict})
    # submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    vis.client.sio.sleep(1)  # give some time to pick up the job and run it
    # assert shared_dict == {"parameter": 42}
    # get job status
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "completed"
    assert response_json["worker_id"] == vis.client.sid
    assert response_json["data"] == {"parameter": 42}
    assert response_json["error"] == ""

    shared_dict["raise"] = True
    # submit another job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
        json={"data": {"parameter": 43}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    vis.client.sio.sleep(1)  # give some time to pick up the job and run it
    # get job status
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "failed"
    assert response_json["worker_id"] == vis.client.sid
    assert response_json["data"] == {"parameter": 43}
    assert response_json["error"] == "Test error"

    # check worker state
    response = requests.get(f"{server}/api/workers/{vis.client.sid}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json == {
        "idle": True,
        "currentJob": None,
    }

    # submit two jobs to finish
    shared_dict["raise"] = False
    for idx in range(2):
        response = requests.post(
            f"{server}/api/rooms/{room}/extensions/modifiers/{mod.__name__}/submit",
            json={"data": {"parameter": 44 + idx}, "userId": user},
        )
        assert response.status_code == 200
        response_json = response.json()
        jobId = response_json.pop("jobId")
        assert response_json == {
            "queuePosition": 0 + idx,
            "status": "success",
        }

    vis.client.sio.sleep(3)  # give some time to pick up the job and run it
    # get all jobs
    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 4
    completed_jobs = [job for job in response_json if job["status"] == "completed"]
    failed_jobs = [job for job in response_json if job["status"] == "failed"]
    assert len(completed_jobs) == 3
    assert len(failed_jobs) == 1


def test_celery_task(server):
    room = "testroom"
    user = "testuser"
    mod_name = next(iter(modifiers.keys()))

    # Create a celery job by calling a server-side modifier extension
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/modifiers/{mod_name}/submit",
        json={"data": {}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }

    # Verify the job was created with celery provider
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "queued"
    assert response_json["provider"] == "celery"

    # Have celery-worker fetch the job
    response = requests.post(
        f"{server}/api/rooms/{room}/jobs/next", json={"workerId": "celery-worker"}
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId
    assert response_json["category"] == "modifiers"
    assert response_json["extension"] == mod_name
    assert response_json["data"] == {}
    assert response_json["status"] == "running"

    # Finish the job
    response = requests.put(
        f"{server}/api/rooms/{room}/jobs/{jobId}/status",
        json={"status": "completed", "workerId": "celery-worker"},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "success"}

    response = requests.get(f"{server}/api/rooms/{room}/jobs")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, list)
    assert len(response_json) == 1
    assert response_json[0]["status"] == "completed"


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
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys | {mod.__name__}

    # submit a job
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    response.raise_for_status()
    jobId = response.json().pop("jobId")
    # check job status
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "queued"

    # disconnect client
    requests.post(f"{server}/api/disconnect/{vis.client.sid}").raise_for_status()

    # check job status again
    response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["status"] == "queued"

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys | {mod.__name__}

    # get available workers for the job
    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == []
    assert response_json["totalWorkers"] == 0

    # try submit a job -> should not fail
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 1,
        "status": "success",
    }
    # reconnect client
    assert not vis.client.sio.connected
    vis.client.connect()
    assert vis.client.sio.connected
    # vis.client.sio.sleep(5)  # give some time to process reconnect

    response = requests.get(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers"
    )
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == [vis.client.sid]
    assert response_json["totalWorkers"] == 1

    # disconnect, change auto_pickup_jobs to True and reconnect
    # requests.post(f"{server}/api/disconnect/{vis.client.sid}").raise_for_status()
    # vis.client.auto_pickup_jobs = True
    # vis.client.connect()
    # assert vis.client.sio.connected
    # vis.client.sio.sleep(10)  # give some time to process reconnect and pickup
    # # check that both jobs are completed
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == "completed"
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == "completed"



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
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys | {mod.__name__}
    # disconnect client
    vis.client.sio.disconnect()
    vis.client.sio.sleep(1)  # give some time to process disconnect

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys # no more jobs in queue, so the modifier should be gone

    # try submit a job -> should fail
    response = requests.post(
        f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/submit",
        json={"data": {"parameter": 42}, "userId": user},
    )
    assert response.status_code == 400
    assert response.json() == {"error": f"No workers available for extension {mod.__name__}"}

# TODO: test job in queue, modifier disconnected
# - submit a new job -> should be ok
# get schema should include the modifier
# TODO: test, no job in queue, modifier disconnected
# - get schema should not include the modifier
# - submit a new job -> should fail
