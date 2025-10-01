from zndraw.extensions import Extension, ExtensionType
from zndraw.zndraw import ZnDraw
from zndraw.extensions.modifiers import modifiers
from zndraw.extensions.selections import selections
from zndraw.utils import atoms_to_dict, update_colors_and_radii
import rdkit2ase
import time
import requests
import json
import pytest

class ModifierExtension(Extension):
    category = ExtensionType.MODIFIER

    parameter: int

class SelectionExtension(Extension):
    category = ExtensionType.SELECTION

    parameter: int


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
    vis = ZnDraw(url=server, room=room, user=user)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys

    for default_mod in default_keys:
        response = requests.get(f"{server}/api/rooms/{room}/extensions/{category}/{default_mod}/workers")
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

    response = requests.get(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers")
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
    vis1 = ZnDraw(url=server, room=room, user=user)
    vis2 = ZnDraw(url=server, room=room, user=user)
    vis1.register_extension(mod)
    vis2.register_extension(mod)

    response = requests.get(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers")
    assert response.status_code == 200
    response_json = response.json()
    assert set(response_json["idleWorkers"]) == {vis1.client.sid, vis2.client.sid}
    assert response_json["totalWorkers"] == 2

    vis1.client.sio.disconnect()
    response = requests.get(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers")
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["idleWorkers"] == [vis2.client.sid]
    assert response_json["totalWorkers"] == 1

    vis2.client.sio.disconnect()
    response = requests.get(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}/workers")
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
    vis = ZnDraw(url=server, room=room, user=user)
    vis.register_extension(mod)

    # /api/rooms/${roomId}/extensions/${category}/${extension}?userId=${userId}
    # post request with data = {}
    response = requests.post(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}?userId={user}", json={"parameter": 42})
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
    response = requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid})
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

    # now we emulate completing the job by the worker /api/rooms/<string:room_id>/jobs/<string:job_id>/complete"
    response = requests.post(f"{server}/api/rooms/{room}/jobs/{jobId}/complete", json={"workerId": vis.client.sid})
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
    response = requests.post(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}?userId={user}", json={"parameter": 43})
    assert response.status_code == 200
    response_json = response.json()
    jobId2 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 0,
        "status": "success",
    }
    response = requests.post(f"{server}/api/rooms/{room}/extensions/{category}/{mod.__name__}?userId={user}", json={"parameter": 44})
    assert response.status_code == 200
    response_json = response.json()
    jobId3 = response_json.pop("jobId")
    assert response_json == {
        "queuePosition": 1,
        "status": "success",
    }

    # pick up next job
    response = requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid})
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
    response = requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid})
    assert response.status_code == 400
    assert response.json() == {"error": "Worker is not idle"}

    # complete job 2
    response = requests.post(f"{server}/api/rooms/{room}/jobs/{jobId2}/complete", json={"workerId": vis.client.sid})
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
    response = requests.post(f"{server}/api/rooms/{room}/jobs/next", json={"workerId": vis.client.sid})
    assert response.status_code == 200
    response_json = response.json()
    assert response_json["jobId"] == jobId3
    assert response_json["data"] == {"parameter": 44}
    assert response_json["category"] == category
    assert response_json["extension"] == mod.__name__ 
    # finish the job
    response = requests.post(f"{server}/api/rooms/{room}/jobs/{jobId3}/complete", json={"workerId": vis.client.sid})
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



    # # pick up job 2
    # response = requests.post(f"{server}/api/rooms/{room}/jobs/{jobId2}/start", json={"workerId": vis.client.sid})
    # assert response.status_code == 200
    # assert response.json() == {"status": "success"}
    # # check job status again
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId2}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == "running"
    # assert response_json["worker_id"] == vis.client.sid
    # # check job 3 status again
    # response = requests.get(f"{server}/api/rooms/{room}/jobs/{jobId3}")
    # assert response.status_code == 200
    # response_json = response.json()
    # assert response_json["status"] == "queued"
    # assert response_json == {}
    # # assert response_json["queuePosition"] == 0  # position should have updated

    
