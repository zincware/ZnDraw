import time

import pytest
import requests

from zndraw.zndraw import ZnDraw


@pytest.fixture
def promoted_template(server, s22):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis1.append(s22[0])
    requests.post(
        f"{server}/api/rooms/s22-0/promote",
        json={"name": "My Template", "description": "A custom template"},
    ).raise_for_status()
    return "s22-0"


def test_rest_list_rooms(server):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1")

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    assert len(rooms) == 2
    assert rooms[0]["id"] == "s22-0"
    assert rooms[0]["template"] == "empty"
    assert rooms[1]["id"] == "s22-1"
    assert rooms[1]["template"] == "empty"


def test_rest_list_templates(server):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1")

    response = requests.get(f"{server}/api/templates")
    assert response.status_code == 200
    templates = response.json()
    assert len(templates) == 1
    assert templates == [
        {
            "id": "empty",
            "name": "Empty Room Template",
            "description": "Empty room template",
        }
    ]


def test_rest_promote_template(server):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1")

    response = requests.post(
        f"{server}/api/rooms/s22-0/promote",
        json={"name": "My Template", "description": "A custom template"},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}

    response = requests.get(f"{server}/api/templates")
    assert response.status_code == 200
    templates = response.json()
    assert len(templates) == 2
    assert templates[0] == {
        "id": "empty",
        "name": "Empty Room Template",
        "description": "Empty room template",
    }
    assert templates[1] == {
        "id": "s22-0",
        "name": "My Template",
        "description": "A custom template",
    }

    with pytest.raises(RuntimeError):  # can't edit the room
        with vis1.lock:
            pass


def test_rest_default_template(server):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1")

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {
        "id": "empty",
        "name": "Empty Room Template",
        "description": "Empty room template",
    }

    response = requests.post(
        f"{server}/api/rooms/s22-0/promote",
        json={"name": "My Template", "description": "A custom template"},
    )
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {
        "id": "empty",
        "name": "Empty Room Template",
        "description": "Empty room template",
    }

    # Set new default template
    response = requests.put(
        f"{server}/api/templates/default", json={"template_id": "s22-0"}
    )
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {
        "id": "s22-0",
        "name": "My Template",
        "description": "A custom template",
    }


def test_rest_promote_template_invalid_room(server):
    response = requests.post(
        f"{server}/api/rooms/invalid-room/promote",
        json={"name": "My Template", "description": "A custom template"},
    )
    assert response.status_code == 404
    assert response.json() == {"error": "Room not found"}


def test_rest_default_template_invalid_template(server):
    response = requests.put(
        f"{server}/api/templates/default", json={"template_id": "nonexistent-template"}
    )
    assert response.status_code == 404
    assert response.json() == {"error": "Template 'nonexistent-template' not found"}


def test_create_from_template(server, promoted_template, s22):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1", template=promoted_template)

    assert len(vis1) == 1
    assert len(vis2) == 1
    assert vis2[0] == s22[0]

    vis2.append(s22[1])
    assert len(vis2) == 2
    assert vis2[0] == s22[0]
    assert vis2[1] == s22[1]

    assert len(vis1) == 1
    assert vis1[0] == s22[0]

    # # Assert the API metadata
    response = requests.get(f"{server}/api/rooms/s22-1")
    assert response.json()["template"] == promoted_template


# Test 2: Setting and using the default template
def test_default_template(server, promoted_template, s22):
    # Set the default
    requests.put(
        f"{server}/api/templates/default", json={"template_id": promoted_template}
    ).raise_for_status()

    # Create a room without specifying a template
    vis4 = ZnDraw(url=server, room="s22-3", user="user1")

    assert len(vis4) == 1
    assert vis4[0] == s22[0]


    response = requests.get(f"{server}/api/rooms/s22-3")
    assert response.json()["template"] == promoted_template


# Test 3: Creating a blank room overrides the default
def test_create_blank_overrides_default(server, promoted_template, s22):
    # Set the default
    requests.put(
        f"{server}/api/templates/default", json={"template_id": promoted_template}
    ).raise_for_status()

    vis5 = ZnDraw(url=server, room="s22-4", user="user1", template=None)
    assert len(vis5) == 0
    response = requests.get(f"{server}/api/rooms/s22-4")
    assert response.json()["template"] == "empty"


def test_default_is_unchanged_after_promote(server, promoted_template):
    """
    Tests that promoting a new template does not change the system default.
    """
    # Create a room without specifying a template. It should use the original default.
    vis3 = ZnDraw(url=server, room="s22-2", user="user1")

    # Assert the data is empty (assuming 'empty' template is blank)
    # TODO: a new room should be empty!
    assert len(vis3) == 0

    # Assert the API metadata
    response = requests.get(f"{server}/api/rooms/s22-2")
    response.raise_for_status()
    assert response.json()["template"] == "empty"


def test_create_from_non_existing_template(server, promoted_template):
    vis2 = ZnDraw(
        url=server, room="s22-1", user="user1", template="nonexistent-template"
    )

    rooms = requests.get(f"{server}/api/rooms/s22-1")
    assert rooms.status_code == 200
    assert rooms.json()["template"] == "empty"  # falls back to empty


def test_rest_get_room(server):
    vis1 = ZnDraw(url=server, room="s22-0", user="user1")
    vis2 = ZnDraw(url=server, room="s22-1", user="user1")

    response = requests.get(f"{server}/api/rooms/s22-0")
    assert response.status_code == 200
    room = response.json()
    assert room["id"] == "s22-0"
    assert room["template"] == "empty"

    response = requests.get(f"{server}/api/rooms/s22-1")
    assert response.status_code == 200
    room = response.json()
    assert room["id"] == "s22-1"
    assert room["template"] == "empty"

    response = requests.get(f"{server}/api/rooms/nonexistent-room")
    assert response.status_code == 404
    assert response.json() == {"error": "Room not found"}
