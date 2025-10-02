import time
from zndraw import Client
import requests
import pytest

def test_list_rooms(server):
    vis1 = Client(url=server, room="s22-0", user="user1")
    vis1.connect()
    vis2 = Client(url=server, room="s22-1", user="user1")
    vis2.connect()

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    assert len(rooms) == 2
    assert rooms[0]["id"] == "s22-0"
    assert rooms[1]["id"] == "s22-1"

def test_list_templates(server):
    vis1 = Client(url=server, room="s22-0", user="user1")
    vis1.connect()
    vis2 = Client(url=server, room="s22-1", user="user1")
    vis2.connect()

    response = requests.get(f"{server}/api/templates")
    assert response.status_code == 200
    templates = response.json()
    assert len(templates) == 1
    assert templates == [{"id": "empty", "name": "Empty Room Template", "description": "Empty room template"}]

def test_promote_template(server):
    vis1 = Client(url=server, room="s22-0", user="user1")
    vis1.connect()
    vis2 = Client(url=server, room="s22-1", user="user1")
    vis2.connect()

    response = requests.post(f"{server}/api/rooms/s22-0/promote", json={"name": "My Template", "description": "A custom template"})
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}
    
    response = requests.get(f"{server}/api/templates")
    assert response.status_code == 200
    templates = response.json()
    assert len(templates) == 2
    assert templates[0] == {"id": "empty", "name": "Empty Room Template", "description": "Empty room template"}
    assert templates[1] == {"id": "s22-0", "name": "My Template", "description": "A custom template"}
    
    with pytest.raises(RuntimeError): # can't edit the room
        with vis1.lock:
            pass

def test_default_template(server):
    vis1 = Client(url=server, room="s22-0", user="user1")
    vis1.connect()
    vis2 = Client(url=server, room="s22-1", user="user1")
    vis2.connect()

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {"id": "empty", "name": "Empty Room Template", "description": "Empty room template"}

    response = requests.post(f"{server}/api/rooms/s22-0/promote", json={"name": "My Template", "description": "A custom template"})
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {"id": "empty", "name": "Empty Room Template", "description": "Empty room template"}

    # Set new default template
    response = requests.put(f"{server}/api/templates/default", json={"template_id": "s22-0"})
    assert response.status_code == 200
    assert response.json() == {"status": "ok"}

    response = requests.get(f"{server}/api/templates/default")
    assert response.status_code == 200
    assert response.json() == {"id": "s22-0", "name": "My Template", "description": "A custom template"}

def test_promote_template_invalid_room(server):
    response = requests.post(f"{server}/api/rooms/invalid-room/promote", json={"name": "My Template", "description": "A custom template"})
    assert response.status_code == 404
    assert response.json() == {"error": "Room not found"}

def test_default_template_invalid_template(server):
    response = requests.put(f"{server}/api/templates/default", json={"template_id": "nonexistent-template"})
    assert response.status_code == 404
    assert response.json() == {"error": "Template 'nonexistent-template' not found"}
