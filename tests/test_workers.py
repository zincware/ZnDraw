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

class SelectionExtension(Extension):
    category = ExtensionType.SELECTION


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
    # connect extension
    vis.register_extension(mod)
    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys | {mod.__name__}
    # disconnect extension via client disconnect
    vis.client.sio.disconnect()

    response = requests.get(f"{server}/api/rooms/{room}/schema/{category}")
    assert response.status_code == 200
    response_json = response.json()
    assert isinstance(response_json, dict)
    assert set(response_json.keys()) == default_keys

