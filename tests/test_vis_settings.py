import time

from zndraw import ZnDraw


def test_vis_settings_same_room_same_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user1")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.8

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.5

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.9
    assert client2.settings.studio_lighting.key_light == 0.9


def test_vis_settings_different_room_same_user(server):
    client1 = ZnDraw(url=server, room="room1", user="user1")
    client2 = ZnDraw(url=server, room="room2", user="user1")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.9


def test_vis_settings_same_room_different_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user2")

    client1.settings.studio_lighting.key_light = 0.8
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.8
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client1.settings.studio_lighting.key_light = 0.5
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.7  # Default value

    client2.settings.studio_lighting.key_light = 0.9
    time.sleep(0.1)
    assert client1.settings.studio_lighting.key_light == 0.5
    assert client2.settings.studio_lighting.key_light == 0.9
