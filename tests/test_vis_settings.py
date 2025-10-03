import time

from zndraw import ZnDraw


def test_vis_settings_same_room_same_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user1")

    client1.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is True
    assert client2.settings.scene.show_floor is True

    client1.settings.scene.show_floor = False
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is False
    assert client2.settings.scene.show_floor is False

    client2.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is True
    assert client2.settings.scene.show_floor is True


def test_vis_settings_different_room_same_user(server):
    client1 = ZnDraw(url=server, room="room1", user="user1")
    client2 = ZnDraw(url=server, room="room2", user="user1")

    client1.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is True
    assert client2.settings.scene.show_floor is False

    client1.settings.scene.show_floor = False
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is False
    assert client2.settings.scene.show_floor is False

    client2.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is False
    assert client2.settings.scene.show_floor is True


def test_vis_settings_same_room_different_user(server):
    client1 = ZnDraw(url=server, room="testroom", user="user1")
    client2 = ZnDraw(url=server, room="testroom", user="user2")

    client1.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is True
    assert client2.settings.scene.show_floor is False

    client1.settings.scene.show_floor = False
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is False
    assert client2.settings.scene.show_floor is False

    client2.settings.scene.show_floor = True
    time.sleep(0.1)
    assert client1.settings.scene.show_floor is False
    assert client2.settings.scene.show_floor is True
