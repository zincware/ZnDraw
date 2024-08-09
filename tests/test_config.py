from zndraw import ZnDraw
from zndraw.config import ZnDrawConfig
import znsocket
import os
import redis


def test_config_defaults(server):
    vis = ZnDraw(url=server, token="test_token")
    ref_config = ZnDrawConfig(vis=None)

    vis.config.to_dict() == ref_config.to_dict()


def test_config_modify_arrows(server):
    room = "test_config_arrows"
    vis = ZnDraw(url=server, token=room)
    vis.config.arrows.normalize = False
    assert vis.config.arrows.normalize is False

    vis.config.arrows.normalize = True
    assert vis.config.arrows.normalize is True

    r = redis.Redis.from_url(os.environ["FLASK_STORAGE"])
    key = f"room:{room}:config"
    config = znsocket.Dict(r, key)
    assert config["arrows"]["normalize"] is True


def test_config_modify_scene(server):
    room = "test_config_scene"
    vis = ZnDraw(url=server, token=room)

    vis.config.scene.fps = 30
    assert vis.config.scene.fps == 30
    vis.config.scene.fps = 60
    assert vis.config.scene.fps == 60

    r = redis.Redis.from_url(os.environ["FLASK_STORAGE"])
    key = f"room:{room}:config"
    config = znsocket.Dict(r, key)
    assert config["scene"]["fps"] == 60
    assert config["scene"]["material"] == "MeshStandardMaterial"
