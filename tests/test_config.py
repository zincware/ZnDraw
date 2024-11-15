from zndraw import ZnDraw

# def test_config_defaults(server):
#     vis = ZnDraw(url=server, token="test_token")
#     ref_config = ZnDrawConfig(vis=None)

#     vis.config == ref_config.to_dict()


def test_config_modify_arrows(server):
    room = "test_config_arrows"
    vis = ZnDraw(url=server, token=room)

    vis.config["arrows"]["normalize"] = False
    assert vis.config["arrows"]["normalize"] is False

    vis.config["arrows"]["normalize"] = True
    assert vis.config["arrows"]["normalize"] is True

    # check other defaults
    assert vis.config["arrows"]["opacity"] == 1.0


# def test_config_replace_znsocket(server):
#     room = "test_config_znsocket"
#     vis = ZnDraw(url=server, token=room)

#     vis.config["scene"] = None
#     assert isinstance(vis.config["scene"], znsocket.Dict)


def test_config_modify_scene(server):
    room = "test_config_scene"
    vis = ZnDraw(url=server, token=room)

    vis.config["scene"]["fps"] = 30
    assert vis.config["scene"]["fps"] == 30
    vis.config["scene"]["fps"] = 60
    assert vis.config["scene"]["fps"] == 60

    # check other defaults
    assert vis.config["scene"]["material"] == "MeshStandardMaterial"
