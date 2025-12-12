from zndraw import ZnDraw


def test_vis_append_v1v2(server, s22):
    v1 = ZnDraw(url=server, room="testroom", user="testuser")
    v2 = ZnDraw(url=server, room="testroom", user="testuser2")
    assert len(v1) == 0
    assert len(v2) == 0

    for atoms in s22:
        v1.append(atoms)

    assert len(v1) == 22
    assert len(v2) == 22

    for _ in range(5):
        del v2[0]

    assert len(v1) == 17
    assert len(v2) == 17
