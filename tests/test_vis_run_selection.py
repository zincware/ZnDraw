from zndraw import ZnDraw
from zndraw.extensions import selections

CELERY_TIMEOUT = 2

def test_vis_run_none_selection(server, s22, celery_worker):
    vis = ZnDraw(url=server, room="test", user="tester")
    vis.extend(s22)
    vis.selection = [0, 1, 2]
    assert vis.selection == frozenset([0, 1, 2])
    selections.NoneSelection().run(vis)
    assert vis.selection == frozenset()
    # run via celery
    vis.selection = [0, 1, 2]
    assert vis.selection == frozenset([0, 1, 2])
    vis.run(selections.NoneSelection())
    vis.socket.sio.sleep(CELERY_TIMEOUT)
    assert vis.selection == frozenset()