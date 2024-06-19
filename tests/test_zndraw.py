import pytest

from zndraw import ZnDraw
import socketio.exceptions


@pytest.fixture
def lst():
    return []


@pytest.fixture
def vis(server):
    return ZnDraw(url=server, token="test_token")


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_zndraw(ref, request):
    """Test the server fixture."""
    obj = request.getfixturevalue(ref)
    assert len(obj) == 0

def test_zndraw_no_connection():
    with pytest.raises(socketio.exceptions.ConnectionError, match="Unable to connect to ZnDraw server at 'ws://localhost:8080'"):
        ZnDraw(url="http://localhost:8080", token="test_token")


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_extend(ref, s22, request):
    """Test the server fixture."""
    obj = request.getfixturevalue(ref)
    obj.extend(s22)
    assert len(obj) == len(s22)


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_getitem(ref, s22, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    assert vis[0] == s22[0]
    assert vis[10] == s22[10]

    with pytest.raises(IndexError):
        vis[100]


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_insert(ref, s22, water, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    vis.insert(0, water)
    assert len(vis) == len(s22) + 1
    assert vis[0] == water
    vis.insert(10, water)
    assert len(vis) == len(s22) + 2
    assert vis[10] == water


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_setitem(ref, s22, water, request):
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    vis[0] = water
    assert vis[0] == water
    vis[10] = water
    assert vis[10] == water

    assert len(vis) == len(s22)


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_setitem_slice(ref, s22, water, request):
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    vis[0:5] = [water] * 5
    assert vis[0:5] == [water] * 5
    assert vis[5:] == s22[5:]


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_delitem(ref, s22, request):
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    del vis[0]
    assert len(vis) == len(s22) - 1
    del vis[10]
    assert len(vis) == len(s22) - 2


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_delitem_slice(ref, server, s22, request):
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    del vis[0:5]
    assert len(vis) == len(s22) - 5
    del vis[5:10]
    assert len(vis) == len(s22) - 10
    del vis[5:]
    assert len(vis) == 5


@pytest.mark.parametrize("ref", ["lst", "vis"])
def test_append(ref, s22, water, request):
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    vis.append(water)
    assert len(vis) == len(s22) + 1
    assert vis[-1] == water
