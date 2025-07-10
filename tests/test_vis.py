import pytest
import redis
import znjson
from ase.build import molecule

from zndraw import ZnDraw, ZnDrawLocal
from zndraw.converter import ASEConverter


@pytest.fixture
def full(server):
    return ZnDraw(url=server, token="test_token")


@pytest.fixture
def local(server):
    r = redis.Redis.from_url("redis://localhost:6379/0")

    return ZnDrawLocal(url=server, token="test_token", r=r)


@pytest.mark.parametrize("ref", ["full", "local"])
def test_append_atoms(ref, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    water = molecule("H2O")
    vis.append(water)

    assert vis[-1] == water


@pytest.mark.parametrize("ref", ["full", "local"])
def test_append_dump(ref, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)

    water = molecule("H2O")
    vis.append(water)

    assert vis[-1] == water


@pytest.mark.parametrize("ref", ["full", "local"])
def test_append_faulty(ref, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)

    water = molecule("H2O")
    data = ASEConverter().encode(water)
    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.append(data)
    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.extend(3.14)


@pytest.mark.parametrize("ref", ["full", "local"])
def test_setitem_atoms(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    water = molecule("H2O")
    vis[0] = water
    assert vis[0] == water

    vis[[1, 2]] = [water, water]
    assert vis[[0, 1, 2]] == [water, water, water]


@pytest.mark.parametrize("ref", ["local", "full"])
def test_setitem_dump(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    water = molecule("H2O")
    vis[0] = water
    assert vis[0] == molecule("H2O")

    vis[[1, 2]] = [water, water]
    assert vis[[0, 1, 2]] == [water, water, water]


@pytest.mark.parametrize("ref", ["full", "local"])
def test_setitem_faulty(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)

    water = molecule("H2O")
    data = ASEConverter().encode(water)
    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis[0] = data
    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.extend(3.14)


@pytest.mark.parametrize("ref", ["full", "local"])
def test_extend_atoms(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    assert vis[:] == s22


@pytest.mark.parametrize("ref", ["local", "full"])
def test_extend_dump(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)
    assert vis[:] == s22


@pytest.mark.parametrize("ref", ["full", "local"])
def test_extend_faulty(ref, request, s22):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)
    vis.extend(s22)

    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.extend(znjson.dumps(s22, cls=znjson.ZnEncoder.from_converters(ASEConverter)))

    data = [ASEConverter().encode(s) for s in s22]

    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.extend(data)

    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.extend(3.14)
