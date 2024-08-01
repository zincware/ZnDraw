import pytest
import redis
import znjson
from ase.build import molecule

from zndraw import ZnDraw, ZnDrawLocal
from zndraw.utils import ASEConverter


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
    vis.append(znjson.dumps(water, cls=znjson.ZnEncoder.from_converters(ASEConverter)))

    assert vis[-1] == water


@pytest.mark.parametrize("ref", ["full", "local"])
def test_append_faulty(ref, request):
    """Test the server fixture."""
    vis = request.getfixturevalue(ref)

    water = molecule("H2O")
    data = ASEConverter().encode(water)
    with pytest.raises(ValueError, match="Unable to parse provided data object"):
        vis.append(data)
