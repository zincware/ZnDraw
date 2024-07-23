from ase.build import molecule
import numpy.testing as npt
import numpy as np
import pytest
import znjson
from ase.calculators.singlepoint import SinglePointCalculator

from zndraw.utils import ASEConverter
from zndraw import ZnDraw


def test_vectors(server):
    vis = ZnDraw(url=server, token="test_token")

    water = molecule("H2O")
    water.info["vectors"] = [[[1, 0, 0], [0, 1, 0]]]

    vis.append(water)

    assert vis[0].info["vectors"] == [[[1, 0, 0], [0, 1, 0]]]

def test_vectors_list_numpy(server):
    vis = ZnDraw(url=server, token="test_token")

    water = molecule("H2O")
    water.info["vectors"] = [np.array([[1, 0, 0], [0, 1, 0]])]

    vis.append(water)

    assert vis[0].info["vectors"] == [[[1, 0, 0], [0, 1, 0]]]

def test_vectors_numpy(server):
    vis = ZnDraw(url=server, token="test_token")

    water = molecule("H2O")
    water.info["vectors"] = np.array([[[1, 0, 0], [0, 1, 0]]])

    vis.append(water)

    assert vis[0].info["vectors"] == [[[1, 0, 0], [0, 1, 0]]]

def test_vectors_format(server):
    vis = ZnDraw(url=server, token="test_token")

    water = molecule("H2O")
    water.info["vectors"] = [1, 2, 3]

    with pytest.raises(ValueError):
        vis.append(water)
    
    water.info["vectors"] = [[1, 2, 3], [4, 5, 6]]

    with pytest.raises(ValueError):
        vis.append(water)
