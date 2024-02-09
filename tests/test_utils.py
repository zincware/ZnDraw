import pytest
from ase.build import molecule
from znframe import Frame

from zndraw.data import RoomGetData
from zndraw.utils import (
    estimate_max_batch_size_for_socket,
    typecast_kwargs,
    wrap_and_check_index,
)


def test_typecast_kwargs():
    @typecast_kwargs
    def func(**data: RoomGetData):
        return data

    data = func(points=True, bookmarks=True)
    assert isinstance(data, RoomGetData)
    assert data.points == True
    assert data.bookmarks == True
    assert data.step == False


def test_typecast_kwargs_exception():
    @typecast_kwargs
    def func(text: str, **data: RoomGetData):
        return data

    with pytest.raises(TypeError):
        data = func(text="Hello World", points=True, bookmarks=True)


@pytest.mark.parametrize(
    "inp, expected",
    [
        ((0, 10), [0]),
        ((-1, 10), [9]),
        (([0, 1, 2], 10), [0, 1, 2]),
        ((slice(0, 10), 10), list(range(0, 10))),
        ((slice(0, 10, 2), 10), list(range(0, 10, 2))),
    ],
)
def test_wrap_and_check_index(inp, expected):
    assert wrap_and_check_index(*inp) == expected


@pytest.mark.parametrize(
    "inp",
    [
        ((10, 10)),
        ((-11, 10)),
        (([0, 1, 2], 2)),
    ],
)
def test_wrap_and_check_index_raises_exceptions(inp):
    with pytest.raises(IndexError):
        wrap_and_check_index(*inp)


def test_estimate_max_chunk_size():
    small = [Frame.from_atoms(molecule("H2O"))]
    big = [Frame.from_atoms(molecule("C60"))]
    assert estimate_max_batch_size_for_socket(
        small
    ) > estimate_max_batch_size_for_socket(big)
    assert estimate_max_batch_size_for_socket(
        small + big
    ) == estimate_max_batch_size_for_socket(big)
