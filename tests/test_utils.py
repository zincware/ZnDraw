import pytest

from zndraw.data import RoomGetData
from zndraw.utils import typecast_kwargs


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
