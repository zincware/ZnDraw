import znsocket

from zndraw.frames.utils import dict_to_nested_znsocket


def test_dict_to_nested_znsocket(znsclient):
    data = {"a": 1, "b": {"c": 2, "d": {"e": 3}}}

    nested_dict = dict_to_nested_znsocket(data, "test_key", znsclient)
    assert znsclient.get(nested_dict.key) == {"a": "1", "b": '"znsocket.Dict:test_key.b"'}
    assert znsclient.get(nested_dict["b"].key) == {
        "c": "2",
        "d": '"znsocket.Dict:test_key.b.d"',
    }
    assert znsclient.get(nested_dict["b"]["d"].key) == {"e": "3"}

    assert nested_dict["a"] == 1
    assert nested_dict["b"]["c"] == 2
    assert nested_dict["b"]["d"]["e"] == 3

    assert isinstance(nested_dict, znsocket.Dict)
    assert isinstance(nested_dict["b"], znsocket.Dict)
    assert isinstance(nested_dict["b"]["d"], znsocket.Dict)

    assert nested_dict == {"a": 1, "b": {"c": 2, "d": {"e": 3}}}
