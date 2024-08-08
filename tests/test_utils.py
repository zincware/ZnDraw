import numpy as np
import numpy.testing as npt

from zndraw.utils import (
    convert_url_to_http,
    direction_to_euler,
    euler_to_direction,
    parse_url,
)


def test_parse_url():
    assert parse_url("http://example.com") == ("http://example.com", None)
    assert parse_url("http://example.com/path") == ("http://example.com", "path")
    assert parse_url("http://localhost:5000") == ("http://localhost:5000", None)
    assert parse_url("http://localhost:5000/path") == ("http://localhost:5000", "path")
    assert parse_url("http://localhost:5000/path/") == ("http://localhost:5000", "path")


def test_conversion_utils():
    """Test conversion functions"""
    direction = np.array([1, 2, 3])
    direction = direction / np.linalg.norm(direction)
    euler = direction_to_euler(direction)
    new_direction = euler_to_direction(euler)

    npt.assert_allclose(direction, new_direction, atol=1e-6)

    direction = np.array([1, 0, 0])
    euler = direction_to_euler(direction)
    new_direction = euler_to_direction(euler)

    npt.assert_allclose(direction, new_direction, atol=1e-6)
    npt.assert_allclose(euler, [0, 0, 0], atol=1e-6)

    direction = np.array([0, -1, 0])
    euler = direction_to_euler(direction)
    new_direction = euler_to_direction(euler)

    npt.assert_allclose(direction, new_direction, atol=1e-6)


def test_url():
    # safe url
    before_url = "ws://localhost:8000/token/1234"
    url = convert_url_to_http(before_url)
    assert url == "http://localhost:8000/token/1234"

    # unsafe url containing ws in token
    before_url = "ws://localhost:8000/token/eNwsdW5k"
    url = convert_url_to_http(before_url)
    assert url == "http://localhost:8000/token/eNwsdW5k"

    # Now with https
    before_url = "wss://localhost:8000/token/1234"
    url = convert_url_to_http(before_url)
    assert url == "https://localhost:8000/token/1234"

    # unsafe url containing ws in token
    before_url = "wss://localhost:8000/token/eNwsdW5k"
    url = convert_url_to_http(before_url)
    assert url == "https://localhost:8000/token/eNwsdW5k"
