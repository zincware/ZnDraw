from zndraw.utils import parse_url


def test_parse_url():
    assert parse_url("http://example.com") == ("http://example.com", None)
    assert parse_url("http://example.com/path") == ("http://example.com", "path")
    assert parse_url("http://localhost:5000") == ("http://localhost:5000", None)
    assert parse_url("http://localhost:5000/path") == ("http://localhost:5000", "path")
    assert parse_url("http://localhost:5000/path/") == ("http://localhost:5000", "path")
