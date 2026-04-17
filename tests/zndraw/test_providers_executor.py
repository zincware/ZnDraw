"""Unit tests for InternalProviderExecutor."""

import json
from pathlib import Path

import pytest
from httpx import MockTransport, Response

from zndraw.providers.executor import InternalProviderExecutor
from zndraw.providers.filesystem import FilesystemRead


@pytest.fixture
def seeded_dir(tmp_path: Path) -> Path:
    (tmp_path / "a.xyz").touch()
    (tmp_path / "b.xyz").touch()
    return tmp_path


async def test_executor_posts_result_for_filesystem(seeded_dir):
    """Executor resolves fsspec handler, runs read(), POSTs bytes."""
    captured: dict = {}

    def handler(request):
        captured["url"] = str(request.url)
        captured["headers"] = dict(request.headers)
        captured["body"] = request.content
        return Response(204)

    transport = MockTransport(handler)

    executor = InternalProviderExecutor(
        base_url="http://test",
        filebrowser_path=str(seeded_dir),
        _transport=transport,
    )

    await executor(
        FilesystemRead,
        params_json=json.dumps({"path": "/"}),
        provider_id="11111111-1111-1111-1111-111111111111",
        request_id="abc",
        token="tok",
    )

    assert captured["url"].endswith(
        "/v1/joblib/providers/11111111-1111-1111-1111-111111111111/results"
    )
    assert captured["headers"]["authorization"] == "Bearer tok"
    assert captured["headers"]["x-request-hash"] == "abc"
    body = json.loads(captured["body"])
    names = {item["name"] for item in body}
    assert names == {"a.xyz", "b.xyz"}
