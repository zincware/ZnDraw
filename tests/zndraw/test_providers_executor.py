"""Unit tests for InternalProviderExecutor."""

import json
import time
from pathlib import Path

import httpx
import pytest
from httpx import MockTransport, Response

from zndraw import ZnDraw
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


def _poll_read_provider(
    vis: ZnDraw,
    full_name: str,
    params: dict[str, str],
    *,
    timeout: float = 15,
    interval: float = 0.5,
) -> tuple[int, bytes]:
    """Poll a provider read until a non-504 response arrives, return (status, body)."""
    from urllib.parse import urlencode

    deadline = time.time() + timeout
    while time.time() < deadline:
        qs = "?" + urlencode(params)
        resp = vis.api.http.get(
            f"{vis.api.base_url}/v1/joblib/rooms/{vis.room}/providers/{full_name}{qs}",
            headers=vis.api.get_headers(),
        )
        if resp.status_code != 504:
            return resp.status_code, resp.content
        time.sleep(interval)
    pytest.fail(f"Provider read for {params} kept timing out after {timeout}s")


def test_internal_filesystem_provider_surfaces_error(server_factory):
    """A provider read that raises must surface as a 4xx RFC 9457 problem+json,
    not a 504 timeout.
    """
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": ".",
        }
    )
    vis = ZnDraw(url=instance.url)
    try:
        status, body = _poll_read_provider(
            vis,
            "@internal:filesystem:FilesystemRead",
            {"path": "/this/does/not/exist/9k"},
        )
        assert status in (400, 404, 422), f"got {status}: {body!r}"
        parsed = json.loads(body)
        # Executor now emits RFC 9457 ProblemDetail, not ad-hoc {"error","type"}
        assert "title" in parsed, parsed
        assert "status" in parsed, parsed
        assert "detail" in parsed, parsed
        assert "FileNotFoundError" in parsed["detail"], parsed
    finally:
        vis.disconnect()


def test_executor_timeout_from_settings(server_factory):
    """Server boots with a custom executor timeout."""
    instance = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": ".",
            "ZNDRAW_SERVER_PROVIDER_EXECUTOR_TIMEOUT": "5",
        }
    )
    with httpx.Client(base_url=instance.url, timeout=10.0) as client:
        r = client.get("/v1/health")
    assert r.status_code == 200
