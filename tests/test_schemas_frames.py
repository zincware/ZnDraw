"""Tests for frame schemas."""

from zndraw.schemas import (
    FrameBulkResponse,
    FrameCreateRequest,
    FrameResponse,
    FrameUpdateRequest,
)


class TestFrameSchemas:
    """Tests for frame Pydantic schemas."""

    def test_frame_response(self) -> None:
        """FrameResponse validates frame data."""
        response = FrameResponse(
            index=0,
            data={"numbers": [1, 2], "positions": [[0, 0, 0], [1, 0, 0]]},
        )
        assert response.index == 0
        assert response.data["numbers"] == [1, 2]

    def test_frame_bulk_response(self) -> None:
        """FrameBulkResponse contains multiple frames."""
        response = FrameBulkResponse(
            frames=[
                {"numbers": [1], "positions": [[0, 0, 0]]},
                {"numbers": [2], "positions": [[1, 0, 0]]},
            ],
            total=10,
            start=0,
            stop=2,
        )
        assert len(response.frames) == 2
        assert response.total == 10

    def test_frame_create_request(self) -> None:
        """FrameCreateRequest validates frame data for creation."""
        request = FrameCreateRequest(
            frames=[{"numbers": [1, 2], "positions": [[0, 0, 0], [1, 0, 0]]}]
        )
        assert len(request.frames) == 1

    def test_frame_update_request(self) -> None:
        """FrameUpdateRequest validates frame update."""
        request = FrameUpdateRequest(
            data={"numbers": [1, 2, 3], "positions": [[0, 0, 0], [1, 0, 0], [2, 0, 0]]}
        )
        assert request.data["numbers"] == [1, 2, 3]

    def test_frame_response_model_dump(self) -> None:
        """FrameResponse can be serialized."""
        response = FrameResponse(index=5, data={"key": "value"})
        dumped = response.model_dump()
        assert dumped == {"index": 5, "data": {"key": "value"}}

    def test_frame_bulk_response_model_dump(self) -> None:
        """FrameBulkResponse can be serialized."""
        response = FrameBulkResponse(
            frames=[{"a": 1}],
            total=100,
            start=0,
            stop=1,
        )
        dumped = response.model_dump()
        assert dumped["total"] == 100
        assert dumped["frames"] == [{"a": 1}]
