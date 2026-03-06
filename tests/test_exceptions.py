import pytest

from zndraw.exceptions import (
    BookmarkNotFound,
    FigureNotFound,
    InvalidCredentials,
    InvalidPayload,
    NotAuthenticated,
    NotInRoom,
    ProblemDetail,
    ProblemException,
    ProblemType,
    RoomNotFound,
    UsernameExists,
    UserNotFound,
    problem_responses,
)


class TestProblemTypeCreate:
    """Tests for ProblemType.create() method."""

    def test_create_with_detail(self) -> None:
        """Test creating a ProblemDetail with a detail message."""
        problem = InvalidCredentials.create(detail="Invalid username or password")

        assert problem.type == "/v1/problems/invalid-credentials"
        assert problem.title == "Unauthorized"
        assert problem.status == 401
        assert problem.detail == "Invalid username or password"
        assert problem.instance is None

    def test_create_with_instance(self) -> None:
        """Test creating a ProblemDetail with an instance URI."""
        problem = UserNotFound.create(
            detail="User not found",
            instance="/v1/users/123",
        )

        assert problem.type == "/v1/problems/user-not-found"
        assert problem.title == "Not Found"
        assert problem.status == 404
        assert problem.detail == "User not found"
        assert problem.instance == "/v1/users/123"

    def test_create_without_optional_fields(self) -> None:
        """Test creating a ProblemDetail without optional fields."""
        problem = NotAuthenticated.create()

        assert problem.type == "/v1/problems/not-authenticated"
        assert problem.title == "Unauthorized"
        assert problem.status == 401
        assert problem.detail is None
        assert problem.instance is None

    @pytest.mark.parametrize(
        ("problem_type", "expected_id", "expected_status"),
        [
            (InvalidCredentials, "invalid-credentials", 401),
            (UsernameExists, "username-exists", 409),
            (NotAuthenticated, "not-authenticated", 401),
            (UserNotFound, "user-not-found", 404),
            (InvalidPayload, "invalid-payload", 400),
            (NotInRoom, "not-in-room", 412),
        ],
    )
    def test_all_problem_types(
        self,
        problem_type: type[ProblemType],
        expected_id: str,
        expected_status: int,
    ) -> None:
        """Test that all problem types create valid ProblemDetails."""
        problem = problem_type.create()

        assert problem.type == f"/v1/problems/{expected_id}"
        assert problem.status == expected_status
        assert isinstance(problem, ProblemDetail)


class TestProblemTypeException:
    """Tests for ProblemType.exception() method."""

    def test_exception_creates_problem_exception(self) -> None:
        """Test that exception() returns a ProblemException."""
        exc = InvalidCredentials.exception(detail="Test error")

        assert isinstance(exc, ProblemException)
        assert exc.problem.type == "/v1/problems/invalid-credentials"
        assert exc.problem.detail == "Test error"

    def test_exception_with_headers(self) -> None:
        """Test that exception() accepts headers."""
        headers = {"WWW-Authenticate": "Bearer"}
        exc = NotAuthenticated.exception(
            detail="Token required",
            headers=headers,
        )

        assert isinstance(exc, ProblemException)
        assert exc.headers == headers
        assert exc.problem.detail == "Token required"

    def test_exception_is_raisable(self) -> None:
        """Test that the exception can be raised and caught."""
        with pytest.raises(ProblemException) as exc_info:
            raise InvalidCredentials.exception(detail="Bad credentials")

        assert exc_info.value.problem.status == 401


class TestProblemDetail:
    """Tests for ProblemDetail model."""

    def test_model_dump_excludes_none(self) -> None:
        """Test that model_dump with exclude_none omits None fields."""
        problem = InvalidCredentials.create(detail="Error")
        data = problem.model_dump(exclude_none=True)

        assert "instance" not in data
        assert data == {
            "type": "/v1/problems/invalid-credentials",
            "title": "Unauthorized",
            "status": 401,
            "detail": "Error",
        }

    def test_model_validate_from_dict(self) -> None:
        """Test that ProblemDetail can be validated from a dict."""
        data = {
            "type": "/v1/problems/invalid-credentials",
            "title": "Unauthorized",
            "status": 401,
            "detail": "Test",
        }
        problem = ProblemDetail.model_validate(data)

        assert problem == InvalidCredentials.create(detail="Test")


class TestSocketIOProblemTypes:
    """Tests for Socket.IO specific problem types."""

    def test_invalid_payload_create(self) -> None:
        """Test InvalidPayload problem type."""
        problem = InvalidPayload.create(detail="Missing room_id")

        assert problem.type == "/v1/problems/invalid-payload"
        assert problem.title == "Bad Request"
        assert problem.status == 400
        assert problem.detail == "Missing room_id"

    def test_not_in_room_create(self) -> None:
        """Test NotInRoom problem type."""
        problem = NotInRoom.create(detail="Not in this room")

        assert problem.type == "/v1/problems/not-in-room"
        assert problem.title == "Precondition Failed"
        assert problem.status == 412
        assert problem.detail == "Not in this room"

    def test_invalid_payload_model_dump_excludes_none(self) -> None:
        """Test that InvalidPayload model_dump with exclude_none omits None fields."""
        problem = InvalidPayload.create(detail="Missing field")
        data = problem.model_dump(exclude_none=True)

        assert "instance" not in data
        assert data == {
            "type": "/v1/problems/invalid-payload",
            "title": "Bad Request",
            "status": 400,
            "detail": "Missing field",
        }

    def test_not_in_room_model_dump_excludes_none(self) -> None:
        """Test that NotInRoom model_dump with exclude_none omits None fields."""
        problem = NotInRoom.create(detail="Not in this room")
        data = problem.model_dump(exclude_none=True)

        assert "instance" not in data
        assert data == {
            "type": "/v1/problems/not-in-room",
            "title": "Precondition Failed",
            "status": 412,
            "detail": "Not in this room",
        }


class TestProblemResponses:
    """Tests for the problem_responses() helper."""

    def test_preserves_all_same_status_descriptions(self) -> None:
        """problem_responses preserves info from all types sharing a status code."""
        result = problem_responses(RoomNotFound, FigureNotFound, BookmarkNotFound)
        # All three are 404 — the result must retain all of them, not just the last
        entry = result[404]
        description = entry["description"]
        assert "room" in description.lower()
        assert "figure" in description.lower()
        assert "bookmark" in description.lower()
