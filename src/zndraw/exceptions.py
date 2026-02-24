"""RFC 9457 Problem Details for HTTP APIs."""

import re
from typing import Any, ClassVar, NoReturn

from fastapi import Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from zndraw_joblib.exceptions import ProviderTimeout


def _camel_to_kebab(name: str) -> str:
    """Convert CamelCase to kebab-case."""
    return re.sub(r"(?<!^)(?=[A-Z])", "-", name).lower()


class ProblemType:
    """Base class for defining problem types with documentation."""

    title: ClassVar[str]
    status: ClassVar[int]

    @classmethod
    def problem_id(cls) -> str:
        """Return kebab-case identifier derived from class name."""
        return _camel_to_kebab(cls.__name__)

    @classmethod
    def type_uri(cls) -> str:
        """Return the full type URI for this problem."""
        return f"/v1/problems/{cls.problem_id()}"

    @classmethod
    def openapi_response(
        cls, description: str | None = None
    ) -> dict[int | str, dict[str, Any]]:
        """Generate OpenAPI response entry for this problem type."""
        return {
            cls.status: {
                "model": ProblemDetail,
                "description": description
                or (cls.__doc__.split("\n")[0] if cls.__doc__ else cls.title),
            }
        }

    @classmethod
    def create(
        cls, detail: str | None = None, instance: str | None = None
    ) -> "ProblemDetail":
        """Create a ProblemDetail instance from this type."""
        return ProblemDetail(
            type=cls.type_uri(),
            title=cls.title,
            status=cls.status,
            detail=detail,
            instance=instance,
        )

    @classmethod
    def exception(
        cls,
        detail: str | None = None,
        instance: str | None = None,
        headers: dict[str, str] | None = None,
    ) -> "ProblemException":
        """Create a ProblemException from this type."""
        return ProblemException(cls.create(detail, instance), headers=headers)

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        """Raise appropriate Python exception on the client side.

        Override in subclasses to map to specific exception types
        (e.g. ``IndexError``, ``KeyError``). The base implementation
        raises a generic ``RuntimeError``.
        """
        raise RuntimeError(problem.detail or problem.title)


# =============================================================================
# Problem Type Definitions
# =============================================================================


class InvalidCredentials(ProblemType):
    """The provided username or password is incorrect.

    This error occurs when attempting to authenticate with credentials
    that do not match any existing user account, or when the password
    is incorrect for the given username.
    """

    title: ClassVar[str] = "Unauthorized"
    status: ClassVar[int] = 401

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class UsernameExists(ProblemType):
    """The requested username is already taken.

    This error occurs when attempting to register a new user account
    with a username that already exists in the system. Choose a
    different username and try again.
    """

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class NotAuthenticated(ProblemType):
    """Authentication is required but was not provided.

    This error occurs when accessing a protected resource without
    providing valid authentication credentials.
    """

    title: ClassVar[str] = "Unauthorized"
    status: ClassVar[int] = 401

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class UserNotFound(ProblemType):
    """The requested user does not exist.

    This error occurs when attempting to access a user account
    that has been deleted or never existed.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class RoomNotFound(ProblemType):
    """The requested room does not exist.

    This error occurs when attempting to access a chat room
    that has been deleted or never existed.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class SessionNotFound(ProblemType):
    """The requested session does not exist or does not belong to the user.

    This error occurs when attempting to access a session that either
    does not exist or belongs to a different user.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class NotRoomMember(ProblemType):
    """The user is not a member of this room.

    This error occurs when attempting to access a private room
    without being a member, or when performing member-only actions.
    """

    title: ClassVar[str] = "Forbidden"
    status: ClassVar[int] = 403

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class AlreadyRoomMember(ProblemType):
    """The user is already a member of this room.

    This error occurs when attempting to join a room
    that the user is already a member of.
    """

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class InvalidPayload(ProblemType):
    """The request payload is missing or malformed.

    This error occurs when a required field is missing from the
    request payload or when the payload structure is invalid.
    """

    title: ClassVar[str] = "Bad Request"
    status: ClassVar[int] = 400

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class UnprocessableContent(ProblemType):
    """The request is syntactically valid but semantically invalid.

    This error occurs when a field value is well-formed but not
    an accepted value (e.g. an unknown preset name).
    """

    title: ClassVar[str] = "Unprocessable Content"
    status: ClassVar[int] = 422

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class NotInRoom(ProblemType):
    """The operation requires being in a room the user hasn't joined.

    This error occurs when attempting to perform an operation
    that requires the user to be in a specific Socket.IO room,
    but the user is not currently in that room.
    """

    title: ClassVar[str] = "Precondition Failed"
    status: ClassVar[int] = 412


class FrameNotFound(ProblemType):
    """The requested frame does not exist.

    This error occurs when attempting to access a frame by index
    that is out of bounds or does not exist in the room's frame storage.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise IndexError(problem.detail or "Frame not found")


class GeometryNotFound(ProblemType):
    """The requested geometry does not exist.

    This error occurs when attempting to access a geometry by key
    that does not exist in the room's geometry storage.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class BookmarkNotFound(ProblemType):
    """The requested bookmark does not exist.

    This error occurs when attempting to access a bookmark at a frame
    index that has no bookmark set.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class SelectionGroupNotFound(ProblemType):
    """The requested selection group does not exist.

    This error occurs when attempting to access a selection group by name
    that does not exist in the room.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class FigureNotFound(ProblemType):
    """The requested figure does not exist.

    This error occurs when attempting to access a figure by key
    that does not exist in the room.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class RoomLocked(ProblemType):
    """The room is locked and cannot be modified.

    This error occurs when attempting to modify a room that is either
    admin-locked or has an active edit lock held by another user.
    """

    title: ClassVar[str] = "Locked"
    status: ClassVar[int] = 423

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        from zndraw.client import RoomLockedError

        raise RoomLockedError(problem.detail or problem.title)


class MessageNotFound(ProblemType):
    """The requested message does not exist.

    This error occurs when trying to access or edit a message
    that doesn't exist or has been deleted.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class ProgressNotFound(ProblemType):
    """The requested progress tracker does not exist.

    This error occurs when attempting to update or complete a
    progress tracker that does not exist or has already completed.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class NotMessageOwner(ProblemType):
    """Only the message owner can edit this message.

    This error occurs when a user attempts to edit a message
    that was posted by another user.
    """

    title: ClassVar[str] = "Forbidden"
    status: ClassVar[int] = 403

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class Forbidden(ProblemType):
    """Admin privileges are required for this operation.

    This error occurs when a non-admin user attempts to perform
    an operation that requires admin privileges, such as user
    management, public job registration, or server shutdown.
    """

    title: ClassVar[str] = "Forbidden"
    status: ClassVar[int] = 403

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class ScreenshotNotFound(ProblemType):
    """The requested screenshot does not exist.

    This error occurs when attempting to access a screenshot by ID
    that does not exist in the room.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)


class ScreenshotTooLarge(ProblemType):
    """The uploaded screenshot exceeds the maximum allowed size.

    This error occurs when the file size exceeds 10 MB.
    """

    title: ClassVar[str] = "Content Too Large"
    status: ClassVar[int] = 413

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class InvalidScreenshotFormat(ProblemType):
    """The uploaded screenshot has an unsupported format.

    Supported formats are: png, jpeg, webp.
    """

    title: ClassVar[str] = "Unprocessable Content"
    status: ClassVar[int] = 422

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class NoFrontendSession(ProblemType):
    """The specified session is not an active frontend session.

    This error occurs when requesting a screenshot capture from a
    session that is not currently connected.
    """

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class ScreenshotNotPending(ProblemType):
    """The screenshot is not in pending status.

    This error occurs when attempting to complete a screenshot
    that has already been completed.
    """

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


class RoomReadOnly(ProblemType):
    """The room has a mounted source and cannot be modified.

    This occurs when attempting to write (append, delete, set) to a room
    that has a provider-backed frame source mounted.
    """

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise PermissionError(problem.detail or problem.title)


class StepOutOfBounds(ProblemType):
    """The step value is outside the valid frame range.

    This error occurs when attempting to set the step to a value
    that is negative or exceeds the total number of frames.
    """

    title: ClassVar[str] = "Unprocessable Content"
    status: ClassVar[int] = 422

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ValueError(problem.detail or problem.title)


# Registry of all problem types
PROBLEM_TYPES: dict[str, type[ProblemType]] = {
    cls.problem_id(): cls
    for cls in [
        InvalidCredentials,
        UsernameExists,
        NotAuthenticated,
        UserNotFound,
        RoomNotFound,
        SessionNotFound,
        NotRoomMember,
        AlreadyRoomMember,
        InvalidPayload,
        UnprocessableContent,
        NotInRoom,
        FrameNotFound,
        GeometryNotFound,
        BookmarkNotFound,
        SelectionGroupNotFound,
        FigureNotFound,
        RoomLocked,
        MessageNotFound,
        ProgressNotFound,
        NotMessageOwner,
        Forbidden,
        StepOutOfBounds,
        ScreenshotNotFound,
        ScreenshotTooLarge,
        InvalidScreenshotFormat,
        NoFrontendSession,
        ScreenshotNotPending,
        RoomReadOnly,
    ]
}


def register_problem_types(*types: type[ProblemType]) -> None:
    """Register additional problem types (e.g., from worker system)."""
    for cls in types:
        PROBLEM_TYPES[cls.problem_id()] = cls


register_problem_types(ProviderTimeout)  # type: ignore[arg-type]


def problem_responses(
    *problem_types: type[ProblemType],
) -> dict[int | str, dict[str, Any]]:
    """Combine multiple ProblemType classes into a responses dict for OpenAPI.

    When multiple problem types share the same status code, their
    descriptions are joined with `` | `` so no information is lost.
    """
    result: dict[int | str, dict[str, Any]] = {}
    for pt in problem_types:
        new = pt.openapi_response()
        for code, entry in new.items():
            if code in result:
                existing_desc = result[code].get("description", "")
                new_desc = entry.get("description", "")
                result[code]["description"] = f"{existing_desc} | {new_desc}"
            else:
                result[code] = entry
    return result


# =============================================================================
# Pydantic Models and Exception Handler
# =============================================================================


class ProblemDetail(BaseModel):
    """RFC 9457 Problem Details."""

    type: str = "about:blank"
    title: str
    status: int
    detail: str | None = None
    instance: str | None = None


class ProblemException(Exception):
    """Exception that carries a ProblemDetail for RFC 9457 responses."""

    def __init__(
        self,
        problem: ProblemDetail,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.problem = problem
        self.headers = headers
        super().__init__(problem.title)


async def problem_exception_handler(request: Request, exc: Exception) -> JSONResponse:
    """Convert ProblemException to RFC 9457 compliant JSON response."""
    assert isinstance(exc, ProblemException)  # for type checker
    if exc.problem.instance is None:
        exc.problem.instance = request.url.path
    return JSONResponse(
        status_code=exc.problem.status,
        content=exc.problem.model_dump(exclude_none=True),
        media_type="application/problem+json",
        headers=exc.headers,
    )
