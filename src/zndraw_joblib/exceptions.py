# src/zndraw_joblib/exceptions.py
"""RFC 9457 Problem Details for HTTP APIs."""

import re
from typing import ClassVar, NoReturn

from fastapi import Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel


def _camel_to_kebab(name: str) -> str:
    """Convert CamelCase to kebab-case."""
    return re.sub(r"(?<!^)(?=[A-Z])", "-", name).lower()


class ProblemDetail(BaseModel):
    """RFC 9457 Problem Details."""

    MEDIA_TYPE: ClassVar[str] = "application/problem+json"

    type: str = "about:blank"
    title: str
    status: int
    detail: str | None = None
    instance: str | None = None


class ProblemError(Exception):
    """Exception that carries a ProblemDetail for RFC 9457 responses."""

    def __init__(
        self,
        problem: ProblemDetail,
        headers: dict[str, str] | None = None,
    ) -> None:
        self.problem = problem
        self.headers = headers
        super().__init__(problem.title)


class ProviderTimeoutError(Exception):
    """Raised on the client side when a provider does not respond in time."""


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
    def create(
        cls, detail: str | None = None, instance: str | None = None
    ) -> ProblemDetail:
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
    ) -> ProblemError:
        """Create a ProblemError from this type."""
        return ProblemError(cls.create(detail, instance), headers=headers)


async def problem_exception_handler(_request: Request, exc: Exception) -> JSONResponse:
    """Convert ProblemError to RFC 9457 compliant JSON response."""
    if not isinstance(exc, ProblemError):
        raise TypeError(f"Expected ProblemError, got {type(exc).__name__}")
    return JSONResponse(
        status_code=exc.problem.status,
        content=exc.problem.model_dump(exclude_none=True),
        media_type=ProblemDetail.MEDIA_TYPE,
        headers=exc.headers,
    )


# Problem Types


class JobNotFound(ProblemType):
    """The requested job does not exist."""

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404


class SchemaConflict(ProblemType):
    """Job schema differs from existing registration."""

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409


class InvalidCategory(ProblemType):
    """Job category is not in the allowed list."""

    title: ClassVar[str] = "Bad Request"
    status: ClassVar[int] = 400


class WorkerNotFound(ProblemType):
    """The requested worker does not exist."""

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        """Raise KeyError for client-side 404 handling."""
        raise KeyError(problem.detail or problem.title)


class TaskNotFound(ProblemType):
    """The requested task does not exist."""

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        """Raise KeyError for client-side 404 handling."""
        raise KeyError(problem.detail or problem.title)


class InvalidTaskTransition(ProblemType):
    """Invalid task status transition."""

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409


class InvalidRoomId(ProblemType):
    """Room ID contains invalid characters (@ or :)."""

    title: ClassVar[str] = "Bad Request"
    status: ClassVar[int] = 400


class Forbidden(ProblemType):
    """Admin privileges required for this operation."""

    title: ClassVar[str] = "Forbidden"
    status: ClassVar[int] = 403


class InternalJobNotConfigured(ProblemType):
    """Internal job is registered but no executor is available."""

    title: ClassVar[str] = "Service Unavailable"
    status: ClassVar[int] = 503


class NoWorkersAvailable(ProblemType):
    """Job has no connected workers to process the task."""

    title: ClassVar[str] = "Conflict"
    status: ClassVar[int] = 409


class ProviderNotFound(ProblemType):
    """The requested provider does not exist."""

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404


class ProviderTimeout(ProblemType):
    """The provider did not respond within the requested wait time."""

    title: ClassVar[str] = "Gateway Timeout"
    status: ClassVar[int] = 504

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise ProviderTimeoutError(problem.detail or problem.title)


class ProviderExecutionFailed(ProblemType):
    """The provider failed to execute the requested read."""

    title: ClassVar[str] = "Bad Request"
    status: ClassVar[int] = 400
