"""RFC 9457 Problem Type documentation endpoints."""

from fastapi import APIRouter
from fastapi.responses import PlainTextResponse

from zndraw.exceptions import PROBLEM_TYPES, ProblemDetail, ProblemException

router = APIRouter(prefix="/v1/problems", tags=["problems"])


@router.get("")
def list_problem_types() -> dict[str, str]:
    """List all registered problem types with their URIs."""
    return {
        problem_id: problem_type.type_uri()
        for problem_id, problem_type in PROBLEM_TYPES.items()
    }


@router.get("/{problem_id}", response_class=PlainTextResponse)
def get_problem_type(problem_id: str) -> str:
    """Get documentation for a specific problem type."""
    problem_type = PROBLEM_TYPES.get(problem_id)
    if problem_type is None:
        raise ProblemException(
            ProblemDetail(
                type="about:blank",
                title="Not Found",
                status=404,
                detail=f"Problem type '{problem_id}' not found",
            )
        )

    doc = problem_type.__doc__ or "No documentation available."
    return (
        f"# {problem_type.__name__}\n\n"
        f"**Status:** {problem_type.status}\n"
        f"**Title:** {problem_type.title}\n\n"
        f"{doc}"
    )
