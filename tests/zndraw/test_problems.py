"""Tests for Problem documentation endpoints."""

import pytest
from httpx import AsyncClient

from zndraw.exceptions import ProblemDetail


@pytest.mark.asyncio
async def test_list_problem_types(client: AsyncClient) -> None:
    """Test listing all registered problem types."""
    response = await client.get("/v1/problems")
    assert response.status_code == 200

    data = response.json()
    assert "invalid-credentials" in data
    assert "username-exists" in data
    assert "not-authenticated" in data
    assert "user-not-found" in data

    # Verify URIs are correct
    assert data["invalid-credentials"] == "/v1/problems/invalid-credentials"


@pytest.mark.asyncio
async def test_get_problem_type_documentation(client: AsyncClient) -> None:
    """Test getting documentation for a specific problem type."""
    response = await client.get("/v1/problems/invalid-credentials")
    assert response.status_code == 200
    assert response.headers["content-type"] == "text/plain; charset=utf-8"

    content = response.text
    assert "# InvalidCredentials" in content
    assert "**Status:** 401" in content
    assert "**Title:** Unauthorized" in content
    assert "username or password is incorrect" in content


@pytest.mark.asyncio
async def test_get_problem_type_not_found(client: AsyncClient) -> None:
    """Test that unknown problem type returns 404 with Problem JSON."""
    response = await client.get("/v1/problems/unknown-problem")
    assert response.status_code == 404
    assert response.headers["content-type"] == "application/problem+json"

    problem = ProblemDetail.model_validate(response.json())
    assert problem.title == "Not Found"
    assert problem.status == 404
    assert problem.detail is not None
    assert "unknown-problem" in problem.detail


@pytest.mark.asyncio
async def test_unhandled_exception_returns_problem_json(client: AsyncClient) -> None:
    """InternalServerError is registered in problem types and documented."""
    response = await client.get("/v1/problems/internal-server-error")
    assert response.status_code == 200
    content = response.text
    assert "# InternalServerError" in content
    assert "**Status:** 500" in content
    assert "**Title:** Internal Server Error" in content
