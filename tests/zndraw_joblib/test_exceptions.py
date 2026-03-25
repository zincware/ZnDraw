# tests/test_exceptions.py
from zndraw_joblib.exceptions import (
    InvalidCategory,
    InvalidTaskTransition,
    JobNotFound,
    ProblemException,
    SchemaConflict,
    TaskNotFound,
    WorkerNotFound,
)


def test_problem_type_creates_problem_detail():
    problem = JobNotFound.create(detail="Job xyz not found")
    assert problem.type == "/v1/problems/job-not-found"
    assert problem.title == "Not Found"
    assert problem.status == 404
    assert problem.detail == "Job xyz not found"


def test_problem_type_creates_exception():
    exc = SchemaConflict.exception(detail="Schema mismatch")
    assert isinstance(exc, ProblemException)
    assert exc.problem.status == 409


def test_all_problem_types_have_correct_status():
    assert JobNotFound.status == 404
    assert SchemaConflict.status == 409
    assert InvalidCategory.status == 400
    assert WorkerNotFound.status == 404
    assert TaskNotFound.status == 404
    assert InvalidTaskTransition.status == 409


def test_internal_job_not_configured_problem():
    from zndraw_joblib.exceptions import InternalJobNotConfigured

    problem = InternalJobNotConfigured.create(
        detail="Job '@internal:modifiers:Rotate' not configured"
    )
    assert problem.status == 503
    assert problem.title == "Service Unavailable"
    assert problem.type == "/v1/problems/internal-job-not-configured"
    assert "@internal:modifiers:Rotate" in problem.detail
