"""Problem-type registry presence tests."""

from zndraw.exceptions import PROBLEM_TYPES


def test_new_problem_types_registered() -> None:
    expected = {
        "group-not-found",
        "group-name-taken",
        "not-group-member",
        "not-group-admin",
        "last-group-admin",
        "group-has-rooms",
        "transfer-target-invalid",
        "share-link-not-found",
        "share-link-invalid",
    }
    missing = expected - set(PROBLEM_TYPES.keys())
    assert not missing, f"Missing problem types: {missing}"


def test_retired_problem_types_removed() -> None:
    retired = {"not-room-member", "already-room-member"}
    present = retired & set(PROBLEM_TYPES.keys())
    assert not present, f"Retired types still present: {present}"
