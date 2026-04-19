"""Static audit: no long-polling handler (one that uses SessionMakerDep directly)
may transitively depend on a yield-based SessionDep. Because yielded FastAPI
deps are kept alive for the full request, holding a SessionDep while a handler
tries to open its own session deadlocks under SQLite serialization.

The "long-polling handler" marker is that ``get_session_maker`` appears as a
*direct* dependency of the route's handler function (not transitively via the
normal ``get_session`` → ``get_session_maker`` wrapper).  Such handlers must
never also have ``get_session`` anywhere in their dep tree, because that yield
session would stay open for the entire request and block any factory-session
writes.
"""

from __future__ import annotations

from fastapi.routing import APIRoute

from zndraw.app import app
from zndraw_auth.db import get_session, get_session_maker


def _flatten(dep) -> list:
    out = [dep]
    for sub in dep.dependencies:
        out.extend(_flatten(sub))
    return out


def test_no_long_polling_handler_holds_full_request_session():
    """Long-polling handlers (direct SessionMakerDep) must not hold a SessionDep.

    A route is a "long-polling handler" when ``get_session_maker`` is one of the
    *direct* dependencies of the handler function itself (first-level only).
    Such routes MUST NOT have ``get_session`` anywhere in their transitive dep
    tree, or SQLite's serialised write lock will deadlock.
    """
    offenders: list[tuple[str, list[str]]] = []
    for route in app.routes:
        if not isinstance(route, APIRoute):
            continue
        dep_tree = route.dependant  # codespell:ignore dependant

        # Handler-level direct deps only (first level, not recursive)
        handler_direct_fns = [
            d.call for d in dep_tree.dependencies if d.call is not None
        ]

        # Only audit routes where the handler itself takes SessionMakerDep
        if get_session_maker not in handler_direct_fns:
            continue

        # Recursively walk the full dep tree to find any SessionDep usage
        all_fns = [d.call for d in _flatten(dep_tree) if d.call is not None]
        if get_session in all_fns:
            offenders.append(
                (
                    f"{route.methods} {route.path}",
                    [getattr(fn, "__name__", str(fn)) for fn in all_fns],
                )
            )

    assert not offenders, (
        "Long-polling handlers (direct SessionMakerDep) also hold a full-request "
        "SessionDep — deadlocks under SQLite serialization:\n"
        + "\n".join(f"  {route}: {names}" for route, names in offenders)
    )
