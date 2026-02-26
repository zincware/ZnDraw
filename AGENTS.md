# Guidelines

Use KISS, DRY, SOLID and YAGNI principles.
Do clean implementations — no patches for backwards compatibility, remove and refactor instead.
The frontend data always originates from the backend — no other data source edge cases.
You MUST NEVER `@pytest.mark.xfail` or similar — all tests must pass!
DO NOT IMPLEMENT WORKAROUNDS! Remove deprecated code instead of marking it.
Redis is available for testing. You know redis, so don't do `... if isinstance(obj, bytes) else ...`.

### Default Values
All defaults must be defined in the **Pydantic model** (single source of truth).
No `camera.near = data.near ?? 0.1` — just `camera.near = data.near`.

# Context7
Always use context7 for code generation, setup, configuration, or library docs.

### REST API
- https://context7.com/microsoft/api-guidelines
- https://context7.com/zalando/restful-api-guidelines
- https://context7.com/websites/fastapi_tiangolo
Use RFC 9457 Problem JSON (application/problem+json) for all 4xx/5xx responses.

### Socket.io / WebSockets
- https://context7.com/socketio/socket.io
- https://context7.com/miguelgrinberg/python-socketio
Use Pydantic models for all socket responses.

### Workers
- https://context7.com/websites/taskiq-python_github_io
- https://context7.com/taskiq-python/taskiq-redis

### Redis
- https://context7.com/redis/redis-py
- https://context7.com/redis/redis-doc
No LUA scripts!

### SQLModel / Databases
- https://context7.com/websites/sqlmodel_tiangolo
- https://context7.com/fastapi/sqlmodel

### Testing
- https://context7.com/pytest-dev/pytest
- https://context7.com/pytest-dev/pytest-asyncio

### Utility Libraries
- https://context7.com/pydantic/pydantic

# Development
## Dependencies
- `uv add <package>` / `uv remove <package>` / `uv add --dev <package>`

## Commands
- `uv run pytest tests/` — run tests
- `uv run pyright .` — type checking
- `uv run ruff format .` — format
- `uv run ruff check --select I --fix .` — fix import order

# Coding Guidelines
1. **Think first** — state assumptions, surface tradeoffs, ask when unclear.
2. **Simplicity** — minimum code, no speculative features, no abstractions for one-time use.
3. **Surgical changes** — touch only what you must, match existing style.
4. **Goal-driven** — define success criteria, loop until verified.

# collections.abc
Implement collections.abc interfaces (MutableMapping, MutableSequence) where sensible.

# Testing
Tests can run up to 15 minutes — be patient!
Use `pytest.mark.parametrize`. Each test must be a function, not a class method.
Tests should be specific and test one thing. Avoid complex setups.

# Documentation
Numpy-style docstrings, concise. Use type hints everywhere.
Use `list[int|float] | None` not `t.Optional[t.List[int|float]]`.
Imports at top of file unless they affect startup time (lazy load).
