This is a new application and you must not consider migrations or backwards compatibility.
Design all new features with maintainability and performance in mind.
Always do a clean implementation, do not try to patch existing code for backwards compatibility, but remove and refactor it instead!
Use KISS, DRY, SOLID and YAGNI principles.
When refactoring, you can break backwards compatibility.
Always consider a better design approach compared to the existing one.
Consider multiple approaches, review them against the principles above, the existing methods and the overall architecture - and choose the best one.
When in doubt, ask for a review of your design approach before implementing it.
You MUST NEVER run `git add` / `git commit` / `git push`!
The frontend data always originates from the backend! No other data source has to be considered / edge cases beyond the backend data must not be handled!
You MUST NEVER `@pytest.mark.xfail` or similar - all tests must pass!
Always review, if the new code will be a workaround or a proper solution before implementing it.
DO NOT IMPLEMENT WORKAROUNDS!
Remove deprecated code instead of just marking it as deprecated!

### Handling Default Values

**Principle: Single Source of Truth**
All default values (e.g., `camera`, `particles`) must be defined exclusively within the **Pydantic model**. Do not scatter fallback logic throughout the codebase.

**Anti-Pattern (Hardcoded Fallbacks)**
Do not perform null checks combined with hardcoded literals.

```js
camera.near = sessionCameraData.near ?? 0.1;
```

**Best Practice (Schema-Driven)**
Rely entirely on the schema to populate default values during initialization. Ensure the data model is fully validated before usage.

```js
camera.near = sessionCameraData.near;
```


# Context7
Always use context7 when I need code generation, setup or configuration steps, or
library/API documentation. This means you should automatically use the Context7 MCP
tools to resolve library id and get library docs without me having to explicitly ask.

# Redis
You can not use LUA scripts!

# collections.abc
If sensible, implement collections.abc interfaces for your classes, such as MutableMapping or MutableSequence.

# Testing

To run tests use
- `uv run pytest tests/`

To run python code use
- `uv run python ...`

To format / lint code use exactly
- `uvx prek run --all-files`

To build the frontend use
- `uv sync --reinstall`

Tests can run for up to 15 minutes - be patient!

When designing new tests, read the old tests first to understand the existing patterns.
Use `pytest.mark.parametrize` to avoid code duplication.
Tests should be very specific and test only one thing.
Avoid complex test setups.
Each test must be a function, not a method of a class!

# Documentation
Use numpy style docstrings.
Docstrings must be concise and to the point.
Use type hints wherever possible. `import typing as t` if necessary, but use `list[int|float] | None` instead of `t.Optional[t.List[int|float]]`!

Do not create markdown files, unless specifically requested!
Do not create files you are not asked for. This includes files in /tmp

Imports should always be at the top of the file - unless they affect startup time of ZnDraw and can be lazy loaded.
