# ZnDraw Development Guidelines

> **Audience**: LLM-driven engineering agents and human developers

ZnDraw is a visualization and editor for atomistic data.
It features a rich plugin-infrastructure and multi-user support.
It can run locally via CLI or can run in a production docker environment.


## Required Development Workflow

**CRITICAL**: Always run these commands in sequence before committing:

### Python
```bash
uv sync                              # Install dependencies
uvx pre-commit run --all-files       # Linting and Code formatting
uv run pytest                        # Run full test suite
```
### TypeScript
```bash
bun install                              # Install dependencies
bun vite build                           # Build the frontend
bun run lint                             # Linting and Code formatting
bun run test                             # Run full test suite
```

**All must pass** - this is enforced by CI.

## Repository Structure
| Path              | Purpose                                                |
|-------------------|--------------------------------------------------------|
| `zndraw/`         | Python App                                             |
| `zndraw_app/`     | Python CLI app                                         |
| `app/`            | Typescript App                                         |
| `tests/`          | Unit and integration tests for the IPSuite library     |
| `docs/`           | Documentation for the IPSuite library                  |


## Code Standards
### Python
- Python ≥ 3.10 with full type annotations
- Follow existing patterns and maintain consistency
- **Prioritize readable, understandable code** - clarity over cleverness
- Avoid obfuscated or confusing patterns even if they're shorter
### TypeScript
- Follow existing patterns and maintain consistency
- **Prioritize readable, understandable code** - clarity over cleverness
