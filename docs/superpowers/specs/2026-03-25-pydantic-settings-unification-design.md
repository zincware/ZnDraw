# Pydantic-Settings Unification

**Date:** 2026-03-25
**Scope:** zndraw-fastapi (primary), zndraw-auth + zndraw-joblib (pyproject.toml support)
**Goal:** Unified, consistent settings resolution across server, client, and CLI using pydantic-settings custom sources.

## Problem

Settings resolution is inconsistent across entry points:

1. **Server-side** (`Settings`, `AuthSettings`, `JobLibSettings`) uses pydantic-settings with env vars only. No `pyproject.toml` support, so project-scoped server config requires env vars.

2. **Python client** (`ZnDraw` dataclass) ignores env vars entirely. `ZNDRAW_URL`, `ZNDRAW_USER`, `ZNDRAW_PASSWORD` only work via Typer's `envvar` in the CLI, not when constructing `ZnDraw()` in Python.

3. **CLI** (`zndraw-cli`) wires env vars through Typer's `envvar` parameter, completely separate from pydantic-settings. Duplicates resolution logic in `cli_agent/connection.py`.

4. **Auto-discovery** uses PID files (`~/.zndraw/server-{port}.pid`) for URL and a separate `tokens.json` for stored credentials. Two files, two formats, no unified state.

5. **`zndraw-cli auth login`** stores tokens in `~/.zndraw/tokens.json` but doesn't set a default server. Users must pass `--url` every time after login, despite the server being known.

6. **`shutdown_token`** is written to PID files and sent as `X-Shutdown-Token` header, but the server's `/v1/admin/shutdown` endpoint only checks `AdminUserDep` -- the token is never validated server-side.

7. **Local server admin** requires configuring `ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL` + `ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD` or relying on dev mode (all users are superusers). The user who started the server should be admin automatically.

## Design

### Approach: Composition

Keep `ZnDraw` as a `@dataclass` (it implements `MutableSequence`). Add a `ClientSettings(BaseSettings)` that resolves connection parameters through the full source chain. `ZnDraw.__post_init__` creates `ClientSettings` from its constructor args to fill in missing values.

Server-side classes each gain `PyprojectTomlConfigSettingsSource` independently.

### Two Phases

- **Phase 1 -- Server-side**: Add `pyproject.toml` support to `Settings`, `AuthSettings`, `JobLibSettings`.
- **Phase 2 -- Client/CLI-side**: New `ClientSettings`, `StateFile`, refactor `ZnDraw` + CLI connection logic, local admin token.

No shared code between phases. Phase 1 can ship independently.

## Phase 1: Server-Side pyproject.toml Support

### Changes

Each of the three settings classes adds `PyprojectTomlConfigSettingsSource` in `settings_customise_sources()`:

| Class | Table Header | Env Prefix |
|---|---|---|
| `Settings` | `[tool.zndraw.server]` | `ZNDRAW_` |
| `AuthSettings` | `[tool.zndraw.auth]` | `ZNDRAW_AUTH_` |
| `JobLibSettings` | `[tool.zndraw.joblib]` | `ZNDRAW_JOBLIB_` |

Note: `Settings` uses `[tool.zndraw.server]` (not `[tool.zndraw]`) to avoid field name clashes with `ClientSettings` (Phase 2), which reads from `[tool.zndraw]`. Env prefix stays `ZNDRAW_` for backward compatibility.

### Priority (all server-side classes)

```
init args (Settings(port=9000))              highest
    |
env vars (ZNDRAW_PORT=9000)
    |
pyproject.toml ([tool.zndraw.server])        lowest
```

### Example pyproject.toml

Each settings class reads from its own sub-table. No shared tables, no `extra="ignore"` needed.

```toml
[tool.zndraw]
# Client-side (read by ClientSettings, Phase 2)
url = "https://zndraw.icp.uni-stuttgart.de"
room = "my-project-room"
user = "researcher@uni-stuttgart.de"
password = "secret"  # user's responsibility if committed

[tool.zndraw.server]
# Server-side (read by Settings)
database_url = "postgresql+asyncpg://user:pass@db/zndraw"
redis_url = "redis://redis:6379"
storage = "/data/frames.lmdb"
host = "0.0.0.0"
port = 8000

[tool.zndraw.auth]
secret_key = "production-secret"
token_lifetime_seconds = 7200
default_admin_email = "admin@example.com"
default_admin_password = "secure-password"

[tool.zndraw.joblib]
worker_timeout_seconds = 120
```

### Implementation

Each class overrides `settings_customise_sources()`:

```python
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic_settings.main import PydanticBaseSettingsSource

class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        env_nested_delimiter="__",
        pyproject_toml_table_header=("tool", "zndraw", "server"),
    )

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        from pydantic_settings import PyprojectTomlConfigSettingsSource

        return (
            init_settings,
            env_settings,
            PyprojectTomlConfigSettingsSource(settings_cls),
        )
    # ... fields unchanged ...
```

Same pattern for `AuthSettings` (header: `("tool", "zndraw", "auth")`) and `JobLibSettings` (header: `("tool", "zndraw", "joblib")`).

## Phase 2: Client/CLI-Side Unification

### Component Overview

```
ClientSettings(BaseSettings)
    resolves: url, room, user, password, token
    sources: init > env > pyproject.toml > StateFileSource > guest fallback

StateFile
    replaces: ServerInfo + TokenStore
    file: ~/.zndraw/state.json (mode 0600)
    manages: unified server registry (local + remote), stored tokens

StateFileSource(PydanticBaseSettingsSource)
    custom pydantic-settings source
    reads StateFile for url (health-check-based discovery) and token (stored credentials)

Local Admin Token
    server writes local_token to state.json on start
    server accepts local_token as superuser auth
    ClientSettings uses local_token when connecting to local server
```

### ClientSettings

New `BaseSettings` subclass for client-side connection parameters:

```python
class ClientSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        pyproject_toml_table_header=("tool", "zndraw"),
    )

    url: str | None = None
    room: str | None = None
    user: str | None = None
    password: SecretStr | None = None
    token: str | None = None
```

`ClientSettings` reads from `[tool.zndraw]` (the top-level table). Server-side `Settings` reads from `[tool.zndraw.server]`. No table overlap, no `extra="ignore"` needed.

### Priority Chain

```mermaid
flowchart TD
    A["init args / CLI flags<br/><i>ZnDraw(url=...) or --url</i>"] --> B["env vars<br/><i>ZNDRAW_URL, ZNDRAW_TOKEN, ...</i>"]
    B --> C["pyproject.toml<br/><i>[tool.zndraw] url = ...</i>"]
    C --> D["StateFileSource<br/><i>~/.zndraw/state.json</i>"]
    D --> E{{"Guest fallback<br/><i>POST /v1/auth/guest</i>"}}

    style A fill:#2d5016,stroke:#4a8c2a,color:#fff
    style B fill:#1a4a6b,stroke:#2980b9,color:#fff
    style C fill:#4a3560,stroke:#7d3c98,color:#fff
    style D fill:#6b3a1a,stroke:#d35400,color:#fff
    style E fill:#4a4a4a,stroke:#888,color:#fff
```

Higher sources override lower sources. Guest fallback is not a settings source -- it runs in `ClientSettings` post-validation when `token` is still `None`.

### StateFile

Replaces both `ServerInfo` + `TokenStore` with a single `~/.zndraw/state.json` (mode 0600):

```json
{
  "servers": {
    "http://localhost:8000": {
      "pid": 12345,
      "version": "0.5.0",
      "local_token": "random-per-start",
      "added_at": "2026-03-25T10:00:00Z"
    },
    "https://zndraw.icp.uni-stuttgart.de": {
      "added_at": "2026-03-24T09:00:00Z"
    }
  },
  "tokens": {
    "https://zndraw.icp.uni-stuttgart.de": {
      "access_token": "eyJ...",
      "email": "user@uni.de",
      "stored_at": "2026-03-25T10:00:00Z"
    }
  }
}
```

Both `servers` and `tokens` use the full URL as key. Local servers have `pid`, `version`, and `local_token`. Remote servers only have `added_at`. All entries have `added_at` for ordering.

Key changes from current state:
- **Merged**: PID files + tokens.json into one file.
- **Unified server registry**: Both local and remote servers live in `servers`. No separate `default_url` — server preference is determined by health checks and recency (see URL Resolution below).
- **`local_token`**: Replaces `shutdown_token`. Generated per server start. Grants superuser access to the local server. Only present on local server entries.
- **`added_at`**: Timestamp for ordering. Newest healthy server wins within each category (localhost first, then remote).
- **Atomic writes**: Write to temp file, then rename (prevents corruption on crash). Note: atomic rename prevents partial reads but does not prevent lost updates from concurrent read-modify-write cycles (e.g., `server start` and `auth login` running simultaneously). This is acceptable given the low frequency of these operations — both are user-initiated CLI commands unlikely to overlap. If needed later, `fcntl.flock` can be added.
- **Migration**: On first access, if old `server-*.pid` or `tokens.json` files exist, migrate them to `state.json` and delete the old files.

### StateFileSource: URL Resolution

The `StateFileSource` resolves `url` from `state.json` using health checks against `/v1/health`. This works for both local and remote servers uniformly.

```mermaid
flowchart TD
    A["Read state.json → servers"] --> B["Partition by localhost<br/>vs non-localhost"]
    B --> C["Sort each partition by<br/>added_at descending (newest first)"]
    C --> D["Try localhost servers first"]

    D --> E{"Has pid?"}
    E -->|yes| F{"os.kill(pid, 0)<br/>alive?"}
    E -->|no| G["GET /v1/health<br/>(timeout 2s)"]
    F -->|yes| G
    F -->|no| H["Skip (dead PID)"]
    G -->|200 OK| I["Use this server's URL"]
    G -->|fail| H

    H --> J{"More localhost<br/>servers?"}
    J -->|yes| D
    J -->|no| K["Try remote servers"]

    K --> L["GET /v1/health<br/>(timeout 2s)"]
    L -->|200 OK| I
    L -->|fail| M{"More remote<br/>servers?"}
    M -->|yes| K
    M -->|no| N["url stays None<br/>(higher source or error)"]

    style I fill:#2d5016,stroke:#4a8c2a,color:#fff
    style N fill:#6b1a1a,stroke:#c0392b,color:#fff
```

**How it works:**
1. **Partition** servers into localhost (URL contains `localhost` or `127.0.0.1`) and non-localhost.
2. **Sort** each partition by `added_at` descending — newest first.
3. **Try localhost first** — for entries with a `pid`, do a fast `os.kill(pid, 0)` pre-filter to skip dead processes without a network call. Then hit `/v1/health` to confirm.
4. **Then try remote** — hit `/v1/health` directly.
5. **First healthy server wins.**

**Cleanup during resolution:**
- **Local entries with dead PID**: removed immediately — `os.kill(pid, 0)` is definitive, there is no "temporarily dead PID" scenario. Safe to clean up on every access.
- **Remote entries with failed health check**: kept — a failed health check is ambiguous (network blip, server restarting). Remote entries are only removed explicitly by `auth logout --url ...`.

This eliminates the need for `default_url`. The `auth login --url remote` command adds the remote to `servers`. Starting a local server adds to `servers`. The resolution algorithm always picks the best available server: localhost preferred (active intent), then remote (passive preference), newest first within each category.

### StateFileSource: Token Resolution

```mermaid
flowchart TD
    A["url resolved<br/>(from any source)"] --> B{"url points to<br/>local server?"}
    B -->|yes| C["Use local_token<br/>from state.json → servers"]
    B -->|no| D{"Token stored for<br/>this url?"}
    D -->|yes| E["Use stored access_token<br/>from state.json → tokens"]
    D -->|no| F["token stays None<br/>(guest fallback later)"]

    C --> G["Return token"]
    E --> G
    F --> G

    style C fill:#2d5016,stroke:#4a8c2a,color:#fff
    style E fill:#4a3560,stroke:#7d3c98,color:#fff
    style F fill:#4a4a4a,stroke:#888,color:#fff
```

"Points to local server" means the resolved URL matches a localhost entry in `state.json → servers` that passed the health check. In that case, `local_token` is used directly — it grants superuser access, bypassing the normal auth chain.

**Implementation detail:** Token resolution depends on the already-resolved `url`. In `StateFileSource.__call__`, the source accesses `self.current_state["url"]` (populated by higher-priority sources or the source's own URL resolution) to determine which token to return. If `url` is still `None` in `current_state`, no token is returned.

### Local Admin Token

The server generates a `local_token` on every start and writes it to `state.json`. The server accepts this token as superuser authentication.

#### Server-Side

On startup (`cli.py`):
1. Generate `local_token = secrets.token_urlsafe(32)` (replaces `shutdown_token`).
2. Write to `state.json → servers → {url} → local_token`.
3. Store `local_token` in `app.state.local_token`.

Auth endpoint or middleware:
- Accept `Authorization: Bearer <local_token>` as superuser. When the token matches `app.state.local_token`, return a synthetic superuser identity.
- This replaces the unused `X-Shutdown-Token` mechanism.
- Normal JWT auth continues to work for all other tokens (guests, logged-in users, remote clients).

`shutdown_server()` (in `server_manager.py`) is updated to read `local_token` from `state.json` and send it as `Authorization: Bearer <local_token>` instead of the old `X-Shutdown-Token` header. This means local shutdown now uses the same auth mechanism as all other local admin operations.

#### Client-Side

When `StateFileSource` resolves `url` to a local server, it also sets `token` to that server's `local_token`. The client connects as superuser with zero configuration.

#### Security Properties

- Token regenerated every server start -- no persistence across restarts.
- `state.json` is mode 0600 -- only the owner can read.
- Only applies to local servers -- remote servers use normal auth.
- Remote clients connecting to the same local server use normal auth (they can't read the state file).

### Full Resolution Flow

End-to-end flow showing how `ClientSettings` resolves all fields:

```mermaid
sequenceDiagram
    participant User as User Code
    participant ZD as ZnDraw.__post_init__
    participant CS as ClientSettings
    participant Init as Init Source
    participant Env as Env Source
    participant Toml as Pyproject Source
    participant SF as StateFile Source
    participant Server as ZnDraw Server

    User->>ZD: ZnDraw() or ZnDraw(url=...)
    ZD->>CS: ClientSettings(url=..., token=..., ...)
    CS->>Init: Read constructor args
    CS->>Env: Read ZNDRAW_* env vars
    CS->>Toml: Read [tool.zndraw] from pyproject.toml
    CS->>SF: Read ~/.zndraw/state.json

    Note over SF: URL: healthy localhost (newest) > healthy remote (newest)
    Note over SF: Token: local_token (if localhost) > stored token

    CS->>CS: Merge (init > env > toml > state file)
    CS-->>ZD: Resolved settings (url, token may be None)

    alt url is None
        ZD->>ZD: Raise ConnectionError
    end

    alt token is None
        ZD->>Server: POST /v1/auth/guest
        Server-->>ZD: guest token
    end

    ZD->>ZD: Create APIManager, verify room

    alt token expired (on first API call)
        Server-->>ZD: 401
        ZD->>ZD: Delete stored token from state.json
        Note over ZD: Log warning: "Stored token expired,<br/>run 'zndraw-cli auth login'<br/>to re-authenticate"
        ZD->>Server: POST /v1/auth/guest
        Server-->>ZD: guest token
    end

    ZD-->>User: Ready ZnDraw instance
```

**Note:** `ClientSettings` only resolves field values from sources (no HTTP calls). Guest fallback, token validation, and expiry handling happen in `ZnDraw.__post_init__` after `ClientSettings` returns.

### ZnDraw Integration

`ZnDraw.__post_init__` creates `ClientSettings` from constructor args:

```python
@dataclass
class ZnDraw(MutableSequence[ase.Atoms]):
    url: str | None = None
    room: str | None = None
    user: str | None = None
    password: SecretStr | str | None = None
    token: str | None = None
    # ... other fields ...

    def __post_init__(self) -> None:
        # Normalize password
        if isinstance(self.password, str):
            self.password = SecretStr(self.password)

        # Resolve via pydantic-settings chain
        resolved = ClientSettings(
            url=self.url,
            room=self.room,
            user=self.user,
            password=self.password,
            token=self.token,
        )

        # url is required — if still None after full chain, error out
        if resolved.url is None:
            raise ConnectionError(
                "No ZnDraw server found. Pass url=, set ZNDRAW_URL, "
                "add [tool.zndraw] url to pyproject.toml, or start a local server."
            )
        self.url = resolved.url
        self.room = resolved.room or str(uuid.uuid4())

        # token: guest fallback if still None
        if resolved.token is None:
            self.token = self._guest_fallback(self.url)
        else:
            self.token = resolved.token

        # ... rest of init (APIManager, SocketManager, etc.) ...
```

**Validation boundaries:** `ClientSettings` itself does not validate that `url` or `token` are non-None — all fields are `Optional`. The "url is required" check and guest fallback live in `ZnDraw.__post_init__` (and equivalently in `cli_agent/connection.py`'s `get_connection`), because they require error handling specific to the caller (e.g., `ConnectionError` for Python client, `SystemExit` for CLI).

Constructor args that are not `None` pass through `init_settings` as highest priority. `None` values fall through to env > pyproject > state file.

### CLI Integration

`zndraw-cli` commands replace manual `resolve_url()` / `resolve_token()` with `ClientSettings`:

```python
# Before (cli_agent/connection.py)
def get_connection(url, token, user, password):
    base_url = resolve_url(url)          # manual PID file logic
    resolved_token = resolve_token(...)   # manual 3-tier fallback
    return Connection(base_url, resolved_token)

# After
def get_connection(url, token, user, password):
    settings = ClientSettings(url=url, token=token, user=user, password=password)
    return Connection(settings.url, settings.token)
```

### auth login / logout Flow

```mermaid
flowchart TD
    subgraph "zndraw-cli auth login --url URL"
        L1["Resolve URL"] --> L2["Device-code flow"]
        L2 --> L3["Store token in<br/>state.json → tokens[url]"]
        L3 --> L4["Add url to<br/>state.json → servers[url]<br/>with added_at = now"]
    end

    subgraph "zndraw-cli auth logout --url URL"
        O1["Resolve URL"] --> O2["Delete token from<br/>state.json → tokens[url]"]
        O2 --> O3["Remove url from<br/>state.json → servers[url]"]
    end

    subgraph "zndraw server start"
        S1["Generate local_token"] --> S2["Write to state.json → servers[url]<br/>with pid, local_token, added_at"]
        S2 --> S3["Store local_token in app.state"]
    end

    subgraph "zndraw server stop"
        T1["Remove state.json → servers[url]"]
    end
```

### Token Expiry

When a stored token (from `auth login`) expires:

1. `url` survives -- the user chose that server, the intent doesn't expire.
2. `token` is used optimistically -- server returns 401.
3. Expired token is deleted from `state.json → tokens`.
4. `log.warning("Stored token expired for %s. Falling back to guest. Run 'zndraw-cli auth login' to re-authenticate.", url)`.
5. Guest fallback provides a working (but unprivileged) session.

## Migration

### State File Migration

On first `StateFile` access:

1. Check for old `~/.zndraw/server-*.pid` files -- migrate to `state.json → servers`.
2. Check for old `~/.zndraw/tokens.json` -- migrate to `state.json → tokens`.
3. Delete old files after successful migration.
4. Old `shutdown_token` values are dropped (not migrated) since they were unused server-side.

### Backward Compatibility

- `ZnDraw(url=..., token=...)` continues to work exactly as before (init args highest priority).
- Env vars (`ZNDRAW_URL`, `ZNDRAW_TOKEN`, etc.) continue to work (env source).
- Old PID files and tokens.json are auto-migrated on first access.
- No breaking changes to the public API.

## Files Changed

### Phase 1
| File | Change |
|---|---|
| `src/zndraw/config.py` | Add `settings_customise_sources()` with `PyprojectTomlConfigSettingsSource` |
| `zndraw-auth: settings.py` | Add `settings_customise_sources()` with `PyprojectTomlConfigSettingsSource` |
| `zndraw-joblib: settings.py` | Add `settings_customise_sources()` with `PyprojectTomlConfigSettingsSource` |
| `tests/test_config.py` | Add tests for pyproject.toml source |

### Phase 2
| File | Change |
|---|---|
| `src/zndraw/client/settings.py` | **New**: `ClientSettings(BaseSettings)` |
| `src/zndraw/state_file.py` | **New**: `StateFile` class (replaces `ServerInfo` + `TokenStore`) |
| `src/zndraw/settings_sources.py` | **New**: `StateFileSource(PydanticBaseSettingsSource)` |
| `src/zndraw/client/core.py` | Refactor `__post_init__` to use `ClientSettings` |
| `src/zndraw/cli_agent/connection.py` | Replace `resolve_url`/`resolve_token` with `ClientSettings` |
| `src/zndraw/auth_utils.py` | Simplify -- guest fallback only (token resolution moved to settings chain) |
| `src/zndraw/server_manager.py` | Remove `ServerInfo`, `TokenStore`, PID file functions. Keep `wait_for_server_ready`, `is_process_running`, `shutdown_server` (updated to use `StateFile`) |
| `src/zndraw/cli.py` | Write `local_token` to `state.json` on server start. Remove entry on stop. |
| `src/zndraw/cli_agent/auth.py` | `auth login` adds remote to `servers`. `auth logout` removes it. |
| `src/zndraw/routes/admin.py` | Accept `local_token` as superuser auth (new dependency or middleware) |
| `src/zndraw/dependencies.py` | Add `LocalTokenOrAdminDep` for local admin token support |
| `tests/` | New tests for `ClientSettings`, `StateFile`, `StateFileSource`, local admin token, migration |

## Testing

### Phase 1
- `pyproject.toml` values loaded correctly for each settings class.
- Env vars override `pyproject.toml` values.
- Init args override env vars.
- Missing `pyproject.toml` is silently ignored (no error).

### Phase 2
- `ClientSettings` resolves through full chain (init > env > toml > state file > guest).
- `StateFile` read/write/migration from old format.
- Health-check-based discovery prefers localhost (newest) over remote (newest).
- `local_token` grants superuser access to local server.
- `local_token` is rejected for remote servers.
- `auth login` adds remote to `servers` and stores token.
- `auth logout` removes from `servers` and clears token.
- Token expiry falls back to guest with warning.
- Concurrent access to `state.json` (server + CLI) uses atomic writes.
- `ZnDraw()` with zero args connects to local server as admin.
- `ZnDraw(url="remote")` with stored token connects as authenticated user.
