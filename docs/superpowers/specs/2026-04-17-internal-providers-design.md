# Internal Providers & Default Filesystem — Design Spec

**Date:** 2026-04-17
**Status:** Draft (companion to [`2026-04-17-dockview-ui-fixes-design.md`](2026-04-17-dockview-ui-fixes-design.md))

## Motivation

`ISSUES.md` line 17 asks for a default filesystem provider rooted at the directory where the server was started, configurable via `ZNDRAW_FILEBROWSER_PATH="."` with a `"NONE"` sentinel.

Today, providers are hosted exclusively by external Python clients: the client calls `vis.register_provider(cls)`, holds an open Socket.IO connection, and responds in real time when the API dispatches a read. Extensions have two halves — remote (`vis.register_job(cls)`) and internal (`@internal:...`, executed by the taskiq worker from `src/zndraw/broker.py`). Providers only have the remote half, so "a provider that ships with the server" has no natural home today.

This spec adds the missing half: **`@internal` providers executed by the taskiq worker**. The default filesystem provider is the first consumer.

## Goals

- Symmetry with extensions: any `Provider` subclass can be registered at `@internal`, discovered on server startup, executed by the taskiq worker, and appears automatically in every room's provider list (because the existing provider list query already includes `@internal`-ish scopes by name match; widens to include `@internal` explicitly — see "Backend changes" below).
- Default local filesystem provider available from first load, configured with one env var.
- `docker-compose`-safe with `zndraw` at `replicas: 3`: the API replicas don't host anything; the existing `taskiq-worker` service (single replica by default) is the host. Scaling the taskiq worker to >1 hits the same broadcast issue that remote providers hit — out of scope, same problem as before.

## Non-Goals

- Migrating existing remote providers (user-registered `vis.register_provider`) to the new pattern. They keep working exactly as today.
- A dedicated "filesystem worker" container. Unnecessary — `taskiq-worker` is the host.
- Per-room internal provider registration. Internal providers are always `@global`-visible.
- Reworking the `LoadFile` extension. It already runs via taskiq; unchanged.

## Architecture

**Today (remote provider path):**

```
API.read_provider()
  → Redis pub/sub emit → Socket.IO room providers:<full_name>
  → external Python client executes Provider.read(handler)
  → POST /v1/providers/{id}/results
  → API long-poll unblocks
```

**New (`@internal` provider path):**

```
API.read_provider()
  → detect provider.room_id == "@internal"
  → enqueue taskiq task: @internal:<category>:<Name>.read(params)
  → taskiq-worker picks up queue
  → executes Provider.read(handler) in-process
  → POST /v1/providers/{id}/results  (same endpoint)
  → API long-poll unblocks
```

Symmetric to the extension split — workers scale via taskiq, retries via taskiq, no open socket from the server to itself.

## Backend changes

### `zndraw_joblib`

- **New discovery helper** in `src/zndraw_joblib/registry.py` (or a sibling file, e.g. `provider_registry.py`): `register_internal_providers(broker, providers: list[type[Provider]], executor)`. Mirrors `register_internal_tasks` for extensions.
- Each `@internal` provider class becomes a taskiq task named `@internal:<category>:<ClassName>` with a handler that: (a) resolves the `handler` argument (for filesystem: `fsspec.filesystem("file")` rooted at the configured path), (b) calls `cls(**params).read(handler)`, (c) POSTs the result to `/v1/providers/{id}/results`.
- **New DB invariant:** on startup, ensure a `ProviderRecord` row exists for every `@internal` provider, analogous to `ensure_internal_jobs`. Idempotent; safe across multiple replicas (PK upsert on `(room_id, category, name)`).
- **Dispatch fork** in `read_provider` (`router.py:1146-1215`): if `provider.room_id == "@internal"`, use the taskiq broker instead of the Socket.IO emit. Same cache / long-poll / timeout handling.
- **Visibility widening:** `_room_provider_filter` (`router.py:970-974`) currently admits `@global` and the room's own providers. Widen to admit `@internal` too, so the existing `GET /v1/rooms/{room_id}/providers?category=filesystem` lists the default filesystem from every room. No separate "internal provider list" endpoint.

### `zndraw` (server)

- **New setting** in `src/zndraw/config.py`:
  ```python
  filebrowser_path: str = "."
  # "." = server cwd. Any non-"none" value roots the default filesystem provider at that path.
  # "none" (case-insensitive) disables the default filesystem entirely.
  ```
- **Provider discovery collector** parallel to `_collect_extensions` — returns the list of `@internal` `Provider` classes bundled with the server. For v1: just `FilesystemRead`.
- **`broker.py`:** after `register_internal_tasks`, additionally call `register_internal_providers(broker, _collect_providers(), internal_provider_executor)`. The executor resolves the filesystem handler from `settings.filebrowser_path` (unless `"none"`, in which case the provider is skipped entirely — no DB row, no task registration).

### `FilesystemRead` provider

- `src/zndraw/providers/filesystem.py` unchanged. It's already a `Provider` subclass with a `read(handler)` method. The `handler` it expects (`handler.glob`, `handler.info`, `handler.ls`) is exactly what `fsspec.filesystem("file")` provides.

## Frontend changes

The existing `GET /v1/rooms/{room_id}/providers?category=filesystem` query in `FilesystemPanel.tsx` will start returning a provider row once this spec lands. Additionally:

- **Hide the activity-bar icon when the filesystem providers list is empty.** Single change in `ActivityBar.tsx` (or at the registry level): query the providers endpoint once on app load; if the `filesystem` category returns zero providers for the current room, the `filesystem` icon is not rendered in any activity bar. Re-queries on `providers_invalidate` socket events (already fired by the backend when a provider is registered/unregistered).
- This matches the behavior for `filebrowser_path="none"` — the icon disappears entirely rather than showing a dead panel.

## Docker-compose guidance

No compose changes required.

- **Single `taskiq-worker` replica** (current default): hosts the filesystem provider.
- **Multiple taskiq-worker replicas:** hits the same Socket.IO broadcast problem that user-registered remote providers hit with multiple clients holding the same provider name. Existing limitation; this spec doesn't invent it.

If a user scales `taskiq-worker` to >1 in the future, the recommended workaround is the same as today: run a single dedicated taskiq-worker replica pinned to the filesystem-hosting role via a separate compose service.

## Settings summary

| `ZNDRAW_SERVER_FILEBROWSER_PATH` | Effect |
|---|---|
| `"."` (default) | `FilesystemRead` registered at `@internal:filesystem:FilesystemRead`, rooted at the taskiq-worker's cwd. |
| `"/abs/path"` | Same, rooted at `/abs/path`. Path resolved to absolute at startup. |
| `"none"` (case-insensitive) | No registration. No DB row. No taskiq task. Frontend's providers list comes back empty. Filesystem activity-bar icon hidden. |

## Migration & validation

- Typecheck + lint: zero errors.
- `uv run pytest`: existing tests pass; add:
  - Unit test: `register_internal_providers` creates the expected `ProviderRecord` rows.
  - Integration test: `GET /v1/rooms/{id}/providers?category=filesystem` returns the `@internal` provider by default.
  - Integration test: `ZNDRAW_SERVER_FILEBROWSER_PATH=none` → empty provider list.
  - Integration test: with the `@internal` provider registered, reading a directory via the existing `GET /v1/rooms/{id}/providers/{provider_name}` endpoint returns file listings (taskiq dispatch happens, result comes back).
- Manual: start `uv run zndraw`, open a room, open the filesystem panel — expect a listing rooted at cwd.

## Open questions

- **`@internal` provider naming convention:** `@internal:<category>:<ClassName>` mirrors extensions. Unambiguous; decided.
- **Should `@internal` providers be user-deletable?** Extensions aren't. Default: same for providers — `DELETE /v1/providers/{id}` checks user ownership; `@internal`-scoped providers have no owner, so the check refuses (same behavior as today's 403 for provider-owned-by-different-user path).
