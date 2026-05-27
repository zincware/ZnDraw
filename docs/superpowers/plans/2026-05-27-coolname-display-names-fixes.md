# Coolname Display Names — Post-Review Fixes

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Predecessor plan:** `docs/superpowers/plans/2026-05-27-coolname-display-names.md`
**PR:** [zincware/ZnDraw#932](https://github.com/zincware/ZnDraw/pull/932)
**Branch:** `worktree-coolname-display-names`

**Goal:** Resolve every issue raised against PR #932 — CI failures (32 tests across Python 3.11/3.12/3.13), local code-review findings (3 Critical / 7 Important / 8 Minor), and CodeRabbit inline + outside-diff comments — so the branch is mergeable.

**Architectural direction:** the display-name rollout itself stays. Fixes here are *corrective*: backfill what was missed, repair what regressed, and document the few deviations that were intentional.

---

## Source inventory

### CI failures (3 jobs, all Python versions)

| # | Test | Root cause |
|---|------|-----------|
| 1 | `test_broadcast_contract.py::test_every_route_action_emits_consistent_room_scoped_events` | `JobsInvalidate` emits surrogate room_id, not real `room.id` |
| 2-13 | `test_client_source.py::test_*` (12 tests) | Mix of: `copy_from` UUID parser crash, frame_count keys not cleaned for display-name addresses, frames never arrive |
| 14 | `test_internal_worker_sweeper.py::test_sweeper_deletes_remote_worker_but_preserves_internal_provider` | `User(...)` direct ctor omits new NOT NULL `display_name` |
| 15 | `test_joblib_channel.py::test_register_job_emits_to_surrogate_channel` | Joblib emits on display-name channel, expected surrogate UUID channel |
| 16 | `test_pyclient_frames_invalidate.py::test_pyclient_cached_length_updates_on_frame_append` | Test sends `owner_id`; API expects `owner` |
| 17 | `test_scope_e2e.py::test_full_group_workflow` | Schema regex rejects group-owner UUID; needs group-name string |
| 18 | `test_sessions.py::test_cross_user_sees_other_users_sessions` | Reads `SessionItem.email` (renamed to `display_name`) |
| 19-25 | `test_providers.py::test_*` (7 tests) | `User(...)` direct ctor omits `display_name` |
| 26-27 | `test_registry.py::test_*` (2 tests) | Same — direct `User(...)` ctor |
| 28 | `test_provider_dispatch.py::test_read_provider_releases_inflight_on_dispatch_failure` | Same — direct `User(...)` ctor |

### Local code-review findings (from `/superpowers:requesting-code-review`)

| Severity | Topic |
|----------|-------|
| Critical 1 | `routes/rooms.py:449` `copy_from` UUID parser crashes on display-name addresses |
| Critical 2 | Frontend `userLock` writers/readers disagree on shape (UUID vs email vs display_name) |
| Critical 3 | `@pytest.mark.protected` test modified without sign-off |
| Important 4 | `EditLockResponse` / `LockUpdate` still leak only `user_id` (UUID) — UI shows "Locked by `<uuid>`" |
| Important 5 | Joblib `Job.room_id` storage vs `list_jobs` lookup incoherent under new contract |
| Important 6 | `_load_room_by_segment` raises where callers expect `None` |
| Important 7 | (duplicate of #4 on the socket side) |
| Important 8 | Local-admin synthetic `User` omits `display_name` |
| Important 9 | `UserNotFound` raised for missing **groups** — wrong noun |
| Important 10 | Baseline-green CI artifact not attached |
| Minor | 8 minor items (regex duplication, `elif` chain, `local imports`, etc.) |

### CodeRabbit comments (CHILL profile, 2 reviews, 16 inline + 3 outside-diff)

Aggregated by file in the task tables below.

---

## Pre-flight (read once, before Task 1)

1. **Confirm branch.** `git rev-parse --abbrev-ref HEAD` → `worktree-coolname-display-names`.
2. **Baseline failure inventory.** Run the local equivalent of CI to reproduce all 32 failures:
   ```bash
   uv run pytest -x --no-header -q 2>&1 | tail -80
   ```
   If your local count differs from CI's 32 by more than ±2, **stop and report** before touching code (suggests an unrelated regression).
3. **Confirm `@pytest.mark.protected` matches** unchanged by this plan:
   ```bash
   grep -rn "@pytest.mark.protected" tests/ | sort
   ```
   Expected: `tests/zndraw/test_broadcast_contract.py`, `tests/zndraw/test_event_schema_audit.py`. The Phase E task explicitly addresses the broadcast-contract modification; do **not** touch event-schema-audit.
4. **Frontend baseline:** `cd frontend && bunx tsc --noEmit && bun run build`. Should be green (no Phase A-D fix is gated on this).

---

## Phase ordering

1. **Phase A — CI unblock (fixtures + obvious test bugs).** Lowest blast radius, highest test-coverage payoff. Targets 14 of 32 failures.
2. **Phase B — Production bugs that CI surfaced.** Fixes the remaining 18 CI failures.
3. **Phase C — EditLock display_name plumbing.** Resolves local-review Critical #2 + Important #4/#7. No CI signal today (no test exercises the broken path), but blocks user-visible features.
4. **Phase D — CodeRabbit cleanups + local-review minors.** Independent, parallelizable.
5. **Phase E — Protected-test resolution.** Decide: revert + fix root cause, or document deviation. Reaches user.

Phases A and B are sequential (B depends on green test infrastructure from A). C, D, E can run in parallel after B.

---

## Phase A — CI unblock

### Task A1: Backfill `display_name` in legacy `User(...)` test fixtures (11 of 32 failures)

**Files:**
- `tests/zndraw/test_internal_worker_sweeper.py` (lines ~45, ~53)
- `tests/zndraw_joblib/test_providers.py` (12 occurrences: 160, 206, 438, 676, 792, 1013, 1067, 1114, 1162, and 3 others)
- `tests/zndraw_joblib/test_registry.py` (lines 256, 290)
- `tests/zndraw_joblib/test_provider_dispatch.py` (User ctor)

**Symptom:**
```
sqlalchemy.exc.IntegrityError: (sqlite3.IntegrityError) NOT NULL constraint failed: user.display_name
```

**Root cause:** Tests call `User(id=..., email="...", hashed_password="x", ...)` directly, bypassing the `create_test_user_model` helper (`tests/zndraw/helpers.py:40`) which already auto-generates a display_name.

**Decision:** Add `display_name=...` to every direct `User(...)` ctor in failing tests. Do **not** rewrite them to use the helper — they're written this way deliberately to control superuser/active flags and IDs.

- [ ] **Step 1: Patch each direct `User(...)` call.** For each match found by:
  ```bash
  grep -rn "User(" tests/zndraw_joblib/ tests/zndraw/test_internal_worker_sweeper.py | grep -v "create_test_user"
  ```
  add a `display_name=...` kwarg. Use a deterministic, unique-per-test slug to avoid uniqueness collisions:
  ```python
  user = User(
      id=uuid.uuid4(),
      email="int@test",
      hashed_password="x",
      display_name=f"int-test-{uuid.uuid4().hex[:8]}",  # NEW
      is_active=True,
      ...
  )
  ```
  For tests creating multiple users in one transaction, use distinct slugs (`f"int-test-a-{...}"`, `f"int-test-b-{...}"`).

- [ ] **Step 2: Verify.** Run the affected suites:
  ```bash
  uv run pytest tests/zndraw_joblib/test_providers.py tests/zndraw_joblib/test_registry.py tests/zndraw_joblib/test_provider_dispatch.py tests/zndraw/test_internal_worker_sweeper.py -q
  ```
  All 11 tests should now pass.

**Success criteria:** 11 CI failures resolved, no new failures introduced.

---

### Task A2: Rename `SessionItem.email` → `display_name` in test

**File:** `tests/zndraw/test_sessions.py:247`

**Symptom:** `AttributeError: 'SessionItem' object has no attribute 'email'`

**Root cause:** Predecessor plan renamed `SessionItem.email` → `display_name` (`src/zndraw/schemas.py`), this test was missed.

- [ ] **Fix:** Change `items[0].email == "user1@local.test"` to:
  ```python
  assert items[0].display_name is not None  # auto-generated by registration
  ```
  (Or, if exact identity matters, add a `display_name=` to the test's user creation and assert against it.)

- [ ] **Verify:** `uv run pytest tests/zndraw/test_sessions.py::test_cross_user_sees_other_users_sessions -q`.

---

### Task A3: Fix `owner_id` → `owner` in pyclient test

**File:** `tests/zndraw/test_pyclient_frames_invalidate.py:60`

**Symptom:** `httpx.HTTPStatusError: Client error '422 Unprocessable Entity' for url '.../v1/rooms'`

**Root cause:** Test sends `{"owner_id": user_id}` where `user_id` is a UUID; API expects `{"owner": display_name}`.

- [ ] **Fix:**
  ```python
  # before
  "owner_id": user_id,
  # after
  "owner": user_display_name,  # pull from the user fixture
  ```
  If the test fixture doesn't already expose `display_name`, fetch it from the auth response or `create_test_user_in_db` return.

- [ ] **Verify:** `uv run pytest tests/zndraw/test_pyclient_frames_invalidate.py -q`.

---

### Task A4: Fix undefined `display_name` in `test_socketio_rooms.py`

**File:** `tests/zndraw/test_socketio_rooms.py:763-764` (also referenced by CodeRabbit inline #15)

**Symptom:** `NameError: display_name is not defined` in `RoomJoin(owner=display_name, ...)`.

**Affected tests:** `test_rest_rejects_updating_other_users_session_camera`, `test_same_room_frame_append_updates_sidebar`.

- [ ] **Fix:** Define `display_name` locally per test from the user fixture before constructing `RoomJoin`:
  ```python
  display_name = user.display_name  # or from create_test_user_in_db return
  room_join = RoomJoin(owner=display_name, room_name=room_name, client_type="frontend")
  ```

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_socketio_rooms.py::test_rest_rejects_updating_other_users_session_camera tests/zndraw/test_socketio_rooms.py::test_same_room_frame_append_updates_sidebar -q
  ```

---

### Task A5: Add status assertion before `.json()` in `test_auth_endpoints.py`

**File:** `tests/zndraw/test_auth_endpoints.py:239-240` (CodeRabbit inline #14)

**Why:** Currently the test reads `resp.json()` without asserting `resp.status_code == 200`. A 500 / 422 surfaces as a misleading payload assertion.

- [ ] **Fix:**
  ```python
  resp = await client.get("/v1/users/available-display-name")
  assert resp.status_code == 200, resp.text
  body = resp.json()
  ```

- [ ] **Verify:** No test failure expected to flip; this is a quality fix that improves diagnostics for *future* breaks.

---

### Task A6: Clean up `dependency_overrides` in `test_router_task_submit.py`

**File:** `tests/zndraw_joblib/test_router_task_submit.py:88-110` (CodeRabbit 2nd review, only inline)

**Why:** The test sets `app.dependency_overrides[resolve_dispatch_room_address] = _override` but never restores it. Causes flaky cross-test leakage when later tests in the same session use the same `app` fixture.

- [ ] **Fix:** Wrap in try/finally:
  ```python
  app.dependency_overrides[resolve_dispatch_room_address] = _override
  try:
      # ... existing assertions
  finally:
      app.dependency_overrides.pop(resolve_dispatch_room_address, None)
  ```

- [ ] **Verify:** Run the test, then any nearby test that uses `app`, both in isolation and together.

---

## Phase B — Production bugs (18 of 32 CI failures)

### Task B1: Fix `copy_from` UUID parser crash

**File:** `src/zndraw/routes/rooms.py:449-453`
**Sources:** Local-review Critical #1; CodeRabbit outside-diff comment; CI failures `test_copy_from_mounted_room_raises` (500≠409), `test_copy_from_room_copies_bookmarks` (500≠201).

**Symptom:** Frontend sends composed addresses (`happy-blue-rabbit/my-room`) as `copy_from`; `UUID(owner_str)` raises `ValueError`, surfaces as unhandled 500.

- [ ] **Fix:**
  ```python
  # before
  if "/" in copy_from:
      owner_str, _, name_part = copy_from.partition("/")
      source_room = await _load_room_by_address(
          session, UUID(owner_str), name_part
      )
  # after
  if "/" in copy_from:
      owner_str, _, name_part = copy_from.partition("/")
      source_room = await _load_room_by_segment(session, owner_str, name_part)
  else:
      source_room = await session.get(Room, copy_from)
  ```
  Add `_load_room_by_segment` to the existing import block at the top of `routes/rooms.py`.

- [ ] **Audit:** Ensure `_load_room_by_segment`'s exception behavior matches what `create_room`'s outer try/except expects. (Depends on Task D7 outcome; if `_load_room_by_segment` still raises, the wrapper catches `ProblemError` already.)

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_client_source.py::test_copy_from_mounted_room_raises tests/zndraw/test_client_source.py::test_copy_from_room_copies_bookmarks -q
  ```

---

### Task B2: Fix cleanup loop UUID-only parsing in `database.py`

**File:** `src/zndraw/database.py:431-441` (CodeRabbit outside-diff)
**CI failure:** `test_client_source.py::test_provider_disconnect_clears_frame_count` (expected 0 frames after disconnect, got 5).

**Symptom:** Provider-disconnect cleanup loop currently does `UUID(owner_part)` on each room's composed address; every new `display-name/room-name` address is silently dropped, leaving frame-count keys + `FramesInvalidate` broadcasts stale.

- [ ] **Fix:** Remove the UUID-only parsing branch and pass `owner_part` straight to a segment-aware lookup:
  ```python
  # before
  try:
      owner_uuid = UUID(owner_part)
  except ValueError:
      continue
  room = await _load_room_by_address(session, owner_uuid, name_part)
  # after
  from zndraw.dependencies import _load_room_by_segment
  room = await _load_room_by_segment(session, owner_part, name_part)
  if room is None:
      continue
  ```
  (Assumes Task D7 made `_load_room_by_segment` return `None` on miss. If executed before D7, keep a `try/except UserNotFound` around it.)

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_client_source.py::test_provider_disconnect_clears_frame_count -q
  ```

---

### Task B3: Fix `JobsInvalidate` emitting surrogate room_id

**Files:**
- `src/zndraw_joblib/router.py:373` (and `sweeper.py:84, 155`)
- Likely the resolution: pass the **canonical** `Room.id` UUID into `build_room_scoped_emission`, not the joblib surrogate.

**Sources:** Local-review Critical #3 (root cause behind protected-test failure); CI failure `test_broadcast_contract.py::test_every_route_action_emits_consistent_room_scoped_events` ("jobs_invalidate: room_id=UUID(...) != UUID(...)").
**Also:** CI failure `test_joblib_channel.py::test_register_job_emits_to_surrogate_channel` ("expected channel `room:0a17...`, got `room:test-user-ffb18e/test-room`") — note this test expects the **opposite** direction (surrogate UUID channel). The two tests pin different invariants.

**Investigation needed first:**

- [ ] **Step 1: Read both tests fully** to understand which channel the codebase is *supposed* to emit on after the refactor. Run:
  ```bash
  grep -nE "channel|room_id|surrogate" tests/zndraw/test_broadcast_contract.py tests/zndraw/test_joblib_channel.py
  ```
  Reconcile: the broadcast-contract test expects `data["room_id"] == room.id` (canonical UUID). The joblib-channel test expects emission on `room:{surrogate_uuid}` channel.
  These are NOT contradictory if: the **room_id field inside the payload** is the canonical UUID, but the **socket.io channel name** is the surrogate. Verify by reading `build_room_scoped_emission` and `emit`.

- [ ] **Step 2: Locate where `room_id` gets stamped into `JobsInvalidate`.** The router currently passes `room_id` (which is the WritableRoomDep surrogate) into `build_room_scoped_emission`. The fix is to **resolve the canonical room.id from the surrogate** before building the emission payload, while keeping the channel name pegged to the surrogate.

- [ ] **Step 3: Fix.** Likely shape (verify against current code first):
  ```python
  # router.py:373 area
  canonical_room_id = await _canonical_room_id_for_surrogate(session, room_id)
  emission = await build_room_scoped_emission(
      session, JobsInvalidate, canonical_room_id, channel_override=f"room:{room_id}"
  )
  ```
  If `build_room_scoped_emission` doesn't support `channel_override`, the fix may need to extend its signature. **Discuss with maintainer if signature change exceeds 10 LOC.**

- [ ] **Verify both directions:**
  ```bash
  uv run pytest tests/zndraw/test_broadcast_contract.py tests/zndraw/test_joblib_channel.py -q
  ```

---

### Task B4: Fix `get_writable_room_id` partitioning plain surrogate IDs

**File:** `src/zndraw/dependencies.py:497-499` (CodeRabbit inline #6 part 2)

**Symptom:** Function now always does `room_id.partition("/")` regardless of whether `"/"` is in the string. For a plain surrogate UUID, `name_part` becomes `""` and `_load_room_by_segment` is called with empty room name.

- [ ] **Fix:**
  ```python
  if "/" in room_id:
      owner_part, _, name_part = room_id.partition("/")
      room = await _load_room_by_segment(session, owner_part, name_part)
      if room is None:
          raise RoomNotFound.exception(room_id)
      return room.id
  # plain surrogate UUID — return as-is (or validate)
  return UUID(room_id)
  ```

- [ ] **Verify:** Run any joblib test that exercises `get_writable_room_id` with both forms — `test_joblib_channel`, `test_router_task_submit`, and broader joblib router tests.

---

### Task B5: Broaden `owner` Field regex to accept group names

**Files:**
- `src/zndraw/schemas.py:53` — `RoomCreate.owner`, `RoomPatchRequest.new_owner`, `RoomResponse.owner` (CodeRabbit inline #12)
- `src/zndraw/socket_events.py:64` — `RoomJoin.owner`, `RoomLeave.owner`, `TypingStart/Stop.owner` (CodeRabbit inline #13)

**CI failure:** `test_scope_e2e.py::test_full_group_workflow` (422 == 201, group-owned room creation).

**Symptom:** Current regex `^[a-z][a-z0-9-]{2,63}$` rejects valid group names (which may contain uppercase, underscores, or `@` prefixes).

**Decision required first:** What's the canonical form for the `owner` field when it refers to a group?
- Option A: Group `slug` column (add if missing — requires migration).
- Option B: Same regex, applied uniformly; rename group to match the regex on create.
- Option C: Alternation regex accepting either user-display-name or group-name shape.

Recommend **Option C** (least invasive, no schema migration). Group names that don't match should fall back to a derived slug returned by `build_public_address` (covered by Task D8).

- [ ] **Step 1: Decide owner-field regex.** Propose:
  ```python
  OWNER_FIELD_PATTERN = r"^[A-Za-z0-9_][A-Za-z0-9_-]{2,63}$"
  ```
  Or remove the pattern from schemas + socket events entirely (length-only `min_length=3, max_length=64`) and rely on the route-level pattern (which already only validates user-owned routes; group-owned addresses use `build_public_address` output).

- [ ] **Step 2: Update both files.** Apply the chosen pattern.

- [ ] **Step 3: Promote to shared constant.** If keeping a pattern, add to `zndraw_auth.display_names`:
  ```python
  OWNER_FIELD_PATTERN_STR = r"^[A-Za-z0-9_][A-Za-z0-9_-]{2,63}$"
  ```
  Import from both `schemas.py` and `socket_events.py`. (Folds Task D11 partially in.)

- [ ] **Step 4: Rewrite `test_scope_e2e.py:60`.** Use the group's name instead of UUID:
  ```python
  # before
  json={"owner_id": gid, "name": "e2e-grp", "visibility": "group"},
  # after
  json={"owner": group.name, "name": "e2e-grp", "visibility": "group"},
  ```

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_scope_e2e.py::test_full_group_workflow tests/zndraw/test_socketio_rooms.py -q
  ```

---

### Task B6: Investigate remaining `test_client_source.py` failures

**Tests still failing after B1 and B2:**
- `test_fetch_round_trip`, `test_fetch_returns_correct_atom_counts`, `test_mount_read_unmount_lifecycle`, `test_fetched_frame_has_radii`, `test_fetched_frame_has_colors`, `test_fetched_frame_has_connectivity`, `test_fetched_frame_preserves_existing_radii`, `test_remount_serves_new_source_data` — all `TimeoutError: Frame 0 not available within 5.0s`.
- `test_slice_on_mounted_room`, `test_slice_middle_range_on_mounted_room` — `IndexError: Frame indices [...] out of range (0-9)`.

**Hypothesis:** Joblib provider dispatch dispatches `ZnDraw(room=...)` with a composed display-name address (see PR description: *"`zndraw_joblib/router.py:submit_task` passes surrogate UUID to `ZnDraw(room=...)` — blocks executor-dispatch tests; needs UUID → composed display-name resolution before kiq dispatch"*). The `a4f48d7e` commit said it was fixed, but CI says otherwise.

- [ ] **Step 1: Re-read** `src/zndraw_joblib/router.py::submit_task` and `src/zndraw_joblib/dependencies.py::resolve_dispatch_room_address` — the path the dispatched job uses to reach the room.

- [ ] **Step 2: Run one failing test with verbose logs:**
  ```bash
  uv run pytest tests/zndraw/test_client_source.py::test_fetch_round_trip -xvs --log-cli-level=DEBUG 2>&1 | tail -80
  ```
  Look for: "room not found", 404s on frames lookup, joblib emission channel mismatch.

- [ ] **Step 3: Likely fix** — depending on what Step 2 shows, this collapses to one of:
  - Joblib's provider broadcasts `FramesInvalidate` on the wrong channel (related to B3).
  - `ZnDraw(room=...)` in the test gets a surrogate UUID instead of a composed address, and frames endpoints reject it.
  - The frame-count cleanup from B2 was only partial.

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_client_source.py -q
  ```

---

### Task B7: Fix edit-lock refresh persists original `msg`

**File:** `src/zndraw/routes/edit_lock.py:103-116` (CodeRabbit outside-diff)

**Why:** On lock refresh, the response uses `request.msg if request.msg is not None else holder.get("msg")`, but Redis is reset with the *original* `raw` blob. Subsequent `GET /edit-lock` reads the stale message.

- [ ] **Fix:** Update the in-memory holder before persisting:
  ```python
  holder = json.loads(raw)
  if holder["lock_token"] != lock_token:
      raise RoomLocked.exception("Room is being edited by another session")
  if request.msg is not None:
      holder["msg"] = request.msg
  await redis.set(key, json.dumps(holder), ex=settings.edit_lock_ttl)
  ttl = await redis.ttl(key)
  return EditLockResponse(
      locked=True,
      lock_token=holder["lock_token"],
      user_id=holder["user_id"],
      sid=holder.get("sid"),
      msg=holder.get("msg"),
      acquired_at=holder["acquired_at"],
      ttl=max(ttl, 0),
  )
  ```

- [ ] **Add test:** New test in `tests/zndraw/test_routes_edit_lock.py` that acquires → refreshes with new `msg` → reads back via GET → asserts new value.

- [ ] **Verify:** `uv run pytest tests/zndraw/test_routes_edit_lock.py -q`.

---

## Phase C — EditLock display_name plumbing

> Resolves local-review Critical #2 + Important #4 + #7 together.

### Task C1: Add `display_name` to `EditLockResponse` and `LockUpdate`

**Files:**
- `src/zndraw/schemas.py` — `EditLockResponse`
- `src/zndraw/socket_events.py` — `LockUpdate`
- `src/zndraw/routes/edit_lock.py` — populate `display_name` in all response constructions (acquire, refresh, GET, release-response if any)
- `src/zndraw/socketio.py` — populate `display_name` in `LockUpdate` emit sites

- [ ] **Add field:**
  ```python
  class EditLockResponse(BaseModel):
      locked: bool
      lock_token: str | None
      user_id: str | None
      display_name: str | None  # NEW
      ...
  ```

- [ ] **Source `display_name`** from the lock holder. The Redis blob (`holder`) needs to persist it too — add `holder["display_name"] = current_user.display_name` on acquire, and read it everywhere the response is built.

- [ ] **Migration:** Old Redis blobs (already-held locks) will have no `display_name` key. Either:
  - Accept `None` and let the UI fall back to user_id (current behavior).
  - Lookup the user by `user_id` once on first refresh and backfill.

  **Recommended:** `holder.get("display_name")` everywhere — no migration. Old blobs expire on TTL anyway.

- [ ] **Verify:**
  ```bash
  uv run pytest tests/zndraw/test_routes_edit_lock.py -q
  ```

### Task C2: Standardize all frontend `userLock` writers on `display_name`

**Files:**
- `frontend/src/hooks/socketHandlers/connectionHandlers.ts:142` — uses `editLockResponse.user_id`
- `frontend/src/hooks/socketHandlers/roomHandlers.ts` (around `setUserLock`) — uses `LockUpdate.user_id`
- `frontend/src/stores/slices/lockSlice.ts:110-112` — uses `user?.email`
- (Pseudo-imports + types in same files)

- [ ] **Wait for C1 to ship the new field**, then:
  - `connectionHandlers.ts:142` → write `editLockResponse.display_name`
  - `roomHandlers.ts` → write `LockUpdate.display_name`
  - `lockSlice.ts:110-112` (`acquireLock`) → write `user?.display_name`
- [ ] **Rename** `LockSlice.setUserLock(email, ...)` parameter to `setUserLock(displayName, ...)` (`lockSlice.ts:16, 36`).
- [ ] **Verify:** Manual smoke (acquire lock, reload page, verify still in edit mode) + add a frontend unit test if a slice test harness exists.

### Task C3: Align `GeometryGrid` lock comparison

**File:** `frontend/src/components/GeometryGrid.tsx:74, 136`

- [ ] **Fix:** Replace `currentUserEmail` with `currentUserDisplayName` in both lock-holder comparisons. Import via `useStore((s) => selectUserDisplayName(s))` or equivalent.

- [ ] **Audit other readers:** Run
  ```bash
  rg -n "userLock|userEmail" frontend/src/
  ```
  Make sure every consumer compares to display_name now.

---

## Phase D — CodeRabbit cleanups + local-review minors

> All independent; can be picked up in any order or by separate subagents in parallel.

### Task D1: Validate fallback in `display_names.py:39-50`

**File:** `src/zndraw_auth/display_names.py`
**Source:** CodeRabbit inline #4

- [ ] **Fix:** Wrap the final fallback (`f"{coolname.generate_slug(3)}-{secrets.token_hex(2)}"`) in the same validation loop — check `RESERVED_DISPLAY_NAMES`, `DISPLAY_NAME_PATTERN`, and DB uniqueness. Loop a small number of times (e.g. 3) before raising.
- [ ] **Add test** to `tests/zndraw_auth/test_display_names.py` confirming the fallback never produces a reserved-word collision.

### Task D2: Constraint-specific error in `users.py:76-79`

**File:** `src/zndraw_auth/users.py`
**Source:** CodeRabbit inline #5

- [ ] **Fix:** Inspect `IntegrityError.orig` (or `str(exc)`) for the violated constraint name. Raise `UsernameExists` only on the `display_name` unique index; raise the appropriate distinct error (e.g., `EmailExists`) for email-uniqueness violations; re-raise otherwise.

### Task D3: Retry display_name allocation on race in `routes/auth.py`

**File:** `src/zndraw/routes/auth.py:51-59`
**Source:** CodeRabbit inline #9

- [ ] **Fix:** Wrap `user_manager.create(UserCreate(..., display_name=...))` in a small retry loop (3-5 attempts). On `IntegrityError` whose constraint matches `display_name`, regenerate via `generate_unique_display_name` and retry. Re-raise on other errors.

### Task D4: Raise on missing user in `routes/groups.py:415-421`

**File:** `src/zndraw/routes/groups.py`
**Source:** CodeRabbit inline #10

- [ ] **Fix:**
  ```python
  user = await session.get(User, user_id)
  if user is None:
      raise UserNotFound.exception(user_id)
  return GroupMemberResponse(..., display_name=user.display_name)
  ```

### Task D5: Preserve `owner_kind` in `_resolve_room_owner`

**File:** `src/zndraw/routes/rooms.py:290-291`
**Source:** CodeRabbit inline #11

- [ ] **Fix:**
  ```python
  if resolved is None:
      return "", owner_kind, ""
  ```
  Instead of hardcoded `"user"`.

### Task D6: Race-guard suggestion in `RegisterDialog.tsx`

**File:** `frontend/src/components/RegisterDialog.tsx:51-56`
**Source:** CodeRabbit inline #3

- [ ] **Fix:** Track an `isDirty` ref set by the input's `onChange`. In the `fetchSuggestion.then` handler, only `setDisplayName(value)` if `!isDirty.current` and component is still mounted (`isMounted` cleanup flag in the effect).

### Task D7: `get_owner_uuid_from_segment` returns `Optional[UUID]`

**File:** `src/zndraw/dependencies.py:239-258`
**Source:** CodeRabbit inline #6 part 1 + local-review Important #6

- [ ] **Fix:** Change return type to `Optional[UUID]`, return `None` on miss instead of raising `UserNotFound`. Update all callers:
  - `_load_room_by_segment` (`dependencies.py:588`)
  - `verify_room` (`dependencies.py:177`)
  - `_load_access_context` (`dependencies.py:596`)
  - `get_share_context_two_segment` (search for it)
  - Any new caller introduced in Task B1/B2

  Each caller must now interpret `None` and emit the right error (typically `RoomNotFound` since "no such owner" is indistinguishable from "no such room" at the user level).

- [ ] **Verify:** Run `tests/zndraw/test_owner_resolution.py` and any test that pinned the old raise behavior. Update those tests if they asserted the exception path was reached.

### Task D8: `Group.name` slug-safety in `build_public_address`

**File:** `src/zndraw/models.py:112-116`
**Source:** CodeRabbit inline #8

- [ ] **Decide:** Add `Group.slug` column (migration) or apply a `slugify` helper at read time.
- [ ] **Recommended (no-migration path):** Apply existing slug helper inline:
  ```python
  group = await session.get(Group, room.owner_group_id)
  if group is None:
      return f"{room.owner_group_id}/{room.room_name}"
  segment = slugify(group.name)  # lowercase, hyphens
  return f"{segment}/{room.room_name}"
  ```
  Verify the slugified output passes `DISPLAY_NAME_PATTERN` (or the broadened owner-field pattern from B5).

### Task D9: Local-admin synthetic User `display_name` default

**File:** `src/zndraw/dependencies.py:83-85`
**Source:** Local-review Important #8

- [ ] **Fix:** Add `display_name="local-admin"` to the synthesized User.

### Task D10: Unguarded dict access in `cli.py:316`

**File:** `src/zndraw/cli.py:316`
**Source:** Local-review Minor

- [ ] **Fix:** Wrap `resp.json()["display_name"]` in try/except → `typer.BadParameter(...)`.

### Task D11: Promote display-name regex to shared constant

**Files affected (verify with grep):**
```bash
rg -nE "\^\[a-z\]\[a-z0-9-\]\{2,63\}\$" src/ frontend/src/ tests/
```
**Source:** Local-review Minor (recommendation #3)

- [ ] **Fix:** Single source of truth — `src/zndraw_auth/display_names.py::DISPLAY_NAME_REGEX_STR`. Import everywhere it's used (or, for the owner-field broader pattern, see Task B5).
- [ ] On the frontend, define `DISPLAY_NAME_REGEX` once in `src/utils/auth.ts` (or `src/constants.ts`) and import.

### Task D12: Rename `UserNotFound` → `OwnerNotFound` for unknown groups

**File:** `src/zndraw/exceptions.py` + every raise site for the group-miss case.
**Source:** Local-review Important #9

- [ ] **Decide:** Rename, or just adjust the *message* of `UserNotFound` to "Owner '...' not found" when raised from `get_owner_uuid_from_segment`.
- [ ] **Recommended:** Introduce `OwnerNotFound` next to `UserNotFound`. Keep `UserNotFound` for `users/{user_id}` endpoints, `OwnerNotFound` for owner-segment lookups.

### Task D13: Minor — `build_public_address` `elif` chain + `session` type

**File:** `src/zndraw/models.py:97-117`
**Source:** Local-review Minor

- [ ] Convert sequential `if`s to `elif`. Annotate `session` parameter as `AsyncSession`.

### Task D14: Minor — lift repeated local imports in `socketio.py`

**File:** `src/zndraw/socketio.py:218-220` (and 2 other handler bodies)
**Source:** Local-review Minor

- [ ] Move the in-function `from zndraw.access import can_read` / `from zndraw.dependencies import _load_room_by_segment, ...` to module-level. Confirm no circular-import regression by running:
  ```bash
  uv run python -c "import zndraw.socketio"
  ```

### Task D15: Minor — RegisterDialog 409 message specialization

**File:** `frontend/src/utils/auth.ts:115-122`
**Source:** Local-review Minor

- [ ] When the 409 `type` URL ends in `/username-exists`, surface "That display name is taken — try another or click regenerate" instead of the generic "Registration failed".

### Task D16: Minor — admin user listing exposes `display_name`

**File:** `src/zndraw/routes/admin.py:32-39` + `frontend/src/components/AdminPanel.tsx:111-119`
**Source:** Local-review Minor

- [ ] Add `display_name` to `AdminUserResponse`; surface it in the admin panel. Email stays for distinct identity.

---

## Phase E — Protected-test resolution

### Task E1: Resolve `test_broadcast_contract.py` modification

**File:** `tests/zndraw/test_broadcast_contract.py:70-75, 113, 168-171`
**Source:** Local-review Critical #3 + pre-flight protocol violation

**Status after Phase B3:** If Task B3 correctly fixes `JobsInvalidate.room_id` to carry the canonical room UUID, the protected test should pass with its **original** assertion form.

- [ ] **Step 1: Confirm B3 landed and test passes** — `uv run pytest tests/zndraw/test_broadcast_contract.py::test_every_route_action_emits_consistent_room_scoped_events -q`.
- [ ] **Step 2: Revert** the post-Phase-B-unnecessary `build_public_address` substitutions in `test_broadcast_contract.py:70-75, 113, 168-171`. The original UUID-form URL construction and `data["room_address"] == room_address` assertion should both work because the canonical fields are restored.

  - If after B3 the test *still* fails with the original form, the situation has shifted: B3 may have only restored `room_id` correctness; the `room_address` field shape is genuinely different now. In that case **stop and ask the maintainer**: either (a) approve the test rewrite (and document the contract change), or (b) introduce a new protected test alongside the old one for the new contract.

- [ ] **Step 3: Document the decision** in the PR description either way.

---

## Phase F — Sanity & merge prep

### Task F1: Full test suite green

- [ ] **Backend:**
  ```bash
  uv run pytest -x -q
  ```
  Expect: 0 failures. If a test marked `@pytest.mark.protected` is still failing here, do **not** proceed — Phase E missed something.

- [ ] **Frontend:**
  ```bash
  cd frontend && bunx tsc --noEmit && bun run build && bun run test
  ```

### Task F2: Prek green

- [ ] `uvx prek run --all-files`. Address any auto-formatting / lint issues that surface (likely zero given the prior `prek` housekeeping commits).

### Task F3: Manual smoke

- [ ] Register flow: server suggests a display_name, regenerate works, 409 surfaces nicely.
- [ ] Acquire edit lock → reload page → confirm still in edit mode (validates C1+C2+C3).
- [ ] Open a chat room → confirm author labels render as display_names.
- [ ] Open the admin panel (if applicable) → verify display_name shown.
- [ ] Browser console clean across all of the above.

### Task F4: PR update

- [ ] Push commits.
- [ ] Update PR description with:
  - Resolution of each of the 32 CI tests (link to the commit that fixed it).
  - Resolution of each Critical / Important local-review finding.
  - Documented decision for Task E1 (protected-test treatment).
- [ ] Re-request review from the original reviewer / CodeRabbit.

---

## Resolution mapping (cross-reference)

| Issue source | Issue ref | Resolved by task(s) |
|---|---|---|
| CI 3.x | `test_broadcast_contract::test_every_route_action_emits_consistent_room_scoped_events` | B3 + E1 |
| CI 3.x | `test_client_source::test_copy_from_*` (2 tests) | B1 |
| CI 3.x | `test_client_source::test_provider_disconnect_clears_frame_count` | B2 |
| CI 3.x | `test_client_source` — 10 other failures | B6 (depends on B1/B2/B3) |
| CI 3.x | `test_internal_worker_sweeper::*` | A1 |
| CI 3.x | `test_joblib_channel::test_register_job_emits_to_surrogate_channel` | B3 |
| CI 3.x | `test_pyclient_frames_invalidate::*` | A3 |
| CI 3.x | `test_scope_e2e::test_full_group_workflow` | B5 |
| CI 3.x | `test_sessions::test_cross_user_sees_other_users_sessions` | A2 |
| CI 3.x | `test_provider_dispatch::*`, `test_providers::*` (7), `test_registry::*` (2) | A1 |
| CI 3.x | `test_socketio_rooms::test_rest_rejects_...`, `test_same_room_frame_append_...` | A4 |
| Local-review | Critical #1 (copy_from UUID) | B1 |
| Local-review | Critical #2 (userLock shape) | C1 + C2 + C3 |
| Local-review | Critical #3 (protected test) | E1 (root cause B3) |
| Local-review | Important #4 / #7 (lock display_name) | C1 |
| Local-review | Important #5 (joblib room_id) | B3 + B4 |
| Local-review | Important #6 (load_room_by_segment) | D7 |
| Local-review | Important #8 (local-admin User) | D9 |
| Local-review | Important #9 (UserNotFound for groups) | D12 |
| Local-review | Important #10 (baseline artifact) | F1 + F4 |
| Local-review | Minor (regex dup) | D11 |
| Local-review | Minor (build_public_address elif) | D13 |
| Local-review | Minor (socketio.py imports) | D14 |
| Local-review | Minor (admin panel display_name) | D16 |
| Local-review | Minor (cli.py:316 dict access) | D10 |
| Local-review | Minor (409 message) | D15 |
| CodeRabbit inline #1 (plan endpoint mix) | docs/plan only | Skip (cosmetic, plan already executed) |
| CodeRabbit inline #2 (plan git add) | docs/plan only | Skip (already committed) |
| CodeRabbit inline #3 (RegisterDialog race) | RegisterDialog.tsx | D6 |
| CodeRabbit inline #4 (display_names fallback) | display_names.py | D1 |
| CodeRabbit inline #5 (users.py IntegrityError) | users.py | D2 |
| CodeRabbit inline #6 part 1 (get_owner_uuid raise) | dependencies.py:239 | D7 |
| CodeRabbit inline #6 part 2 (get_writable_room_id partition) | dependencies.py:497 | B4 |
| CodeRabbit inline #7 (test_socketio_rooms NameError) | test_socketio_rooms.py | A4 |
| CodeRabbit inline #8 (Group.name slug) | models.py:112 | D8 |
| CodeRabbit inline #9 (auth.py race) | routes/auth.py | D3 |
| CodeRabbit inline #10 (groups.py missing user) | routes/groups.py | D4 |
| CodeRabbit inline #11 (owner_kind hardcode) | routes/rooms.py:290 | D5 |
| CodeRabbit inline #12 (RoomCreate.owner regex) | schemas.py:53 | B5 |
| CodeRabbit inline #13 (socket_events owner regex) | socket_events.py:64 | B5 |
| CodeRabbit inline #14 (test_auth_endpoints assert) | test_auth_endpoints.py | A5 |
| CodeRabbit inline #15 (router_task_submit cleanup) | test_router_task_submit.py | A6 |
| CodeRabbit outside #1 (database.py cleanup) | database.py:431 | B2 |
| CodeRabbit outside #2 (edit_lock refresh) | edit_lock.py:103 | B7 |
| CodeRabbit outside #3 (rooms.py copy_from) | rooms.py:449 | B1 |

---

## Effort estimate

| Phase | Tasks | Estimated effort |
|-------|-------|------------------|
| A | 6 | 1-2h (mostly mechanical) |
| B | 7 | 4-6h (B3 + B6 are the unknowns) |
| C | 3 | 2-3h |
| D | 16 | 4-6h (parallelizable) |
| E | 1 | 30min after B3 lands |
| F | 4 | 1h |
| **Total** | **37** | **~13-19h** (single-agent serial); ~7-10h with Phase D parallelized |

---

## Risk register

| Risk | Likelihood | Mitigation |
|------|------------|------------|
| Task B3 requires `build_room_scoped_emission` signature change | Medium | Confirm with maintainer before extending. Fallback: emit explicitly without the helper for `JobsInvalidate`. |
| Task B5 broader regex breaks an existing test that pinned the narrow form | Medium | Run full `tests/zndraw/test_schemas_room.py` and `tests/zndraw/test_socketio_rooms.py` after; adjust expected values. |
| Task C1 Redis migration of pre-existing `holder` blobs causes lock-holder regressions in deployment | Low | `holder.get("display_name")` returns `None` → UI falls back to user_id (current). |
| Task D7 changing raise-to-None breaks an unrelated caller | Medium | Audit with `rg "get_owner_uuid_from_segment"` and update each callsite. |
| Task E1 turns out Task B3 didn't fully restore the contract | Medium | The fallback (new protected test alongside the old) is acceptable — escalate to maintainer. |

---

## Notes for the implementing agent

- **Run pytest selectively.** Don't repeatedly run the whole suite between every step — it's slow. Run only the test files relevant to the task in progress; do the full sweep at Task F1.
- **Prefer subagent-driven-development for Phase D.** Tasks D1-D16 are independent; dispatching them in parallel as subagents is appropriate (each gets ~150 LOC of context).
- **Don't reflexively delete protected-test markers** even if they're inconvenient. Phase E exists specifically to handle this.
- **Commit per task** with a conventional-commit message tying back to the task number, e.g., `fix(B1): copy_from accepts display-name addresses (#932)`.
