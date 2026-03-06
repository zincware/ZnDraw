---
phase: 1
slug: client-package
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-05
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | pytest 9.0.2 |
| **Config file** | pyproject.toml (convention-based discovery) |
| **Quick run command** | `uv run pytest tests/test_client/ -x -q` |
| **Full suite command** | `uv run pytest tests/ -x` |
| **Estimated runtime** | ~300 seconds (full suite) |

---

## Sampling Rate

- **After every task commit:** Run `uv run pytest tests/test_client/ -x -q && uv run python -c "from zndraw import ZnDraw"`
- **After every plan wave:** Run `uv run pytest tests/ -x`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 300 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 01-01-01 | 01 | 1 | CLNT-01 | smoke | `uv run python -c "from zndraw.client import ZnDraw"` | N/A | ⬜ pending |
| 01-01-02 | 01 | 1 | CLNT-02 | unit | `uv run pytest tests/test_client/test_serialization.py -x` | ❌ W0 | ⬜ pending |
| 01-01-03 | 01 | 1 | CLNT-03 | unit | `uv run pytest tests/test_client/test_exceptions.py -x` | ❌ W0 | ⬜ pending |
| 01-01-04 | 01 | 1 | CLNT-04 | smoke | `uv run python -c "from zndraw.client.lock import ZnDrawLock"` | N/A | ⬜ pending |
| 01-01-05 | 01 | 1 | CLNT-05 | smoke | `uv run python -c "from zndraw.client.api import APIManager"` | N/A | ⬜ pending |
| 01-01-06 | 01 | 1 | CLNT-06 | smoke | `uv run python -c "from zndraw.client.socket import SocketManager"` | N/A | ⬜ pending |
| 01-01-07 | 01 | 1 | CLNT-07 | smoke | `uv run python -c "from zndraw.client.core import ZnDraw"` | N/A | ⬜ pending |
| 01-01-08 | 01 | 1 | CLNT-08 | regression | `uv run python -c "from zndraw import ZnDraw; print(ZnDraw)"` | N/A | ⬜ pending |
| 01-02-01 | 02 | 1 | CLNT-09 | regression | `uv run pytest tests/ -x` | ✅ | ⬜ pending |
| 01-02-02 | 02 | 1 | CLNT-10 | unit | `uv run pytest tests/test_client/ -x` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tests/test_client/__init__.py` — empty file for pytest discovery
- [ ] `tests/test_client/test_serialization.py` — covers CLNT-02, CLNT-10
- [ ] `tests/test_client/test_exceptions.py` — covers CLNT-03, CLNT-10

*Existing infrastructure covers framework and fixture requirements.*

---

## Manual-Only Verifications

*All phase behaviors have automated verification.*

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 300s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
