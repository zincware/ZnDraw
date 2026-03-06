---
phase: 3
slug: socket-handler-modules
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-06
---

# Phase 3 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Playwright 1.58.x + tsc (TypeScript compiler) |
| **Config file** | `frontend/playwright.config.ts` |
| **Quick run command** | `cd frontend && npx tsc --noEmit` |
| **Full suite command** | `cd frontend && npx playwright test` |
| **Estimated runtime** | ~10s (tsc), ~120s (E2E) |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx tsc --noEmit`
- **After every plan wave:** Run `cd frontend && npx tsc --noEmit && npx vite build`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 03-01-01 | 01 | 1 | SOCK-08 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-01-02 | 01 | 1 | SOCK-08 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-01 | 02 | 1 | SOCK-01 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-02 | 02 | 1 | SOCK-02 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-03 | 02 | 1 | SOCK-03 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-04 | 02 | 1 | SOCK-04 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-05 | 02 | 1 | SOCK-05 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-06 | 02 | 1 | SOCK-06 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-02-07 | 02 | 1 | SOCK-06 | compilation | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-03-01 | 03 | 2 | SOCK-07 | compilation + manual | `cd frontend && npx tsc --noEmit` | N/A | ⬜ pending |
| 03-03-02 | 03 | 2 | SOCK-09 | e2e | `cd frontend && npx playwright test` | ✅ 12 specs | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements.

- TypeScript compilation (`tsc --noEmit`) verifies SOCK-01 through SOCK-08 (structural correctness + type safety)
- Existing Playwright E2E specs verify SOCK-09 (behavioral equivalence)
- No new test files needed

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Orchestrator ~150 lines | SOCK-07 | Line count metric | `wc -l frontend/src/hooks/useSocketManager.ts` — target ≤ 180 |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
