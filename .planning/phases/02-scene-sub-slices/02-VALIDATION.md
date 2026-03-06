---
phase: 2
slug: scene-sub-slices
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-06
---

# Phase 2 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Playwright 1.58.0 + TypeScript compiler |
| **Config file** | `frontend/playwright.config.ts` |
| **Quick run command** | `cd frontend && bunx tsc --noEmit` |
| **Full suite command** | `cd frontend && bunx playwright test` |
| **Estimated runtime** | ~60 seconds (tsc), ~120 seconds (E2E) |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && bunx tsc --noEmit`
- **After every plan wave:** Run `cd frontend && bunx playwright test`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 120 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 02-01-01 | 01 | 1 | SCEN-01 | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) | pending |
| 02-01-02 | 01 | 1 | SCEN-02 | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) | pending |
| 02-01-03 | 01 | 1 | SCEN-03 | compile + E2E | `cd frontend && bunx playwright test e2e/editing.spec.ts` | editing.spec.ts | pending |
| 02-01-04 | 01 | 1 | SCEN-04 | compile + E2E | `cd frontend && bunx playwright test e2e/geometry-drawing.spec.ts` | geometry-drawing.spec.ts | pending |
| 02-01-05 | 01 | 1 | SCEN-05 | compile + E2E | `cd frontend && bunx playwright test e2e/camera-session.spec.ts` | camera-session.spec.ts | pending |
| 02-01-06 | 01 | 1 | SCEN-06 | compile | `cd frontend && bunx tsc --noEmit` | N/A (type check) | pending |
| 02-01-07 | 01 | 1 | SCEN-07 | compile + E2E | `cd frontend && bunx playwright test` | All 12 specs | pending |

*Status: pending / green / red / flaky*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements.

TypeScript compiler (`tsc --noEmit`) validates structural type composition (SCEN-01 through SCEN-07). Playwright E2E specs validate behavioral correctness (SCEN-03 through SCEN-05, SCEN-07). No new test files needed for this structural refactor.

---

## Manual-Only Verifications

All phase behaviors have automated verification.

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 120s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
