# Dockview Performance Audit — Design Spec

**Date:** 2026-04-17
**Status:** Draft (companion to [`2026-04-17-dockview-ui-fixes-design.md`](2026-04-17-dockview-ui-fixes-design.md))

## Motivation

`ISSUES.md` line 10: *"Feels slower (compared to main). In general the performance of the UI is not great!"*

**Framing:** the goal is not parity with `main`. The goal is a UI that is objectively fast — the audit profiles the current branch to find what's slow, and fixes follow. If the resulting app outperforms `main`, good; if it merely matches, the audit was not aggressive enough.

Preliminary empirical signals collected while browsing the running app on branch `spec/dockview-ui-redesign`:

- JS heap ~125 MB after landing on a room (dev mode, with HMR — production build is expected to be smaller but the relative delta vs main is what matters).
- 250 network resources resolved by the time the viewer renders (dev again, inflated, but a comparison to main in the same mode is the useful number).
- `frontend/src/panels/registry.tsx` imports every `Panel.component` synchronously: `GeometryPanel` (pulls in `@mui/x-data-grid`), `PlotsBrowserPanel`, `ChatPanel`, `FilesystemPanel`, `SecondaryPanel`, `RoomsPanel`. Every dockview-using page hauls these in before the user clicks a single icon.
- No `React.lazy` / `Suspense` boundary anywhere in the panels layer.
- The MUI `x-data-grid` dependency is in the `mui` manual-chunk in `vite.config.ts`, so it's co-bundled with the rest of MUI instead of lazy-loaded with the panel that uses it.
- `ActivityBar.tsx` adds `window` drag listeners without `{ passive: true }` (three instances of the component run simultaneously).

These are candidate culprits, not a diagnosis. This spec is an **audit plan**, not an implementation plan — the goal is to produce a ranked list of regressions backed by numbers, then to propose and land fixes in follow-up PRs.

## Goals

- Quantify where time and memory go on the current branch in a handful of concrete metrics — without any commit comparison.
- Identify the top N (≈5) contributors by order of magnitude across three surfaces: **initial load**, **chrome interactions** (open/close panels, drag, resize), and **R3F / WebGL viewer** (frame time, draw-call count, memory under common interactions like orbiting, scrubbing, splitter drag).
- Propose a fix for each and estimate the expected delta.
- Land fixes in bounded PRs, each with a measurement before/after.
- Leave a **CI-visible perf budget** in place so future follow-ups can't silently regress.
- Target: the app feels objectively fast. Beating `main` is the floor, not the ceiling.

## Non-Goals

- Fixing anything in this spec. This is the audit.
- Rewriting the panels layer or the dockview integration — the recommendations here feed into a future implementation spec.
- Comparing against `main`. The historical diff is not interesting; the absolute perf of the current branch is.

## Audit methodology (full C path)

### Step 1 — Metric capture (current branch only)

Run both **dev** (`bun run dev` + `uv run zndraw`) and **prod** (`bun run build` + `uv run zndraw` serving the static bundle) modes on `HEAD`. Prod numbers are what ships; dev numbers are what the developer feels.

**Load metrics**, captured via Playwright automation (reusable harness saved to `frontend/tests/perf/`):

| Metric | Source |
|---|---|
| TTFB | `PerformanceNavigationTiming.responseStart - requestStart` |
| FCP (First Contentful Paint) | `PerformanceObserver` on `paint` |
| LCP (Largest Contentful Paint) | `PerformanceObserver` on `largest-contentful-paint` |
| Time-to-interactive proxy | time from `navigationStart` to the first frame where `[data-testid=viewer-view]` renders non-empty content |
| JS heap after landing | `performance.memory.usedJSHeapSize` sampled 3× after a 2 s settle |
| Resource count / bytes | `performance.getEntriesByType("resource")` total + sum of `transferSize` |
| Number of script chunks | `document.querySelectorAll("script[src]").length` |
| Main JS bundle size | `stats.json` from Vite build (`rollup-plugin-visualizer`) |

**Chrome-interaction metrics**, captured over scripted Playwright interactions (open geometries panel, drag sidebar handle, drag icon to bottom bar, open rooms panel, switch rooms):

| Metric | Source |
|---|---|
| Long tasks during interaction | `PerformanceObserver` on `longtask` |
| React commit count per interaction | `Profiler` `onRender` callback (temporary wrapper, removed before commit) |
| Frame time during drag | `requestAnimationFrame` timestamps sampled for the duration of a pointermove drag |
| JS heap delta per interaction | `performance.memory.usedJSHeapSize` before/after |

**R3F / WebGL viewer metrics** — profile the actual 3D path, because chrome fixes alone won't fix viewer-bound slowness:

| Metric | Source |
|---|---|
| Frame time at idle | average frame time with the viewer open, no interaction, measured via `requestAnimationFrame` loop over 5 s |
| Frame time during orbit | same, while dispatching synthetic mouse-drag events on the canvas |
| Frame time during frame scrub | same, while stepping `vis.step` via the store at 30 Hz |
| Frame time during splitter drag | same, while dragging a dockview splitter adjacent to the viewer |
| Draw calls per frame | `renderer.info.render.calls` from the R3F `<Canvas>` (access via a dev-only probe) |
| GPU memory / texture count | `renderer.info.memory.geometries` + `.textures` |
| Dispose-correctness | count geometries/textures before + after room switch; assert no growth on repeat switches |

The WebGL probe is a dev-only component that reads `useThree().gl.info` once per second and sends to the perf harness. Removed before any fix PR merges.

Capture runs: 5 per mode, median reported. Scripts land in `frontend/tests/perf/capture.ts`. Raw JSON saved to `/tmp/perf-audit/<timestamp>.json`.

### Step 2 — Bundle inspection

On the prod build, run `bun run build` with `rollup-plugin-visualizer` enabled and inspect the tree-map. The point is to find chunks that are larger than they need to be on the critical path, regardless of what `main` looks like.

Expected suspects (pre-diagnosis; audit will confirm or reject):

- `@mui/x-data-grid` in the `mui` chunk (used only by `GeometryPanel` and the `/rooms/` `DataGrid`, but loaded on every initial page).
- `dockview-react` + its CSS as eager import — needed eagerly; investigate whether the CSS can be trimmed.
- `plotly.js-dist-min` is already its own chunk (`plotly`) — verify it still lazy-loads and isn't pulled in by a non-plot code path.
- `ketcher` / SMILES stack — already lazy, verify no eager import slipped in.
- `three` chunk — inspect what's being imported from `three` and whether we're pulling in unused loaders (STLLoader, FBXLoader, etc.) when only GLTF/ASE-derived paths are used.
- `@react-three/drei` — known to have broad barrel exports; check that only used helpers land in the chunk.

### Step 3 — Runtime profile

For the Time-to-Interactive window and each scripted chrome interaction:
- Chrome DevTools Performance trace (automated via Playwright `tracing-start` / `tracing-stop`).
- React Profiler on the landing flow + each interaction (`Profiler` wrappers around `DockviewLayout`, `LandingPage`, and each `SidebarZone` for the duration of the audit; removed before commit).
- Count `commits` per component and cumulative `actualDuration` over the interaction.

For the R3F viewer:
- Enable `renderer.debug.checkShaderErrors` + capture one `renderer.info` snapshot per second.
- Three.js `stats.js` / `r3f-perf` drop-in overlay for the duration of a 30 s sampled session covering idle / orbit / scrub / splitter-drag.
- If a specific R3F component is hot in React Profiler, move its subscription pattern from hook-per-component to a store-driven `useFrame` side-effect to avoid re-render cascades.
- Check for classic mistakes: materials created inside `useFrame`, `new THREE.Vector3()` per frame in tight loops, lights re-registered on re-render, `<Canvas>` remounting on theme change.

### Step 4 — Rank the bottlenecks

For each identified contributor, record:

- `surface`: initial load / chrome interaction / R3F viewer.
- `metric`: which of the metrics in step 1 it dominates.
- `measurement`: the observed absolute number (e.g. "420 KB", "85 ms per commit", "18 ms / frame during scrub").
- `root cause`: the specific import / hook / listener / geometry allocation.
- `proposed fix`: the code change (not written yet).
- `expected improvement`: estimate (e.g. "shaves 40 KB from initial JS", "cuts React commits per panel-open from 12 to 3", "drops scrub frame time from 18 ms to 6 ms").
- `risk`: low / medium / high.

Output: `docs/superpowers/specs/2026-04-17-dockview-perf-audit-results.md` (results doc, written after the audit runs).

## Likely fixes (pre-audit, to verify)

The audit might confirm these; they're listed so implementation planning can start in parallel, but none land before the numbers justify them.

### Fix A — Lazy-load `PANELS[*].component`

Wrap every tool panel component in `React.lazy` at the registry level; mount sites (`SidebarZone`, `BottomZone`) wrap in `<Suspense fallback={<PanelSkeleton />}>`. Expected: moves `x-data-grid`, the plotly-touching code paths, the chat-message stack, and the filesystem panel into on-demand chunks. No functional change.

### Fix B — Split the `mui` manual chunk

Today's `vite.config.ts` bundles `@mui/material`, `@mui/icons-material`, `@mui/x-data-grid`, `@mui/x-tree-view` into one chunk. Split: `mui-core` (material + icons, loaded everywhere) vs `mui-grid` (`x-data-grid` + `x-tree-view`, lazy-loaded with the panels that need them).

### Fix C — Passive window listeners

`ActivityBar.tsx` `window.addEventListener("dragstart"/"dragend", …)` — pass `{ passive: true }`. Three instances × two listeners = six non-passive listeners today. Modest but free.

### Fix D — Defer `useColorScheme` re-renders

`DockviewLayout.tsx` reads `mode`, `systemMode` from `useColorScheme()` on every render. Re-renders propagate to the `DockviewReact` theme prop. Verify with the profiler whether these re-renders are cheap or cascading; if cascading, memo the `dockTheme` derivation on `resolvedMode`.

### Fix E — Dev-mode console errors

Per the screenshot session, two `401 /api/v1/users/me` and one `"constraints" key not found` error fire on every landing. These don't regress user-visible perf directly, but they block network pipelining and add flames to the React Profiler trace. The 401s are benign (guest flow); the `constraints` error is a real bug in `frontend/src/myapi/client.ts:40` and should be silenced or fixed.

### Fix F — R3F subscription patterns

Hypothesis: individual R3F components subscribe to zustand slices that fire on every frame scrub, causing a React commit cascade during playback. The fix is to read hot state (current frame, selection, step) inside `useFrame` from the store's `getState()` without subscribing, and let Three.js own the per-frame update directly. Confirm in the profiler before applying.

### Fix G — `<Canvas>` remount on theme flip

If step 3 finds `<Canvas>` remounts when MUI light/dark toggles, that's a WebGL-context-loss event dressed up as a re-render. Fix: stabilize the `Canvas` key + isolate theme-dependent DOM from the R3F tree.

## Minified / code-splitting proposal

Requested explicitly: "consider splitting the minified versions for faster initial loads".

The audit's bundle-diff step will surface the concrete splits. Pre-diagnosis, the expected shape:

- **Initial critical chunk:** React, React-DOM, React Router, Zustand, `dockview-react`, `@mui/material` core, `@mui/icons-material` (used site-wide).
- **Panel chunks (lazy):** one per panel — `panel-geometries` (pulls `x-data-grid`), `panel-plots-browser`, `panel-rooms`, `panel-filesystem`, `panel-chat`, `panel-selections`, `panel-modifiers`, `panel-analysis`.
- **Heavy feature chunks (lazy, already in place):** `three`, `plotly`, `ketcher`.

The `bun run build` target already emits chunks; the change is moving panels from eager imports in `registry.tsx` into `React.lazy` boundaries so Vite / Rollup actually split them.

## CI perf budget

Add a GitHub Actions (or equivalent) job that runs `frontend/tests/perf/capture.ts` against the current branch's prod build and fails if any measured value exceeds a **fixed budget** committed as `frontend/tests/perf/budget.json`:

- Initial JS bytes ≤ budget.
- JS heap after landing ≤ budget.
- LCP ≤ budget.
- Idle viewer frame time ≤ 16.7 ms (60 fps).
- Scrub-interaction frame time ≤ 33 ms (30 fps) as a floor.

Budgets are absolute targets, not deltas from a baseline — raising them requires an explicit commit with a justification in the PR body. A regressing PR cannot slip through by "just below the delta threshold."

On first land, the budget file is seeded from the audit results (the post-fix numbers, not the pre-fix ones).

## Migration & validation

- Audit script lands first as a standalone PR (no behavior change): `frontend/tests/perf/`, including the WebGL / R3F probe.
- Results doc lands second, after the audit runs.
- Fix PRs land third onwards, each with before/after measurements captured by the same harness.
- CI perf budget lands alongside the last fix PR (so the budget captures the target state, not an arbitrary snapshot).

## Open questions

- **Budget storage:** commit `budget.json` in the repo vs store in S3 / a GitHub artifact? Default: commit — small JSON, easy to diff in PRs.
- **Granularity of the CI budget:** per-metric thresholds vs one aggregate score? Default: per-metric (clearer on what regressed).
- **R3F probe packaging:** ship the `useThree`-based probe as a dev-only component (gated on `import.meta.env.DEV`) vs as a separate entry point? Default: dev-only, gated — zero prod cost.
