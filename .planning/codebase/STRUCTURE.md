# Codebase Structure

**Analysis Date:** 2026-03-05

## Directory Layout

```
zndraw-fastapi/
├── src/zndraw/                  # Python backend package
│   ├── __init__.py              # Public API: ZnDraw, Extension, Category
│   ├── _version.py              # Auto-generated version (hatch-vcs)
│   ├── app.py                   # FastAPI app assembly + Socket.IO mount
│   ├── config.py                # Pydantic Settings (ZNDRAW_ env prefix)
│   ├── database.py              # Lifespan: DB, Redis, storage, broker init
│   ├── dependencies.py          # FastAPI DI: RedisDep, StorageDep, WritableRoomDep, etc.
│   ├── models.py                # SQLModel tables: Room, Message, RoomGeometry, etc.
│   ├── schemas.py               # Pydantic request/response schemas
│   ├── socket_events.py         # Pydantic models for Socket.IO events
│   ├── socketio.py              # Socket.IO handlers (connect, room_join, etc.)
│   ├── redis.py                 # RedisKey namespace patterns
│   ├── exceptions.py            # RFC 9457 ProblemType definitions
│   ├── responses.py             # Custom response classes (MessagePackResponse)
│   ├── client.py                # Python ZnDraw client (MutableSequence[ase.Atoms])
│   ├── accessors.py             # Client accessor classes (Selections, Geometries, etc.)
│   ├── executor.py              # InternalExtensionExecutor for built-in extensions
│   ├── broker.py                # TaskIQ broker config
│   ├── result_backends.py       # CompositeResultBackend (Redis + Storage)
│   ├── cli.py                   # Typer CLI (zndraw, zndraw-db commands)
│   ├── server_manager.py        # Server process management helpers
│   ├── connectivity.py          # Atomic bond connectivity (vesin/networkx)
│   ├── enrichment.py            # Frame data enrichment (bond detection, etc.)
│   ├── materials.py             # MaterialProp Pydantic model
│   ├── transformations.py       # Transform/InArrayTransform for dynamic data refs
│   ├── tqdm.py                  # ZnDrawTqdm (server-side progress via Redis)
│   ├── routes/                  # FastAPI router modules
│   │   ├── __init__.py
│   │   ├── admin.py             # User management, server shutdown
│   │   ├── auth.py              # Guest/JWT login, registration
│   │   ├── bookmarks.py         # Frame bookmark CRUD
│   │   ├── chat.py              # Chat message CRUD with pagination
│   │   ├── edit_lock.py         # Distributed edit lock (Redis TTL)
│   │   ├── figures.py           # Plotly figure CRUD
│   │   ├── frames.py            # Frame read/write (msgpack binary)
│   │   ├── geometries.py        # Geometry CRUD + default camera
│   │   ├── presets.py           # Visual preset CRUD + apply
│   │   ├── problems.py          # RFC 9457 problem type docs
│   │   ├── progress.py          # tqdm-style progress trackers
│   │   ├── rooms.py             # Room CRUD, listing, duplication
│   │   ├── screenshots.py       # Screenshot capture/upload/download
│   │   ├── selection_groups.py  # Named selection group CRUD
│   │   ├── server_settings.py   # Default room, global config
│   │   ├── step.py              # Current frame step get/set
│   │   ├── tools.py             # RDKit molecule image generation
│   │   ├── trajectory.py        # File upload/download (ASE formats)
│   │   └── utility.py           # Version endpoint
│   ├── geometries/              # Pydantic geometry models
│   │   ├── __init__.py          # Registry dict + re-exports
│   │   ├── base.py              # BaseGeometry + property type aliases
│   │   ├── arrow.py             # Arrow geometry
│   │   ├── bonds.py             # Bond geometry
│   │   ├── box.py               # Box geometry
│   │   ├── camera.py            # Camera geometry (CameraType enum)
│   │   ├── cell.py              # Unit cell geometry
│   │   ├── circle_curve.py      # CircleCurve geometry
│   │   ├── curve.py             # Curve geometry (CurveMarker)
│   │   ├── floor.py             # Floor plane geometry
│   │   ├── fog.py               # Fog effect
│   │   ├── lights.py            # DirectionalLight, AmbientLight, HemisphereLight
│   │   ├── pathtracing.py       # PathTracing settings
│   │   ├── plane.py             # Plane geometry
│   │   ├── property_inspector.py # PropertyInspector settings
│   │   ├── shape.py             # Shape geometry
│   │   └── sphere.py            # Sphere (particle) geometry
│   ├── extensions/              # Compute extensions
│   │   ├── __init__.py          # Re-exports all extensions
│   │   ├── abc.py               # Extension base class + Category enum
│   │   ├── analysis.py          # Distance, DihedralAngle, Properties1D/2D
│   │   ├── modifiers.py         # Center, Delete, Duplicate, Replicate, etc.
│   │   ├── molecule_building.py # AddFromSMILES, PackBox
│   │   ├── selections.py        # All, Invert, Neighbour, Random, Range, etc.
│   │   └── filesystem.py        # Filesystem provider extension
│   ├── storage/                 # Frame storage backends
│   │   ├── __init__.py          # Re-exports AsebytesStorage, StorageRouter
│   │   ├── asebytes_backend.py  # AsyncBlobIO wrapper (memory/LMDB/MongoDB/H5MD)
│   │   └── router.py            # StorageRouter (provider mount support)
│   ├── providers/               # Data provider system
│   │   ├── __init__.py
│   │   ├── filesystem.py        # Filesystem data provider
│   │   └── frame_source.py      # FrameSource interface
│   ├── cli_agent/               # Agent CLI (zndraw-cli)
│   │   ├── __init__.py          # Typer app with sub-commands
│   │   ├── admin.py             # Admin commands
│   │   ├── auth.py              # Auth commands
│   │   ├── bookmarks.py         # Bookmark commands
│   │   ├── chat.py              # Chat commands
│   │   ├── connection.py        # Connection helpers
│   │   ├── extensions.py        # Extension commands
│   │   ├── figures.py           # Figure commands
│   │   ├── frames.py            # Frame commands
│   │   ├── geometries.py        # Geometry commands
│   │   ├── gif.py               # GIF export commands
│   │   ├── jobs.py              # Job commands
│   │   ├── mount.py             # Mount command
│   │   ├── output.py            # Output formatting helpers
│   │   ├── presets.py           # Preset commands
│   │   ├── rooms.py             # Room commands
│   │   ├── screenshots.py       # Screenshot commands
│   │   ├── selection.py         # Selection commands
│   │   ├── selection_groups.py  # Selection group commands
│   │   ├── sessions.py          # Session commands
│   │   └── step.py              # Step commands
│   ├── presets/                 # Built-in visual presets
│   │   └── __init__.py
│   ├── jobs/                    # Job-related utilities
│   └── static/                  # Built frontend assets (generated by vite build)
│       ├── index.html
│       └── assets/              # Hashed JS/CSS bundles
├── frontend/                    # React frontend
│   ├── package.json             # Dependencies (React 19, Three.js, MUI, Zustand, etc.)
│   ├── tsconfig.json            # TypeScript config
│   ├── vite.config.ts           # Vite build config + dev proxy
│   ├── playwright.config.ts     # Playwright E2E test config
│   ├── bun-env.d.ts             # Bun type declarations
│   ├── src/
│   │   ├── index.html           # HTML shell with #root div
│   │   ├── index.css            # Global CSS
│   │   ├── frontend.tsx         # React entry point (createRoot)
│   │   ├── App.tsx              # App shell: theme, query client, router
│   │   ├── socket.ts            # Socket.IO singleton + connectWithAuth
│   │   ├── store.tsx            # Main Zustand store (AppState = all slices)
│   │   ├── formStore.ts         # Extension panel UI state (immer)
│   │   ├── roomsStore.tsx       # Room list state for overview page
│   │   ├── pages/
│   │   │   ├── landingPage.tsx      # Main room view (3D scene + sidebar)
│   │   │   ├── roomList.tsx         # Room list/overview page
│   │   │   ├── templateSelection.tsx # Landing/template selection page
│   │   │   ├── filesystemBrowser.tsx # Filesystem browser page
│   │   │   └── cliLoginApprove.tsx  # CLI login approval page
│   │   ├── stores/
│   │   │   ├── slices/
│   │   │   │   ├── connectionSlice.ts  # roomId, user, session, connection state
│   │   │   │   ├── playbackSlice.ts    # frame, playing, bookmarks, fps
│   │   │   │   ├── sceneSlice.ts       # geometries, selections, mode, drawing/editing
│   │   │   │   ├── lockSlice.ts        # edit lock state + renewal timers
│   │   │   │   └── uiSlice.ts          # chat, snackbar, progress, screenshot
│   │   │   ├── geometryStore.ts    # Geometry panel mode/selection (immer)
│   │   │   └── windowManagerStore.ts # Floating window positions/sizes
│   │   ├── myapi/
│   │   │   └── client.ts            # REST API client (axios + all endpoint functions)
│   │   ├── hooks/
│   │   │   ├── useSocketManager.ts      # Socket.IO event wiring + room lifecycle
│   │   │   ├── useCameraControls.ts     # OrbitControls integration
│   │   │   ├── useDefaultCamera.ts      # Default camera query
│   │   │   ├── useDragAndDrop.tsx       # File drag-and-drop upload
│   │   │   ├── useFrameBatch.ts         # Frame data fetching with batching
│   │   │   ├── useFrameEditing.ts       # Frame editing mode helpers
│   │   │   ├── useFrameLoadTime.ts      # Frame load time tracking
│   │   │   ├── useFramePrefetcher.ts    # Prefetch adjacent frames
│   │   │   ├── useGeometries.ts         # Geometry data queries
│   │   │   ├── useGeometryCameraSync.ts # Sync camera geometry to Three.js
│   │   │   ├── useGeometryEditing.ts    # Geometry editing helpers
│   │   │   ├── useGeometryPersistence.ts # Auto-save geometry changes
│   │   │   ├── useJobs.ts              # Job listing + schema queries
│   │   │   ├── useKeyboardShortcuts.ts  # Global keyboard shortcuts
│   │   │   ├── useMoleculeImage.ts      # RDKit molecule image generation
│   │   │   ├── usePathtracingMesh.ts    # Path tracing mesh updates
│   │   │   ├── usePropertyInspector.ts  # Property inspector data
│   │   │   ├── usePropertyInspectorSettings.ts # PI settings management
│   │   │   ├── useRegisterFrameKeys.ts  # Register geometry frame data keys
│   │   │   ├── useSchemas.ts            # Extension schema queries
│   │   │   ├── useStepControl.ts        # Frame step synchronization
│   │   │   ├── useChat.ts              # Chat messages (infinite query)
│   │   │   └── useFigures.ts           # Figure data queries
│   │   ├── components/
│   │   │   ├── Canvas.tsx              # R3F Canvas with geometry dispatch
│   │   │   ├── CanvasErrorState.tsx     # Error fallback for canvas
│   │   │   ├── CanvasLoadingState.tsx   # Loading skeleton for canvas
│   │   │   ├── SideBar.tsx             # Left sidebar (nav + panel)
│   │   │   ├── PrimaryDrawer.tsx       # Icon navigation drawer
│   │   │   ├── SecondaryPanel.tsx      # Extension form panel (JSONForms)
│   │   │   ├── SelectionsPanel.tsx     # Selection management panel
│   │   │   ├── SelectionLayer.tsx      # Selection overlay
│   │   │   ├── SelectionGroupsPanel.tsx # Named selection groups
│   │   │   ├── SelectionTrackOverlay.tsx # Selection track visualization
│   │   │   ├── ProgressBar.tsx         # Bottom frame slider
│   │   │   ├── BookmarkLayer.tsx       # Bookmark indicators on slider
│   │   │   ├── FrameReference.tsx      # Frame reference display
│   │   │   ├── FrameSelectionInput.tsx # Frame selection input
│   │   │   ├── ChatWindow.tsx          # Chat window (lazy loaded)
│   │   │   ├── ConnectionDialog.tsx    # Python connection info dialog
│   │   │   ├── ConnectionStatus.tsx    # Connection status indicator
│   │   │   ├── LoginDialog.tsx         # Login form dialog
│   │   │   ├── RegisterDialog.tsx      # Registration form dialog
│   │   │   ├── AdminPanel.tsx          # Admin user management panel
│   │   │   ├── UserProfileDialog.tsx   # User profile/password change
│   │   │   ├── RoomManagementMenu.tsx  # Room actions dropdown
│   │   │   ├── DuplicateRoomDialog.tsx # Room duplication dialog
│   │   │   ├── FigureWindow.tsx        # Plotly figure display
│   │   │   ├── WindowManager.tsx       # Floating window container
│   │   │   ├── AddPlotButton.tsx       # Add plot button
│   │   │   ├── DropOverlay.tsx         # Drag-and-drop overlay
│   │   │   ├── ProgressNotifications.tsx # Progress toast notifications
│   │   │   ├── JobHistoryPanel.tsx     # Task history list
│   │   │   ├── JobListItem.tsx         # Single job list item
│   │   │   ├── JobStatusChips.tsx      # Job status indicator chips
│   │   │   ├── PathTracingRenderer.tsx # Path tracing toggle
│   │   │   ├── PathtracingUpdater.tsx  # Path tracing update trigger
│   │   │   ├── SiMGenButtons.tsx       # SiMGen integration buttons
│   │   │   ├── SiMGenTutorialDialog.tsx # SiMGen tutorial iframe
│   │   │   ├── CameraManager.js       # Camera control manager
│   │   │   ├── three/                  # Three.js geometry components
│   │   │   │   ├── Particles.tsx       # Sphere/atom instanced rendering
│   │   │   │   ├── Bonds.tsx           # Bond instanced rendering
│   │   │   │   ├── Arrow.tsx           # Arrow instanced rendering
│   │   │   │   ├── Box.tsx             # Box instanced rendering
│   │   │   │   ├── Plane.tsx           # Plane instanced rendering
│   │   │   │   ├── Shape.tsx           # Shape instanced rendering
│   │   │   │   ├── Camera.tsx          # Camera viewport rendering
│   │   │   │   ├── Cell.tsx            # Unit cell wireframe
│   │   │   │   ├── Curve.tsx           # Spline curve rendering
│   │   │   │   ├── CircleCurve.tsx     # Circle curve rendering
│   │   │   │   ├── Floor.tsx           # Floor plane
│   │   │   │   ├── Fog.tsx             # Fog effect
│   │   │   │   ├── Lights.tsx          # Light components
│   │   │   │   ├── DrawingIndicator.tsx # Drawing mode indicator
│   │   │   │   ├── EditingIndicator.tsx # Editing mode indicator
│   │   │   │   ├── MultiGeometryTransformControls.tsx # Transform gizmos
│   │   │   │   ├── HoverInfoBox.tsx    # Hover tooltip
│   │   │   │   ├── StaticInfoBox.tsx   # Static info display
│   │   │   │   ├── crosshair.tsx       # Crosshair overlay
│   │   │   │   ├── KeyboardShortcutsHandler.tsx # In-canvas keyboard
│   │   │   │   ├── ScreenshotProvider.tsx # Screenshot capture
│   │   │   │   ├── PathtracingScreenshotCapture.tsx # PT screenshot
│   │   │   │   ├── VirtualCanvas.tsx   # Offscreen canvas
│   │   │   │   ├── materials.tsx       # Material definitions
│   │   │   │   └── GeometryErrorBoundary.tsx # Error boundary
│   │   │   ├── geometry/               # Geometry panel components
│   │   │   │   ├── GeometryPanel.tsx   # Main geometry panel
│   │   │   │   ├── GeometryGrid.tsx    # Geometry list grid
│   │   │   │   ├── GeometryForm.tsx    # Geometry edit form
│   │   │   │   └── DeleteConfirmDialog.tsx # Delete confirmation
│   │   │   ├── filesystem/             # Filesystem browser components
│   │   │   │   ├── index.ts            # Barrel export
│   │   │   │   ├── FilesystemSelector.tsx # FS selector component
│   │   │   │   ├── FileList.tsx        # File listing
│   │   │   │   ├── FileBreadcrumbs.tsx # Path breadcrumbs
│   │   │   │   └── LoadFileDialog.tsx  # Load file dialog
│   │   │   ├── jsonforms-renderers/    # Custom JSONForms renderers
│   │   │   │   ├── ArrayEditorDialog.tsx # Array editing dialog
│   │   │   │   ├── ArrayFieldToolbar.tsx # Array field toolbar
│   │   │   │   ├── CustomColorPicker.tsx # Color picker renderer
│   │   │   │   ├── CustomDynamicEnumWithColorPicker.tsx # Dynamic enum + color
│   │   │   │   ├── CustomRangeSlider.tsx # Range slider renderer
│   │   │   │   ├── CustomSmilesEditor.tsx # SMILES string editor
│   │   │   │   ├── CustomSmilesPackEditor.tsx # SMILES pack editor
│   │   │   │   ├── DynamicEnumRenderer.tsx # Dynamic enum dropdown
│   │   │   │   ├── LightPositionRenderer.tsx # Light position editor
│   │   │   │   ├── MaterialEditor.tsx  # Three.js material selector
│   │   │   │   ├── OwnershipToggleRenderer.tsx # Owner toggle
│   │   │   │   ├── PositionAttachmentRenderer.tsx # Position attachment
│   │   │   │   ├── PropertyInspectorRenderer.tsx # Property inspector
│   │   │   │   ├── SmilesEditDialog.tsx # SMILES edit dialog (Ketcher)
│   │   │   │   ├── StaticValueDisplay.tsx # Static value display
│   │   │   │   ├── StringEnumControl.tsx # String enum control
│   │   │   │   ├── TransformEditor.tsx # Transform editor
│   │   │   │   ├── Vec3Renderer.tsx    # 3D vector input
│   │   │   │   ├── Vertices2DRenderer.tsx # 2D vertices editor
│   │   │   │   └── components/
│   │   │   │       └── PropertySelector.tsx # Property selector
│   │   │   └── shared/                 # Shared UI components
│   │   │       ├── FormLabelWithHelp.tsx # Label with help tooltip
│   │   │       ├── LoadingSkeletons.tsx # Loading skeleton components
│   │   │       └── MoleculePreview.tsx # Molecule image preview
│   │   ├── types/
│   │   │   ├── chat.ts                 # Chat message types
│   │   │   ├── jobs.ts                 # Job/task status types
│   │   │   └── property-inspector.ts   # Property inspector types
│   │   ├── utils/
│   │   │   ├── auth.ts                 # JWT auth (login, register, acquireToken)
│   │   │   ├── msgpack-numpy.ts        # Msgpack encoder/decoder (numpy compat)
│   │   │   ├── cameraUtils.ts          # Camera position resolution
│   │   │   ├── colorUtils.ts           # Color conversion utilities
│   │   │   ├── convertInstancedMesh.ts # Instanced mesh conversion
│   │   │   ├── geometryData.ts         # Geometry data helpers
│   │   │   ├── geometryDefaults.ts     # Geometry default values
│   │   │   ├── geometryEditing.ts      # Geometry editing utilities
│   │   │   ├── formStorage.ts          # Extension form localStorage
│   │   │   ├── jsonforms.ts            # JSONForms renderer config
│   │   │   ├── keyboard.ts             # Keyboard shortcut definitions
│   │   │   ├── jobUtils.ts             # Job utility functions
│   │   │   ├── propertyFormatting.ts   # Property display formatting
│   │   │   ├── roomTracking.ts         # Last visited room tracking
│   │   │   ├── roomValidation.ts       # Room ID validation
│   │   │   ├── screenshot.ts           # Screenshot download helper
│   │   │   ├── threeObjectPools.ts     # Three.js object pooling
│   │   │   ├── transformProcessor.ts   # Transform data processing
│   │   │   ├── versionCompatibility.ts # Client/server version check
│   │   │   ├── bookmarks.ts            # Bookmark utilities
│   │   │   ├── arrayEditor.ts          # Array editing utilities
│   │   │   └── remark-frame-link.js    # Markdown frame link plugin
│   │   ├── constants/
│   │   │   ├── fileTypes.ts            # Accepted file types
│   │   │   └── layout.ts              # Layout constants (heights, widths)
│   │   └── zndraw/
│   │       └── static/                # Static assets (favicon, etc.)
│   └── e2e/                           # Playwright E2E tests
│       ├── helpers.ts                 # Test helper utilities
│       ├── camera-session.spec.ts
│       ├── chat-features.spec.ts
│       ├── constraint-visualization.spec.ts
│       ├── docs-screenshots.spec.ts
│       ├── editing.spec.ts
│       ├── extensions-analysis.spec.ts
│       ├── frame-invalidation.spec.ts
│       ├── frames-navigation.spec.ts
│       ├── geometry-drawing.spec.ts
│       ├── registration.spec.ts
│       ├── socket-sync.spec.ts
│       └── ui-panels-chat.spec.ts
├── tests/                           # Python tests (pytest)
│   ├── conftest.py                  # Shared fixtures
│   ├── test_*.py                    # Test modules (~60 files)
│   ├── test_cli_agent/              # Agent CLI tests
│   └── worker/                      # Worker/task tests
├── docker/                          # Docker configurations
│   ├── production/                  # Multi-container production setup
│   ├── standalone/                  # Single-container setup
│   └── templates/                   # Docker Compose templates
├── docs/                            # Documentation (Sphinx)
│   └── source/                      # RST source files
├── data/                            # Sample data files
│   └── frames.lmdb/                 # LMDB frame storage
├── pyproject.toml                   # Python package config (hatchling)
├── main.py                          # Dev runner (uvicorn with reload)
├── Dockerfile                       # Container build
├── hatch_build.py                   # Custom build hook (frontend build)
├── CLAUDE.md -> AGENTS.md           # AI coding guidelines
├── AGENTS.md                        # AI coding guidelines
└── uv.lock                          # UV lockfile
```

## Directory Purposes

**`src/zndraw/`:**
- Purpose: Python backend package (FastAPI server + client library)
- Contains: All backend code, models, routes, storage, extensions
- Key files: `app.py` (entry), `database.py` (lifespan), `client.py` (Python client)

**`src/zndraw/routes/`:**
- Purpose: FastAPI API routers, one file per resource
- Contains: 18 router modules, each exports a `router` variable
- Pattern: All routes prefixed with `/v1/`

**`src/zndraw/geometries/`:**
- Purpose: Pydantic models for 3D scene objects
- Contains: One file per geometry type + base class + registry
- Key file: `__init__.py` with `geometries` dict mapping type name -> class

**`src/zndraw/extensions/`:**
- Purpose: Compute operations that run on the server
- Contains: Extension base class + category-organized implementations
- Key file: `abc.py` with `Extension` base class and `Category` enum

**`src/zndraw/storage/`:**
- Purpose: Pluggable frame storage backends
- Contains: `AsebytesStorage` (actual storage) + `StorageRouter` (mount/provider layer)

**`src/zndraw/cli_agent/`:**
- Purpose: Structured JSON CLI for LLM agents
- Contains: One Typer sub-app per resource, mirrors REST API

**`src/zndraw/providers/`:**
- Purpose: Data provider interfaces for external frame sources
- Contains: Filesystem provider, frame source base class

**`frontend/src/`:**
- Purpose: React + Three.js frontend application
- Contains: Pages, components, stores, hooks, utilities, types

**`frontend/src/pages/`:**
- Purpose: Top-level page components (one per route)
- Contains: 5 page files matching router paths

**`frontend/src/stores/`:**
- Purpose: Zustand store slices and secondary stores
- Contains: 5 state slices (connection, playback, scene, lock, UI) + geometry store + window manager store

**`frontend/src/hooks/`:**
- Purpose: Custom React hooks for data fetching and state management
- Contains: ~25 hooks for frames, geometries, camera, sockets, keyboard, etc.

**`frontend/src/components/`:**
- Purpose: React UI components
- Contains: Top-level components + subdirectories for three.js, geometry panel, filesystem, jsonforms renderers, shared

**`frontend/src/components/three/`:**
- Purpose: Three.js geometry rendering components (inside R3F Canvas)
- Contains: One component per geometry type + camera, lights, indicators, overlays

**`frontend/src/components/jsonforms-renderers/`:**
- Purpose: Custom JSONForms renderers for specialized form inputs
- Contains: Color pickers, SMILES editors, vec3 inputs, material selectors, etc.

**`frontend/src/myapi/`:**
- Purpose: REST API client wrapper
- Contains: Single `client.ts` file with axios instance + all API functions

**`frontend/src/utils/`:**
- Purpose: Utility functions (auth, encoding, geometry, formatting)
- Contains: ~20 utility modules

**`frontend/src/types/`:**
- Purpose: TypeScript type definitions
- Contains: Chat, jobs, property inspector types

**`frontend/e2e/`:**
- Purpose: Playwright end-to-end tests
- Contains: ~13 spec files testing UI workflows

**`tests/`:**
- Purpose: Python backend tests (pytest)
- Contains: ~60 test files + conftest fixtures
- Subdirectories: `test_cli_agent/`, `worker/`

**`docker/`:**
- Purpose: Docker deployment configurations
- Contains: Production (multi-container), standalone, and template configs

## Key File Locations

**Entry Points:**
- `src/zndraw/app.py`: FastAPI + Socket.IO ASGI app (`socket_app`)
- `main.py`: Dev runner (uvicorn with reload)
- `frontend/src/frontend.tsx`: React entry point (createRoot)
- `frontend/src/App.tsx`: App shell (providers, router)

**Configuration:**
- `src/zndraw/config.py`: Backend settings (env vars)
- `frontend/vite.config.ts`: Vite build + dev proxy config
- `frontend/tsconfig.json`: TypeScript compiler options
- `pyproject.toml`: Python package, dependencies, entry points

**Core Backend Logic:**
- `src/zndraw/database.py`: Application lifespan (resource init/teardown)
- `src/zndraw/dependencies.py`: FastAPI dependency injection
- `src/zndraw/socketio.py`: Socket.IO event handlers
- `src/zndraw/models.py`: SQL data models
- `src/zndraw/schemas.py`: API request/response schemas

**Core Frontend Logic:**
- `frontend/src/store.tsx`: Main Zustand store
- `frontend/src/hooks/useSocketManager.ts`: Socket.IO real-time sync
- `frontend/src/myapi/client.ts`: REST API client
- `frontend/src/socket.ts`: Socket.IO singleton
- `frontend/src/utils/auth.ts`: Authentication logic

**3D Rendering:**
- `frontend/src/components/Canvas.tsx`: R3F Canvas + geometry dispatch
- `frontend/src/components/three/Particles.tsx`: Atom rendering
- `frontend/src/components/three/Bonds.tsx`: Bond rendering

**Testing:**
- `tests/conftest.py`: Shared pytest fixtures (server, client, Redis)
- `frontend/e2e/helpers.ts`: Playwright test helpers
- `frontend/playwright.config.ts`: Playwright configuration

## Naming Conventions

**Files (Backend):**
- `snake_case.py`: All Python modules
- Route files: resource name singular (`frames.py`, `chat.py`, `step.py`)
- Geometry files: geometry name lowercase (`sphere.py`, `camera.py`, `bonds.py`)

**Files (Frontend):**
- `camelCase.tsx` / `camelCase.ts`: All TypeScript/React files
- Hooks: `use{Name}.ts` (e.g., `useSocketManager.ts`, `useFrameBatch.ts`)
- Pages: `camelCase.tsx` in `pages/` (e.g., `landingPage.tsx`, `roomList.tsx`)
- Components: `PascalCase.tsx` (e.g., `Canvas.tsx`, `SideBar.tsx`)
- Store slices: `{name}Slice.ts` (e.g., `connectionSlice.ts`)
- Three.js components: `PascalCase.tsx` matching geometry type (e.g., `Particles.tsx`, `Arrow.tsx`)
- Utils: `camelCase.ts` (e.g., `auth.ts`, `colorUtils.ts`)

**Directories:**
- Backend: `snake_case` (e.g., `cli_agent/`, `selection_groups.py`)
- Frontend: `camelCase` or `kebab-case` (e.g., `jsonforms-renderers/`, `three/`)

**Python Naming:**
- Classes: `PascalCase` (e.g., `RoomGeometry`, `BaseGeometry`, `WritableRoomDep`)
- Functions: `snake_case` (e.g., `get_redis`, `verify_room`)
- FastAPI deps: `PascalCase` aliases with `Dep` suffix (e.g., `RedisDep`, `StorageDep`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `PROBLEM_TYPES`, `WORKER_EMAIL`)

**TypeScript Naming:**
- React components: `PascalCase` (e.g., `MainPage`, `GeometryPanel`)
- Hooks: `usePascalCase` (e.g., `useSocketManager`, `useFrameBatch`)
- Interfaces/Types: `PascalCase` (e.g., `ConnectionSlice`, `AppState`)
- API functions: `camelCase` (e.g., `listGeometries`, `createRoom`)
- Store actions: `camelCase` (e.g., `setRoomId`, `updateGeometry`)
- Constants: `UPPER_SNAKE_CASE` (e.g., `LAYOUT_CONSTANTS`, `TRAJECTORY_ACCEPT`)

## Where to Add New Code

**New REST Endpoint:**
1. Create route module: `src/zndraw/routes/{resource}.py`
2. Define request/response schemas in `src/zndraw/schemas.py`
3. Add ProblemType if needed in `src/zndraw/exceptions.py`
4. Register router in `src/zndraw/app.py`
5. Add tests in `tests/test_{resource}.py`

**New Geometry Type:**
1. Create Pydantic model: `src/zndraw/geometries/{name}.py` (extend `BaseGeometry`)
2. Register in `src/zndraw/geometries/__init__.py` (`geometries` dict + imports)
3. Create Three.js component: `frontend/src/components/three/{Name}.tsx`
4. Register in `frontend/src/components/Canvas.tsx` (component dispatch map)
5. Add custom JSONForms renderers if needed in `frontend/src/components/jsonforms-renderers/`

**New Extension:**
1. Create extension class: add to `src/zndraw/extensions/{category}.py` (extend `Extension`, set `category`)
2. Register in category dict (e.g., `modifiers`, `selections`, `analysis`)
3. Re-export in `src/zndraw/extensions/__init__.py`
4. Frontend auto-discovers via job schema API (no frontend code needed)
5. Add tests in `tests/`

**New Frontend Page:**
1. Create page component: `frontend/src/pages/{pageName}.tsx`
2. Add route in `frontend/src/App.tsx` (router config)
3. Add hooks as needed in `frontend/src/hooks/`

**New Frontend Component:**
1. Create component file: `frontend/src/components/{ComponentName}.tsx`
2. For Three.js components: place in `frontend/src/components/three/`
3. For shared UI: place in `frontend/src/components/shared/`

**New Zustand State:**
1. Add to existing slice in `frontend/src/stores/slices/{slice}Slice.ts`
2. Or create new slice file, compose into `store.tsx`

**New Custom Hook:**
1. Create: `frontend/src/hooks/use{Name}.ts`
2. Pattern: use React Query for data fetching, Zustand for state

**New API Function:**
1. Add to `frontend/src/myapi/client.ts` (follow existing pattern)
2. Define TypeScript interfaces for request/response

**New SQL Model:**
1. Add to `src/zndraw/models.py` (extend `SQLModel, table=True`)
2. Import at module level to register with `SQLModel.metadata`
3. Tables auto-created on startup (`init_db_on_startup=True`)

**New Socket.IO Event:**
1. Define Pydantic model in `src/zndraw/socket_events.py`
2. Add handler in `src/zndraw/socketio.py` or route that broadcasts
3. Add frontend listener in `frontend/src/hooks/useSocketManager.ts`

## Special Directories

**`src/zndraw/static/`:**
- Purpose: Built frontend assets served by FastAPI in production
- Generated: Yes (by `vite build` output to this directory)
- Committed: No (`.gitignore`d, built during package build via `hatch_build.py`)

**`data/`:**
- Purpose: Sample data files for development
- Generated: No
- Committed: Yes (contains `frames.lmdb/`)

**`docker/`:**
- Purpose: Docker deployment configurations
- Generated: No
- Committed: Yes

**`.planning/`:**
- Purpose: GSD planning documents
- Generated: Yes (by GSD commands)
- Committed: Yes

**`zndraw-media/`:**
- Purpose: Screenshot/media storage (per-room directories)
- Generated: Yes (by screenshot capture)
- Committed: No (`.gitignore`d)

**`frontend/node_modules/`:**
- Purpose: NPM/Bun dependencies
- Generated: Yes
- Committed: No

**`tmp/`:**
- Purpose: Temporary files during development
- Generated: Yes
- Committed: No

---

*Structure analysis: 2026-03-05*
