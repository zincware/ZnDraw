# Constraint Visualization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reintroduce constraint visualization via a default `constraints-fixed-atoms` geometry and enhance the TransformEditor with a constraint selector.

**Architecture:** Add InArrayTransform-based default geometry to room creation. Enhance TransformEditor to fetch constraint data from the current frame and present a dropdown (no free-solo) when `source === "constraints"`.

**Tech Stack:** Python (Pydantic, FastAPI), TypeScript (React, MUI, react-query, Playwright)

---

### Task 1: Add `constraints-fixed-atoms` Default Geometry

**Files:**
- Modify: `src/zndraw/routes/rooms.py:66-109`

**Step 1: Add imports**

At the top of `rooms.py`, add:

```python
from zndraw.transformations import InArrayTransform
from zndraw.materials import MeshBasicMaterial
```

**Step 2: Add default geometry entry**

In `_initialize_default_geometries()`, add to the `defaults` dict after `"floor"`:

```python
"constraints-fixed-atoms": (
    "Sphere",
    {
        "active": True,
        "position": InArrayTransform(
            source="constraints",
            path="0.kwargs.indices",
            filter="arrays.positions",
        ),
        "radius": InArrayTransform(
            source="constraints",
            path="0.kwargs.indices",
            filter="arrays.radii",
        ),
        "color": ["#FF0000"],
        "material": MeshBasicMaterial(wireframe=True),
        "scale": [(0.71, 0.71, 0.71)],
        "selecting": {"enabled": False},
        "hovering": {"enabled": False},
    },
),
```

Note: The `Sphere` model accepts `InArrayTransform` and `MeshBasicMaterial` objects directly — Pydantic handles serialization via `model_dump_json()`.

**Step 3: Run existing tests to verify no regressions**

Run: `uv run pytest tests/test_constraints.py tests/test_rooms.py -v`
Expected: All tests pass. No behavioral change to existing rooms.

**Step 4: Commit**

```bash
git add src/zndraw/routes/rooms.py
git commit -m "feat: add constraints-fixed-atoms default geometry"
```

---

### Task 2: Write Backend Test for Default Constraint Geometry

**Files:**
- Modify: `tests/test_constraints.py` (add at end of file)

**Step 1: Write tests**

Add these tests to `tests/test_constraints.py`:

```python
# =============================================================================
# Section 4: Default geometry tests
# =============================================================================


def test_default_constraint_geometry_created_on_room_creation(server: str):
    """New rooms include a constraints-fixed-atoms geometry."""
    import uuid
    from zndraw import ZnDraw

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    geo = client.geometries["constraints-fixed-atoms"]
    assert geo is not None

    # Verify it's a Sphere with InArrayTransform position
    from zndraw.geometries import Sphere
    from zndraw.transformations import InArrayTransform

    assert isinstance(geo, Sphere)
    assert isinstance(geo.position, InArrayTransform)
    assert geo.position.source == "constraints"
    assert geo.position.path == "0.kwargs.indices"
    assert geo.position.filter == "arrays.positions"

    assert isinstance(geo.radius, InArrayTransform)
    assert geo.radius.source == "constraints"
    assert geo.radius.filter == "arrays.radii"

    assert geo.color == ["#FF0000"]
    assert geo.selecting.enabled is False
    assert geo.hovering.enabled is False

    client.disconnect()


def test_constraint_geometry_renders_fixed_atoms(server: str):
    """Constraint geometry filters positions to only fixed atoms."""
    import uuid
    from zndraw import ZnDraw

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    atoms = ase.Atoms("H5", positions=[[i, 0, 0] for i in range(5)])
    atoms.set_constraint(FixAtoms(indices=[1, 3]))
    client.append(atoms)

    # The constraint geometry exists and has correct transform config
    geo = client.geometries["constraints-fixed-atoms"]
    from zndraw.transformations import InArrayTransform

    assert isinstance(geo.position, InArrayTransform)
    # Verify constraint data roundtrip
    retrieved = client[0]
    assert len(retrieved.constraints) == 1
    assert isinstance(retrieved.constraints[0], FixAtoms)
    np.testing.assert_array_equal(retrieved.constraints[0].index, [1, 3])

    client.disconnect()
```

**Step 2: Run tests**

Run: `uv run pytest tests/test_constraints.py -v -k "default_constraint_geometry or constraint_geometry_renders"`
Expected: Both tests pass.

**Step 3: Commit**

```bash
git add tests/test_constraints.py
git commit -m "test: add backend tests for default constraint geometry"
```

---

### Task 3: Enhance TransformEditor with Constraint Selector

**Files:**
- Modify: `frontend/src/components/jsonforms-renderers/TransformEditor.tsx`

**Step 1: Add imports and constraint fetching**

Replace the entire `TransformEditor.tsx` with an enhanced version that:
- Takes `roomId` and `currentFrame` as additional props (or reads from store)
- When `source === "constraints"`, fetches constraint data via `getFrames`
- Shows an MUI `Select` dropdown with constraint entries
- On selection, sets `path` to `{index}.kwargs.indices`

```tsx
import CloseIcon from "@mui/icons-material/Close";
import FilterAltIcon from "@mui/icons-material/FilterAlt";
import {
	Box,
	CircularProgress,
	FormControl,
	IconButton,
	InputLabel,
	MenuItem,
	Paper,
	Select,
	TextField,
	Tooltip,
	Typography,
} from "@mui/material";
import { useQuery } from "@tanstack/react-query";
import { getFrames } from "../../myapi/client";
import { useAppStore } from "../../store";
import type { Transform } from "../../utils/transformProcessor";

interface TransformEditorProps {
	value: Transform;
	label: string;
	required?: boolean;
	onChange: (newValue: Transform) => void;
	onClear: () => void;
}

interface ConstraintEntry {
	index: number;
	name: string;
	kwargs: Record<string, any>;
}

/**
 * Parse constraint data from frame into selectable entries.
 */
function parseConstraints(data: any): ConstraintEntry[] {
	if (!Array.isArray(data)) return [];
	return data
		.map((item: any, index: number) => {
			if (
				typeof item === "object" &&
				item !== null &&
				typeof item.name === "string" &&
				typeof item.kwargs === "object"
			) {
				return { index, name: item.name, kwargs: item.kwargs };
			}
			return null;
		})
		.filter((entry: ConstraintEntry | null): entry is ConstraintEntry =>
			entry !== null
		);
}

/**
 * Format a constraint entry for display in the selector.
 */
function formatConstraintLabel(entry: ConstraintEntry): string {
	const indices = entry.kwargs.indices;
	if (Array.isArray(indices)) {
		const preview =
			indices.length > 5
				? `[${indices.slice(0, 5).join(", ")}, ...]`
				: `[${indices.join(", ")}]`;
		return `#${entry.index}: ${entry.name} — indices: ${preview}`;
	}
	return `#${entry.index}: ${entry.name}`;
}

/**
 * TransformEditor component - inline editor for transform objects.
 *
 * When source is "constraints", provides a dropdown to select from
 * actual constraint entries in the current frame. Otherwise shows
 * text fields for source, path, and filter.
 */
export default function TransformEditor({
	value,
	label,
	required,
	onChange,
	onClear,
}: TransformEditorProps) {
	const transform = value || {
		type: "in_array" as const,
		source: "",
		path: "",
		filter: "",
	};

	const roomId = useAppStore((state) => state.roomId);
	const currentFrame = useAppStore((state) => state.currentFrame);

	// Fetch constraint data when source is "constraints"
	const isConstraintSource = transform.source === "constraints";
	const { data: constraintData, isLoading: isLoadingConstraints } = useQuery({
		queryKey: ["frame", roomId, currentFrame, "constraints"],
		queryFn: ({ signal }) =>
			getFrames(roomId!, currentFrame, ["constraints"], signal),
		enabled: isConstraintSource && !!roomId,
	});

	const constraints = isConstraintSource
		? parseConstraints(constraintData?.constraints)
		: [];

	// Extract current selection index from path (e.g., "0.kwargs.indices" → 0)
	const selectedIndex = (() => {
		const match = transform.path.match(/^(\d+)\.kwargs\.indices$/);
		return match ? Number.parseInt(match[1], 10) : -1;
	})();

	const updateField = (field: string, newValue: string) => {
		onChange({
			...transform,
			[field]: newValue,
		} as Transform);
	};

	const handleConstraintSelect = (index: number) => {
		onChange({
			...transform,
			path: `${index}.kwargs.indices`,
		} as Transform);
	};

	return (
		<Box sx={{ marginBottom: 2 }}>
			<Box
				sx={{
					display: "flex",
					alignItems: "center",
					justifyContent: "space-between",
					marginBottom: 1,
				}}
			>
				<Box sx={{ display: "flex", alignItems: "center", gap: 1 }}>
					<FilterAltIcon fontSize="small" color="primary" />
					<Typography variant="subtitle2" color="primary">
						{label}
						{required && " *"} (Transform)
					</Typography>
				</Box>
				<Tooltip title="Remove transform and use simple key">
					<IconButton size="small" onClick={onClear}>
						<CloseIcon fontSize="small" />
					</IconButton>
				</Tooltip>
			</Box>

			<Paper elevation={1} sx={{ padding: 2, backgroundColor: "#f9f9f9" }}>
				<Box sx={{ display: "flex", flexDirection: "column", gap: 1.5 }}>
					<TextField
						fullWidth
						size="small"
						label="Source"
						value={transform.source || ""}
						onChange={(e) => updateField("source", e.target.value)}
						placeholder="e.g., constraints"
						helperText="Frame data key containing indices"
					/>

					{isConstraintSource ? (
						<FormControl fullWidth size="small">
							<InputLabel>Constraint</InputLabel>
							<Select
								value={selectedIndex >= 0 ? selectedIndex : ""}
								label="Constraint"
								onChange={(e) =>
									handleConstraintSelect(e.target.value as number)
								}
								disabled={isLoadingConstraints}
								endAdornment={
									isLoadingConstraints ? (
										<CircularProgress size={20} sx={{ mr: 2 }} />
									) : undefined
								}
							>
								{constraints.length === 0 && !isLoadingConstraints ? (
									<MenuItem disabled>
										No constraints in current frame
									</MenuItem>
								) : (
									constraints.map((entry) => (
										<MenuItem key={entry.index} value={entry.index}>
											{formatConstraintLabel(entry)}
										</MenuItem>
									))
								)}
							</Select>
						</FormControl>
					) : (
						<TextField
							fullWidth
							size="small"
							label="Path"
							value={transform.path || ""}
							onChange={(e) => updateField("path", e.target.value)}
							placeholder="e.g., 0.kwargs.indices"
							helperText="Dot-separated path to extract indices"
						/>
					)}

					<TextField
						fullWidth
						size="small"
						label="Filter"
						value={transform.filter || ""}
						onChange={(e) => updateField("filter", e.target.value)}
						placeholder="e.g., arrays.positions"
						helperText="Frame data key to filter"
					/>
				</Box>
			</Paper>
		</Box>
	);
}
```

Key changes from original:
- Added `useAppStore` to read `roomId` and `currentFrame`
- Added `useQuery` to fetch constraints when `source === "constraints"`
- Added `parseConstraints()` and `formatConstraintLabel()` helpers
- Replaced `path` TextField with MUI `Select` when in constraint mode
- No free-solo: only selectable options from actual constraint data
- Shows "No constraints in current frame" when empty
- Shows loading spinner while fetching

**Step 2: Verify frontend builds**

Run: `cd frontend && bun run build`
Expected: Build succeeds with no type errors.

**Step 3: Commit**

```bash
git add frontend/src/components/jsonforms-renderers/TransformEditor.tsx
git commit -m "feat: add constraint selector to TransformEditor"
```

---

### Task 4: Write Playwright E2E Tests for Constraint Visualization

**Files:**
- Create: `frontend/e2e/constraint-visualization.spec.ts`

**Step 1: Write the E2E test**

```typescript
import { test, expect } from "@playwright/test";
import { BASE_URL, CLI, PY, waitForScene } from "./helpers";

const ROOM = "test-constraints";

/**
 * Open the geometry panel by clicking the sidebar "Manage geometries" button.
 */
async function openGeometryPanel(page: import("@playwright/test").Page) {
	const manageBtn = page.getByRole("button", { name: "Manage geometries" });
	await manageBtn.click();
	await page.waitForSelector('[role="grid"]', {
		state: "visible",
		timeout: 10000,
	});
}

function setupConstraintRoom() {
	CLI(`rooms create --room-id ${ROOM}`);
	PY(`
from zndraw import ZnDraw
import ase
from ase.constraints import FixAtoms, FixedLine

vis = ZnDraw(url='${BASE_URL}', room='${ROOM}')
del vis[:]

# 5 atoms: fix atoms 0 and 2, constrain atom 1 to a line
atoms = ase.Atoms('H5', positions=[(0,0,0),(2,0,0),(4,0,0),(6,0,0),(8,0,0)])
atoms.set_constraint([
    FixAtoms(indices=[0, 2]),
    FixedLine(1, direction=[1, 0, 0]),
])
vis.append(atoms)
`);
}

test.describe("Constraint Visualization", () => {
	test.describe.configure({ mode: "serial" });

	test.beforeAll(() => {
		setupConstraintRoom();
	});

	test("constraints-fixed-atoms geometry appears in panel", async ({
		page,
	}) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		await openGeometryPanel(page);
		const grid = page.locator('[role="grid"]');
		await expect(
			grid
				.getByText("constraints-fixed-atoms", { exact: true })
				.first(),
		).toBeVisible({ timeout: 5000 });

		await page.screenshot({
			path: "e2e/screenshots/constraint-geometry-panel.png",
		});
	});

	test("constraint geometry renders in scene", async ({ page }) => {
		await page.goto(`${BASE_URL}/rooms/${ROOM}`);
		await waitForScene(page);

		// Wait for the scene to fully render with constraints
		await page.waitForTimeout(3000);

		await page.screenshot({
			path: "e2e/screenshots/constraint-visualization.png",
		});
	});
});
```

**Step 2: Run the E2E test**

Run: `cd frontend && npx playwright test e2e/constraint-visualization.spec.ts`
Expected: Tests pass — constraint geometry visible in panel and scene.

**Step 3: Commit**

```bash
git add frontend/e2e/constraint-visualization.spec.ts
git commit -m "test: add Playwright E2E tests for constraint visualization"
```

---

### Task 5: Add Documentation with Molecular Example

**Files:**
- Modify: `docs/source/python-api.rst` (add constraint section)

**Step 1: Add constraint visualization section**

Add a section to `python-api.rst` (or the most appropriate existing doc page) showing:

```rst
Constraint Visualization
------------------------

ZnDraw automatically visualizes atomic constraints. When you upload atoms with
ASE constraints, a ``constraints-fixed-atoms`` geometry overlays red wireframe
spheres on the constrained atoms.

.. code-block:: python

   import ase
   from ase.constraints import FixAtoms
   from zndraw import ZnDraw

   # Create acetic acid and constrain the methyl group
   atoms = ase.Atoms(
       "C2H4O2",
       positions=[
           [0.0, 0.0, 0.0],   # C (methyl) — constrained
           [1.5, 0.0, 0.0],   # C (carboxyl)
           [-0.5, 1.0, 0.0],  # H — constrained
           [-0.5, -1.0, 0.0], # H — constrained
           [2.5, 1.0, 0.0],   # O (C=O)
           [2.5, -1.0, 0.0],  # O (OH)
       ],
   )
   atoms.set_constraint(FixAtoms(indices=[0, 2, 3]))

   vis = ZnDraw(url="http://localhost:8000")
   vis.append(atoms)

The constrained atoms (methyl C and two H atoms) will appear with a red
wireframe sphere overlay, while the COOH group remains undecorated.

**Customization:** Open the geometry panel and click ``constraints-fixed-atoms``
to change the color, scale, or target a different constraint. The Transform
Editor shows a dropdown of all constraints in the current frame — select one
to switch which atoms are highlighted.
```

**Step 2: Build docs**

Run: `cd docs && make html`
Expected: Docs build without warnings.

**Step 3: Commit**

```bash
git add docs/source/python-api.rst
git commit -m "docs: add constraint visualization guide with molecular example"
```

---

### Task 6: Final Verification

**Step 1: Run full backend test suite**

Run: `uv run pytest tests/ -v --timeout=900`
Expected: All tests pass (including new constraint geometry tests).

**Step 2: Run frontend build**

Run: `cd frontend && bun run build`
Expected: Clean build, no errors.

**Step 3: Run type checking**

Run: `uv run pyright .`
Expected: No new errors.

**Step 4: Run formatting**

Run: `uv run ruff format . && uv run ruff check --select I --fix .`
Expected: Clean.

**Step 5: Final commit (if any formatting changes)**

```bash
git add -A
git commit -m "style: apply formatting"
```
