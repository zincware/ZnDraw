import { vi } from "vitest";

// Mock the panels registry to avoid loading all panel components (which pull
// in Three.js, plotly, etc.) via the activityBarSlice → registry chain.
vi.mock("./panels/registry", () => ({
	PANELS: {},
	TOOL_PANEL_IDS: [],
	SIDEBAR_DEFAULT_PX: 600,
	SIDEBAR_MIN_PX: 200,
	SIDEBAR_MAX_PX: 900,
	BOTTOM_DEFAULT_PX: 260,
	BOTTOM_MIN_PX: 120,
	BOTTOM_MAX_PX: 560,
	getDefaultsForBar: () => [],
}));

// Mock Three.js examples add-ons that reference canvas/window at module load time.
vi.mock("three/examples/jsm/Addons.js", () => ({}));
