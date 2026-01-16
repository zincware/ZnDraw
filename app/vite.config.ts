import react from "@vitejs/plugin-react-swc";
import { defineConfig } from "vite";
import { execSync } from "node:child_process";

/**
 * Gets the Python package version from the installed zndraw package.
 * Falls back to "0.0.0-dev" if unable to retrieve.
 */
function getPythonPackageVersion(): string {
	try {
		// Try to get version from installed package
		const version = execSync(
			'uv run python -c "from zndraw import __version__; print(__version__)"',
			{ encoding: "utf8", stdio: ["pipe", "pipe", "ignore"] },
		).trim();
		return version;
	} catch {
		// Fallback if unable to get version
		return "0.0.0-dev";
	}
}

export default defineConfig({
	plugins: [react()],
	optimizeDeps: {
		// Note: ketcher is NOT included here - it's lazy loaded and should not be pre-bundled
		esbuildOptions: {
			define: {
				global: "globalThis",
			},
		},
	},
	base: "/",
	build: {
		outDir: "../../src/zndraw/static",
		emptyOutDir: true,
		// Disable modulepreload for heavy chunks - they'll load on-demand
		modulePreload: {
			resolveDependencies: (filename, deps) => {
				// Filter out heavy chunks from preload - they should only load when needed
				return deps.filter(
					(dep) =>
						!dep.includes("ketcher") &&
						!dep.includes("plotly") &&
						!dep.includes("ChatWindow"),
				);
			},
		},
		commonjsOptions: {
			transformMixedEsModules: true,
		},
		rollupOptions: {
			output: {
				manualChunks: {
					// Heavy 3D rendering - loaded when entering a room
					three: ["three", "@react-three/fiber", "@react-three/drei"],
					// Note: ketcher is NOT listed here - it's lazy loaded via React.lazy()
					// and will naturally code-split with SmilesEditDialog
					// Charting library - loaded when viewing plots
					plotly: ["plotly.js", "plotly.js-dist-min"],
					// UI framework - used everywhere, but cacheable
					mui: [
						"@mui/material",
						"@mui/icons-material",
						"@mui/x-data-grid",
						"@mui/x-tree-view",
					],
					// React ecosystem - stable, highly cacheable
					vendor: [
						"react",
						"react-dom",
						"react-router-dom",
						"@tanstack/react-query",
						"zustand",
					],
				},
			},
		},
	},
	root: "./src",
	publicDir: "public",
	resolve: {
		dedupe: ["react", "react-dom", "@emotion/react", "@emotion/styled"],
	},
	define: {
		// Required for Ketcher - polyfills for browser environment
		"process.env.NODE_ENV": JSON.stringify(
			process.env.NODE_ENV || "production",
		),
		global: "globalThis",
		// Inject version in dev mode if not already set
		"import.meta.env.VITE_APP_VERSION": JSON.stringify(
			process.env.VITE_APP_VERSION || getPythonPackageVersion(),
		),
	},
	server: {
		proxy: {
			"/api": "http://localhost:5000",
			"/socket.io": {
				target: "ws://localhost:5000",
				ws: true,
			},
		},
	},
});
