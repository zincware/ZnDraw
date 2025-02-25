import react from "@vitejs/plugin-react-swc";
import { defineConfig } from "vite";

// https://vitejs.dev/config/
export default defineConfig({
	plugins: [react()],
	base: "/",
	build: {
		outDir: "../zndraw/templates", // Output directory for templates
		emptyOutDir: true, // Clear the output directory before building
	},
	publicDir: "public", // Directory for static assets
	server: {
		proxy: {
			"/reset": "http://localhost:1234",
			"/token": "http://localhost:1234",
			"/upload": "http://localhost:1234",
			"/download": "http://localhost:1234",
			"/socket.io": {
				target: "ws://localhost:1234",
				ws: true,
			},
		},
	},
});
