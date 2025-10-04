import react from "@vitejs/plugin-react-swc";
import { defineConfig } from "vite";

export default defineConfig({
	plugins: [react()],
	base: "/",
	build: {
		outDir: '../../src/zndraw/static',
		emptyOutDir: true,
	},
	root: "./src",
	publicDir: "public", // Directory for static assets
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
