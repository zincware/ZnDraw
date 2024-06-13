import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  build: {
    outDir: "../zndraw/templates", // Output directory for templates
    emptyOutDir: true, // Clear the output directory before building
  },
  publicDir: "public", // Directory for static assets
  server: {
    proxy: {
      "/reset": "http://localhost:1234",
      "/token": "http://localhost:1234",
      "/socket.io": {
        target: "ws://localhost:1234",
        ws: true,
      },
    },
  },
});
