import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  build: {
    outDir: "../zndraw/templates", // Output directory for templates
    emptyOutDir: true, // Clear the output directory before building
  },
  publicDir: "../zndraw/static", // Optional: directory for static assets
  server: {
    proxy: {
      // string shorthand: http://localhost:5173/foo -> http://localhost:4567/foo
      "/reset": "http://localhost:12134",
      // // Proxying websockets or socket.io: ws://localhost:5173/socket.io -> ws://localhost:5174/socket.io
      "/socket.io": {
        target: "ws://localhost:1234",
        ws: true,
      },
    },
  },
});
