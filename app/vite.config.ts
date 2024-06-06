import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  build: {
    outDir: '../zndraw/templates',  // Output directory for templates
    emptyOutDir: true,  // Clear the output directory before building
  },
  publicDir: '../zndraw/static',  // Optional: directory for static assets
});
