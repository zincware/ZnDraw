import react from "@vitejs/plugin-react-swc";
import { defineConfig } from "vite";

export default defineConfig({
  plugins: [react()],
  optimizeDeps: {
    include: [
      'ketcher-react',
      'ketcher-core',
      'ketcher-standalone'
    ],
    esbuildOptions: {
      define: {
        global: 'globalThis'
      },
    },
  },
  base: "/",
  build: {
    outDir: "../../src/zndraw/static",
    emptyOutDir: true,
    commonjsOptions: {
      transformMixedEsModules: true,
    },
  },
  root: "./src",
  publicDir: "public",
  resolve: {
    dedupe: ['react', 'react-dom', '@emotion/react', '@emotion/styled'],
  },
  define: {
    // Required for Ketcher - polyfills for browser environment
    'process.env.NODE_ENV': JSON.stringify(process.env.NODE_ENV || 'production'),
    global: 'globalThis',
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
