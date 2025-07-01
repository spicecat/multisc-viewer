import { defineConfig } from "vite";
import { VitePluginNode } from "vite-plugin-node";
import tsconfigPaths from "vite-tsconfig-paths";

export default defineConfig({
  server: {
    port: 5000,
  },
  build: {
    target: "ESNext",
  },
  plugins: [
    tsconfigPaths(),
    ...VitePluginNode({
      adapter: "nest",
      appPath: "src/main.ts",
      tsCompiler: "swc",
      swcOptions: {
        configFile: "config/server/.swcrc",
      },
    }),
  ],
});
