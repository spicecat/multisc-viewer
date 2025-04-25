// config/client/vite.config.mts
import { existsSync, readFileSync } from "fs";
import { basename, dirname, resolve } from "path";
import preprocessor, { sveltePreprocess } from "file:///home/ahteh/Developer/plot-viewer/node_modules/svelte-preprocess/dist/index.js";
import { compile, compileModule, preprocess } from "file:///home/ahteh/Developer/plot-viewer/node_modules/svelte/src/compiler/index.js";
import { defineConfig } from "file:///home/ahteh/Developer/plot-viewer/node_modules/vite/dist/node/index.js";
var vite_config_default = defineConfig({
  resolve: {
    conditions: ["svelte"]
  },
  plugins: [
    /* @__PURE__ */ (() => {
      let isSSR = false;
      return {
        name: "svelte",
        enforce: "pre",
        async resolveId(source, importer, options) {
          let matches;
          if (source === "__route.svelte") {
            if (!importer) throw new Error("__route.svelte used as entry");
            if ((matches = /__client\/routes\/(.+)$/.exec(importer)) === null)
              throw new Error("Imported '__route.svelte' not from route");
            return `src/client/routes/${matches[1]}.svelte`;
          }
          if (source === "__app.svelte") {
            if (!importer) throw new Error("__app.svelte used as entry");
            return "src/client/templates/app.svelte";
          }
          if (source === "$meta") {
            if (!importer) throw new Error("$meta used as entry");
            return "src/client/templates/$meta.svelte.js";
          }
          if (options.isEntry) {
            return source;
          } else if (!basename(source).includes(".")) {
            const relative = source.startsWith("./") || source.startsWith("../");
            if (relative) {
              const path = resolve(dirname(importer), source) + ".ts";
              if (existsSync(path)) {
                return path;
              } else {
                return null;
              }
            } else if ((matches = /\$lib\/(.+)/.exec(source)) !== null) {
              const path = `src/client/lib/${matches[1]}.ts`;
              if (existsSync(path)) {
                return path;
              } else {
                return null;
              }
            } else {
              const path = `${source}.ts`;
              if (existsSync(path)) {
                return path;
              } else {
                return null;
              }
            }
          } else if (source.endsWith(".svelte") || source.endsWith(".ts")) {
            const relative = source.startsWith("./") || source.startsWith("../");
            if (relative) {
              return resolve(dirname(importer), source);
            } else if ((matches = /\$lib\/(.+)/.exec(source)) !== null) {
              return `src/client/lib/${matches[1]}`;
            } else {
              return source;
            }
          } else {
            return null;
          }
        },
        async load(id, options) {
          if (/__client\/routes\/(.+)$/.test(id)) {
            if (options?.ssr) {
              return readFileSync("src/client/templates/server.js").toString();
            } else {
              return readFileSync("src/client/templates/client.js").toString();
            }
          } else {
            return null;
          }
        },
        async transform(code, id, options) {
          if (!id.includes(".svelte")) return null;
          if (options?.ssr) isSSR = true;
          if (id.endsWith(".js")) {
            let matches = /.*\/([^/]+\.svelte)/.exec(id);
            if (!matches) throw new Error("Svelte matching error");
            const filename = matches[1];
            const preprocessed = await preprocess(
              code,
              preprocessor({
                typescript: {
                  compilerOptions: {
                    module: "es2020",
                    target: "es2020",
                    verbatimModuleSyntax: true
                  }
                }
              }),
              { filename }
            );
            matches = /src\/client\/(.+)\.svelte\.js/.exec(id);
            const result = compileModule(preprocessed.code, {
              generate: options?.ssr ? "server" : "client",
              filename: matches?.[1].split("/").at(-1) || "unknown"
            });
            return { ...result.js };
          } else {
            let matches = /.*\/([^/]+\.svelte)/.exec(id);
            if (!matches) throw new Error("Svelte matching error");
            const filename = matches[1];
            const preprocessed = await preprocess(
              code,
              sveltePreprocess({
                typescript: {
                  compilerOptions: {
                    module: "es2020",
                    target: "es2020",
                    verbatimModuleSyntax: true
                  }
                }
              }),
              { filename }
            );
            matches = /src\/client\/lib\/components\/(.+)\.svelte/.exec(id);
            const result = compile(preprocessed.code, {
              generate: options?.ssr ? "server" : "client",
              name: matches?.[1].split("/").at(-1) || "App",
              runes: true,
              css: "injected"
            });
            return { ...result.js };
          }
        },
        generateBundle(options, bundle) {
          if (isSSR) {
            Object.entries(bundle).forEach(([, chunk]) => {
              if (chunk.type === "chunk" && chunk.facadeModuleId) {
                const matches = /src\/__client\/routes\/(.+)/.exec(
                  chunk.facadeModuleId
                );
                if (matches) {
                  const route = matches[1];
                  chunk.imports.forEach((id) => {
                    if (id.startsWith("assets/")) {
                      const path = id.split("/").slice(1);
                      const routeNesting = route.split("/").length;
                      const correction = new Array(routeNesting).fill("..").join("/") + "/";
                      const file = path.at(-1);
                      if (!file)
                        this.error(
                          `Failed to correct import ${id} in route ${route}`
                        );
                      const pattern = new RegExp(
                        `import\\s*(?:\\{.*\\}\\s*from\\s*)?("./assets/${file}"|'./assets/${file}')`
                      );
                      chunk.code = chunk.code.replace(
                        pattern,
                        (match, subId) => {
                          return match.replace(
                            subId,
                            subId.replace(`./assets/${file}`, correction + id)
                          );
                        }
                      );
                    }
                  });
                  chunk.fileName = `routes/${route}.svelte.js`;
                }
              }
            });
          } else {
            Object.entries(bundle).forEach(([, chunk]) => {
              if (chunk.type === "chunk" && chunk.facadeModuleId) {
                const matches = /src\/__client\/routes\/(.+)/.exec(
                  chunk.facadeModuleId
                );
                if (matches) {
                  const route = matches[1];
                  chunk.imports.forEach((id) => {
                    if (id.startsWith("assets/")) {
                      const path = id.split("/").slice(1);
                      const importedNesting = path.length;
                      const routeNesting = route.split("/").length;
                      const file = path.at(-1);
                      const correction = routeNesting > importedNesting ? new Array(routeNesting - importedNesting).fill("..").join("/") + "/" : "./";
                      if (!file)
                        this.error(
                          `Failed to correct import ${id} in route ${route}`
                        );
                      const pattern = new RegExp(
                        `import\\s*(?:\\{.*\\}\\s*from\\s*)?("./${file}"|'./${file}')`
                      );
                      chunk.code = chunk.code.replace(
                        pattern,
                        (match, id2) => {
                          return match.replace(
                            id2,
                            id2.replace("./", correction)
                          );
                        }
                      );
                    }
                  });
                  chunk.fileName = `assets/${route}.svelte.js`;
                }
              }
            });
          }
        }
      };
    })()
  ]
});
export {
  vite_config_default as default
};
//# sourceMappingURL=data:application/json;base64,ewogICJ2ZXJzaW9uIjogMywKICAic291cmNlcyI6IFsiY29uZmlnL2NsaWVudC92aXRlLmNvbmZpZy5tdHMiXSwKICAic291cmNlc0NvbnRlbnQiOiBbImNvbnN0IF9fdml0ZV9pbmplY3RlZF9vcmlnaW5hbF9kaXJuYW1lID0gXCIvaG9tZS9haHRlaC9EZXZlbG9wZXIvcGxvdC12aWV3ZXIvY29uZmlnL2NsaWVudFwiO2NvbnN0IF9fdml0ZV9pbmplY3RlZF9vcmlnaW5hbF9maWxlbmFtZSA9IFwiL2hvbWUvYWh0ZWgvRGV2ZWxvcGVyL3Bsb3Qtdmlld2VyL2NvbmZpZy9jbGllbnQvdml0ZS5jb25maWcubXRzXCI7Y29uc3QgX192aXRlX2luamVjdGVkX29yaWdpbmFsX2ltcG9ydF9tZXRhX3VybCA9IFwiZmlsZTovLy9ob21lL2FodGVoL0RldmVsb3Blci9wbG90LXZpZXdlci9jb25maWcvY2xpZW50L3ZpdGUuY29uZmlnLm10c1wiO2ltcG9ydCB7IGV4aXN0c1N5bmMsIHJlYWRGaWxlU3luYyB9IGZyb20gXCJmc1wiO1xuaW1wb3J0IHsgYmFzZW5hbWUsIGRpcm5hbWUsIHJlc29sdmUgfSBmcm9tIFwicGF0aFwiO1xuaW1wb3J0IHByZXByb2Nlc3NvciwgeyBzdmVsdGVQcmVwcm9jZXNzIH0gZnJvbSBcInN2ZWx0ZS1wcmVwcm9jZXNzXCI7XG5pbXBvcnQgeyBjb21waWxlLCBjb21waWxlTW9kdWxlLCBwcmVwcm9jZXNzIH0gZnJvbSBcInN2ZWx0ZS9jb21waWxlclwiO1xuaW1wb3J0IHsgZGVmaW5lQ29uZmlnIH0gZnJvbSBcInZpdGVcIjtcblxuZXhwb3J0IGRlZmF1bHQgZGVmaW5lQ29uZmlnKHtcbiAgcmVzb2x2ZToge1xuICAgIGNvbmRpdGlvbnM6IFtcInN2ZWx0ZVwiXSxcbiAgfSxcbiAgcGx1Z2luczogW1xuICAgICgoKSA9PiB7XG4gICAgICBsZXQgaXNTU1IgPSBmYWxzZTtcblxuICAgICAgcmV0dXJuIHtcbiAgICAgICAgbmFtZTogXCJzdmVsdGVcIixcbiAgICAgICAgZW5mb3JjZTogXCJwcmVcIixcbiAgICAgICAgYXN5bmMgcmVzb2x2ZUlkKHNvdXJjZSwgaW1wb3J0ZXIsIG9wdGlvbnMpIHtcbiAgICAgICAgICBsZXQgbWF0Y2hlczogUmVnRXhwRXhlY0FycmF5IHwgbnVsbDtcbiAgICAgICAgICBpZiAoc291cmNlID09PSBcIl9fcm91dGUuc3ZlbHRlXCIpIHtcbiAgICAgICAgICAgIGlmICghaW1wb3J0ZXIpIHRocm93IG5ldyBFcnJvcihcIl9fcm91dGUuc3ZlbHRlIHVzZWQgYXMgZW50cnlcIik7XG5cbiAgICAgICAgICAgIGlmICgobWF0Y2hlcyA9IC9fX2NsaWVudFxcL3JvdXRlc1xcLyguKykkLy5leGVjKGltcG9ydGVyKSkgPT09IG51bGwpXG4gICAgICAgICAgICAgIHRocm93IG5ldyBFcnJvcihcIkltcG9ydGVkICdfX3JvdXRlLnN2ZWx0ZScgbm90IGZyb20gcm91dGVcIik7XG5cbiAgICAgICAgICAgIHJldHVybiBgc3JjL2NsaWVudC9yb3V0ZXMvJHttYXRjaGVzWzFdfS5zdmVsdGVgO1xuICAgICAgICAgIH1cblxuICAgICAgICAgIGlmIChzb3VyY2UgPT09IFwiX19hcHAuc3ZlbHRlXCIpIHtcbiAgICAgICAgICAgIGlmICghaW1wb3J0ZXIpIHRocm93IG5ldyBFcnJvcihcIl9fYXBwLnN2ZWx0ZSB1c2VkIGFzIGVudHJ5XCIpO1xuXG4gICAgICAgICAgICByZXR1cm4gXCJzcmMvY2xpZW50L3RlbXBsYXRlcy9hcHAuc3ZlbHRlXCI7XG4gICAgICAgICAgfVxuXG4gICAgICAgICAgaWYgKHNvdXJjZSA9PT0gXCIkbWV0YVwiKSB7XG4gICAgICAgICAgICBpZiAoIWltcG9ydGVyKSB0aHJvdyBuZXcgRXJyb3IoXCIkbWV0YSB1c2VkIGFzIGVudHJ5XCIpO1xuXG4gICAgICAgICAgICByZXR1cm4gXCJzcmMvY2xpZW50L3RlbXBsYXRlcy8kbWV0YS5zdmVsdGUuanNcIjtcbiAgICAgICAgICB9XG5cbiAgICAgICAgICBpZiAob3B0aW9ucy5pc0VudHJ5KSB7XG4gICAgICAgICAgICByZXR1cm4gc291cmNlO1xuICAgICAgICAgIH0gZWxzZSBpZiAoIWJhc2VuYW1lKHNvdXJjZSkuaW5jbHVkZXMoXCIuXCIpKSB7XG4gICAgICAgICAgICBjb25zdCByZWxhdGl2ZSA9XG4gICAgICAgICAgICAgIHNvdXJjZS5zdGFydHNXaXRoKFwiLi9cIikgfHwgc291cmNlLnN0YXJ0c1dpdGgoXCIuLi9cIik7XG5cbiAgICAgICAgICAgIGlmIChyZWxhdGl2ZSkge1xuICAgICAgICAgICAgICBjb25zdCBwYXRoID0gcmVzb2x2ZShkaXJuYW1lKGltcG9ydGVyISksIHNvdXJjZSkgKyBcIi50c1wiO1xuXG4gICAgICAgICAgICAgIGlmIChleGlzdHNTeW5jKHBhdGgpKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHBhdGg7XG4gICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgICAgICAgICAgLy8gdGhpcy5lcnJvcihgVW5hYmxlIHRvIHJlc29sdmUgJHtzb3VyY2V9LCBpbXBvcnRlZCBmcm9tICR7aW1wb3J0ZXJ9YCk7XG4gICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0gZWxzZSBpZiAoKG1hdGNoZXMgPSAvXFwkbGliXFwvKC4rKS8uZXhlYyhzb3VyY2UpKSAhPT0gbnVsbCkge1xuICAgICAgICAgICAgICBjb25zdCBwYXRoID0gYHNyYy9jbGllbnQvbGliLyR7bWF0Y2hlc1sxXX0udHNgO1xuXG4gICAgICAgICAgICAgIGlmIChleGlzdHNTeW5jKHBhdGgpKSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIHBhdGg7XG4gICAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgICAgcmV0dXJuIG51bGw7XG4gICAgICAgICAgICAgICAgLy8gdGhpcy5lcnJvcihgVW5hYmxlIHRvIHJlc29sdmUgJHtzb3VyY2V9LCBpbXBvcnRlZCBmcm9tICR7aW1wb3J0ZXJ9YCk7XG4gICAgICAgICAgICAgIH1cbiAgICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICAgIGNvbnN0IHBhdGggPSBgJHtzb3VyY2V9LnRzYDtcblxuICAgICAgICAgICAgICBpZiAoZXhpc3RzU3luYyhwYXRoKSkge1xuICAgICAgICAgICAgICAgIHJldHVybiBwYXRoO1xuICAgICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICAgIHJldHVybiBudWxsO1xuICAgICAgICAgICAgICAgIC8vIHRoaXMuZXJyb3IoYFVuYWJsZSB0byByZXNvbHZlICR7c291cmNlfSwgaW1wb3J0ZWQgZnJvbSAke2ltcG9ydGVyfWApO1xuICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9XG4gICAgICAgICAgfSBlbHNlIGlmIChzb3VyY2UuZW5kc1dpdGgoXCIuc3ZlbHRlXCIpIHx8IHNvdXJjZS5lbmRzV2l0aChcIi50c1wiKSkge1xuICAgICAgICAgICAgY29uc3QgcmVsYXRpdmUgPVxuICAgICAgICAgICAgICBzb3VyY2Uuc3RhcnRzV2l0aChcIi4vXCIpIHx8IHNvdXJjZS5zdGFydHNXaXRoKFwiLi4vXCIpO1xuXG4gICAgICAgICAgICBpZiAocmVsYXRpdmUpIHtcbiAgICAgICAgICAgICAgcmV0dXJuIHJlc29sdmUoZGlybmFtZShpbXBvcnRlciEpLCBzb3VyY2UpO1xuICAgICAgICAgICAgfSBlbHNlIGlmICgobWF0Y2hlcyA9IC9cXCRsaWJcXC8oLispLy5leGVjKHNvdXJjZSkpICE9PSBudWxsKSB7XG4gICAgICAgICAgICAgIHJldHVybiBgc3JjL2NsaWVudC9saWIvJHttYXRjaGVzWzFdfWA7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICByZXR1cm4gc291cmNlO1xuICAgICAgICAgICAgfVxuICAgICAgICAgIH0gZWxzZSB7XG4gICAgICAgICAgICByZXR1cm4gbnVsbDtcbiAgICAgICAgICB9XG4gICAgICAgIH0sXG4gICAgICAgIGFzeW5jIGxvYWQoaWQsIG9wdGlvbnMpIHtcbiAgICAgICAgICBpZiAoL19fY2xpZW50XFwvcm91dGVzXFwvKC4rKSQvLnRlc3QoaWQpKSB7XG4gICAgICAgICAgICBpZiAob3B0aW9ucz8uc3NyKSB7XG4gICAgICAgICAgICAgIC8vIHRoaXMuZXJyb3IoYFNTUiByZW5kZXJpbmcgcm91dGUgJHtpZH0gbm90IGNhcHR1cmVkIGJ5IHJlc29sdmVJZGApO1xuICAgICAgICAgICAgICByZXR1cm4gcmVhZEZpbGVTeW5jKFwic3JjL2NsaWVudC90ZW1wbGF0ZXMvc2VydmVyLmpzXCIpLnRvU3RyaW5nKCk7XG4gICAgICAgICAgICB9IGVsc2Uge1xuICAgICAgICAgICAgICByZXR1cm4gcmVhZEZpbGVTeW5jKFwic3JjL2NsaWVudC90ZW1wbGF0ZXMvY2xpZW50LmpzXCIpLnRvU3RyaW5nKCk7XG4gICAgICAgICAgICB9XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIHJldHVybiBudWxsO1xuICAgICAgICAgIH1cbiAgICAgICAgfSxcbiAgICAgICAgYXN5bmMgdHJhbnNmb3JtKGNvZGUsIGlkLCBvcHRpb25zKSB7XG4gICAgICAgICAgLy8gdGhpcy5pbmZvKGB0cmFuc2Zvcm1pbmcgJHtpZH1gKTtcbiAgICAgICAgICBpZiAoIWlkLmluY2x1ZGVzKFwiLnN2ZWx0ZVwiKSkgcmV0dXJuIG51bGw7XG4gICAgICAgICAgaWYgKG9wdGlvbnM/LnNzcikgaXNTU1IgPSB0cnVlO1xuXG4gICAgICAgICAgaWYgKGlkLmVuZHNXaXRoKFwiLmpzXCIpKSB7XG4gICAgICAgICAgICBsZXQgbWF0Y2hlcyA9IC8uKlxcLyhbXi9dK1xcLnN2ZWx0ZSkvLmV4ZWMoaWQpO1xuICAgICAgICAgICAgaWYgKCFtYXRjaGVzKSB0aHJvdyBuZXcgRXJyb3IoXCJTdmVsdGUgbWF0Y2hpbmcgZXJyb3JcIik7XG4gICAgICAgICAgICBjb25zdCBmaWxlbmFtZSA9IG1hdGNoZXNbMV07XG5cbiAgICAgICAgICAgIGNvbnN0IHByZXByb2Nlc3NlZCA9IGF3YWl0IHByZXByb2Nlc3MoXG4gICAgICAgICAgICAgIGNvZGUsXG4gICAgICAgICAgICAgIHByZXByb2Nlc3Nvcih7XG4gICAgICAgICAgICAgICAgdHlwZXNjcmlwdDoge1xuICAgICAgICAgICAgICAgICAgY29tcGlsZXJPcHRpb25zOiB7XG4gICAgICAgICAgICAgICAgICAgIG1vZHVsZTogXCJlczIwMjBcIixcbiAgICAgICAgICAgICAgICAgICAgdGFyZ2V0OiBcImVzMjAyMFwiLFxuICAgICAgICAgICAgICAgICAgICB2ZXJiYXRpbU1vZHVsZVN5bnRheDogdHJ1ZSxcbiAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgfSksXG4gICAgICAgICAgICAgIHsgZmlsZW5hbWUgfSxcbiAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgIG1hdGNoZXMgPSAvc3JjXFwvY2xpZW50XFwvKC4rKVxcLnN2ZWx0ZVxcLmpzLy5leGVjKGlkKTtcbiAgICAgICAgICAgIGNvbnN0IHJlc3VsdCA9IGNvbXBpbGVNb2R1bGUocHJlcHJvY2Vzc2VkLmNvZGUsIHtcbiAgICAgICAgICAgICAgZ2VuZXJhdGU6IG9wdGlvbnM/LnNzciA/IFwic2VydmVyXCIgOiBcImNsaWVudFwiLFxuICAgICAgICAgICAgICBmaWxlbmFtZTogbWF0Y2hlcz8uWzFdLnNwbGl0KFwiL1wiKS5hdCgtMSkgfHwgXCJ1bmtub3duXCIsXG4gICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgcmV0dXJuIHsgLi4ucmVzdWx0LmpzIH07XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIGxldCBtYXRjaGVzID0gLy4qXFwvKFteL10rXFwuc3ZlbHRlKS8uZXhlYyhpZCk7XG4gICAgICAgICAgICBpZiAoIW1hdGNoZXMpIHRocm93IG5ldyBFcnJvcihcIlN2ZWx0ZSBtYXRjaGluZyBlcnJvclwiKTtcbiAgICAgICAgICAgIGNvbnN0IGZpbGVuYW1lID0gbWF0Y2hlc1sxXTtcblxuICAgICAgICAgICAgY29uc3QgcHJlcHJvY2Vzc2VkID0gYXdhaXQgcHJlcHJvY2VzcyhcbiAgICAgICAgICAgICAgY29kZSxcbiAgICAgICAgICAgICAgc3ZlbHRlUHJlcHJvY2Vzcyh7XG4gICAgICAgICAgICAgICAgdHlwZXNjcmlwdDoge1xuICAgICAgICAgICAgICAgICAgY29tcGlsZXJPcHRpb25zOiB7XG4gICAgICAgICAgICAgICAgICAgIG1vZHVsZTogXCJlczIwMjBcIixcbiAgICAgICAgICAgICAgICAgICAgdGFyZ2V0OiBcImVzMjAyMFwiLFxuICAgICAgICAgICAgICAgICAgICB2ZXJiYXRpbU1vZHVsZVN5bnRheDogdHJ1ZSxcbiAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgfSksXG4gICAgICAgICAgICAgIHsgZmlsZW5hbWUgfSxcbiAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgIG1hdGNoZXMgPSAvc3JjXFwvY2xpZW50XFwvbGliXFwvY29tcG9uZW50c1xcLyguKylcXC5zdmVsdGUvLmV4ZWMoaWQpO1xuICAgICAgICAgICAgY29uc3QgcmVzdWx0ID0gY29tcGlsZShwcmVwcm9jZXNzZWQuY29kZSwge1xuICAgICAgICAgICAgICBnZW5lcmF0ZTogb3B0aW9ucz8uc3NyID8gXCJzZXJ2ZXJcIiA6IFwiY2xpZW50XCIsXG4gICAgICAgICAgICAgIG5hbWU6IG1hdGNoZXM/LlsxXS5zcGxpdChcIi9cIikuYXQoLTEpIHx8IFwiQXBwXCIsXG4gICAgICAgICAgICAgIHJ1bmVzOiB0cnVlLFxuICAgICAgICAgICAgICBjc3M6IFwiaW5qZWN0ZWRcIixcbiAgICAgICAgICAgIH0pO1xuXG4gICAgICAgICAgICByZXR1cm4geyAuLi5yZXN1bHQuanMgfTtcbiAgICAgICAgICB9XG4gICAgICAgIH0sXG4gICAgICAgIGdlbmVyYXRlQnVuZGxlKG9wdGlvbnMsIGJ1bmRsZSkge1xuICAgICAgICAgIGlmIChpc1NTUikge1xuICAgICAgICAgICAgT2JqZWN0LmVudHJpZXMoYnVuZGxlKS5mb3JFYWNoKChbLCBjaHVua10pID0+IHtcbiAgICAgICAgICAgICAgaWYgKGNodW5rLnR5cGUgPT09IFwiY2h1bmtcIiAmJiBjaHVuay5mYWNhZGVNb2R1bGVJZCkge1xuICAgICAgICAgICAgICAgIGNvbnN0IG1hdGNoZXMgPSAvc3JjXFwvX19jbGllbnRcXC9yb3V0ZXNcXC8oLispLy5leGVjKFxuICAgICAgICAgICAgICAgICAgY2h1bmsuZmFjYWRlTW9kdWxlSWQsXG4gICAgICAgICAgICAgICAgKTtcblxuICAgICAgICAgICAgICAgIGlmIChtYXRjaGVzKSB7XG4gICAgICAgICAgICAgICAgICBjb25zdCByb3V0ZSA9IG1hdGNoZXNbMV07XG4gICAgICAgICAgICAgICAgICAvLyB0aGlzLmluZm8oYHJvdXRlICR7cm91dGV9IGltcG9ydHM6ICR7SlNPTi5zdHJpbmdpZnkoY2h1bmsuaW1wb3J0cyl9YCk7XG5cbiAgICAgICAgICAgICAgICAgIGNodW5rLmltcG9ydHMuZm9yRWFjaCgoaWQpID0+IHtcbiAgICAgICAgICAgICAgICAgICAgLy8gdGhpcy5pbmZvKGByb3V0ZSAke3JvdXRlfSBpbXBvcnRpbmcgJHtpZH1gKTtcbiAgICAgICAgICAgICAgICAgICAgaWYgKGlkLnN0YXJ0c1dpdGgoXCJhc3NldHMvXCIpKSB7XG4gICAgICAgICAgICAgICAgICAgICAgY29uc3QgcGF0aCA9IGlkLnNwbGl0KFwiL1wiKS5zbGljZSgxKTtcbiAgICAgICAgICAgICAgICAgICAgICBjb25zdCByb3V0ZU5lc3RpbmcgPSByb3V0ZS5zcGxpdChcIi9cIikubGVuZ3RoO1xuXG4gICAgICAgICAgICAgICAgICAgICAgY29uc3QgY29ycmVjdGlvbiA9XG4gICAgICAgICAgICAgICAgICAgICAgICBuZXcgQXJyYXkocm91dGVOZXN0aW5nKS5maWxsKFwiLi5cIikuam9pbihcIi9cIikgKyBcIi9cIjtcbiAgICAgICAgICAgICAgICAgICAgICBjb25zdCBmaWxlID0gcGF0aC5hdCgtMSk7XG5cbiAgICAgICAgICAgICAgICAgICAgICBpZiAoIWZpbGUpXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVycm9yKFxuICAgICAgICAgICAgICAgICAgICAgICAgICBgRmFpbGVkIHRvIGNvcnJlY3QgaW1wb3J0ICR7aWR9IGluIHJvdXRlICR7cm91dGV9YCxcbiAgICAgICAgICAgICAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgICAgICAgICAgICBjb25zdCBwYXR0ZXJuID0gbmV3IFJlZ0V4cChcbiAgICAgICAgICAgICAgICAgICAgICAgIGBpbXBvcnRcXFxccyooPzpcXFxcey4qXFxcXH1cXFxccypmcm9tXFxcXHMqKT8oXCIuL2Fzc2V0cy8ke2ZpbGV9XCJ8Jy4vYXNzZXRzLyR7ZmlsZX0nKWAsXG4gICAgICAgICAgICAgICAgICAgICAgKTtcbiAgICAgICAgICAgICAgICAgICAgICBjaHVuay5jb2RlID0gY2h1bmsuY29kZS5yZXBsYWNlKFxuICAgICAgICAgICAgICAgICAgICAgICAgcGF0dGVybixcbiAgICAgICAgICAgICAgICAgICAgICAgIChtYXRjaCwgc3ViSWQ6IHN0cmluZykgPT4ge1xuICAgICAgICAgICAgICAgICAgICAgICAgICByZXR1cm4gbWF0Y2gucmVwbGFjZShcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdWJJZCxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgICBzdWJJZC5yZXBsYWNlKGAuL2Fzc2V0cy8ke2ZpbGV9YCwgY29ycmVjdGlvbiArIGlkKSxcbiAgICAgICAgICAgICAgICAgICAgICAgICAgKTtcbiAgICAgICAgICAgICAgICAgICAgICAgIH0sXG4gICAgICAgICAgICAgICAgICAgICAgKTtcbiAgICAgICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICAgICAgfSk7XG5cbiAgICAgICAgICAgICAgICAgIGNodW5rLmZpbGVOYW1lID0gYHJvdXRlcy8ke3JvdXRlfS5zdmVsdGUuanNgO1xuICAgICAgICAgICAgICAgIH1cbiAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgfSk7XG4gICAgICAgICAgfSBlbHNlIHtcbiAgICAgICAgICAgIE9iamVjdC5lbnRyaWVzKGJ1bmRsZSkuZm9yRWFjaCgoWywgY2h1bmtdKSA9PiB7XG4gICAgICAgICAgICAgIGlmIChjaHVuay50eXBlID09PSBcImNodW5rXCIgJiYgY2h1bmsuZmFjYWRlTW9kdWxlSWQpIHtcbiAgICAgICAgICAgICAgICBjb25zdCBtYXRjaGVzID0gL3NyY1xcL19fY2xpZW50XFwvcm91dGVzXFwvKC4rKS8uZXhlYyhcbiAgICAgICAgICAgICAgICAgIGNodW5rLmZhY2FkZU1vZHVsZUlkLFxuICAgICAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgICAgICBpZiAobWF0Y2hlcykge1xuICAgICAgICAgICAgICAgICAgY29uc3Qgcm91dGUgPSBtYXRjaGVzWzFdO1xuICAgICAgICAgICAgICAgICAgLy8gdGhpcy5pbmZvKGByb3V0ZSAke3JvdXRlfSBpbXBvcnRzOiAke0pTT04uc3RyaW5naWZ5KGNodW5rLmltcG9ydHMpfWApO1xuXG4gICAgICAgICAgICAgICAgICBjaHVuay5pbXBvcnRzLmZvckVhY2goKGlkKSA9PiB7XG4gICAgICAgICAgICAgICAgICAgIC8vIHRoaXMuaW5mbyhgcm91dGUgJHtyb3V0ZX0gaW1wb3J0aW5nICR7aWR9YCk7XG4gICAgICAgICAgICAgICAgICAgIGlmIChpZC5zdGFydHNXaXRoKFwiYXNzZXRzL1wiKSkge1xuICAgICAgICAgICAgICAgICAgICAgIGNvbnN0IHBhdGggPSBpZC5zcGxpdChcIi9cIikuc2xpY2UoMSk7XG4gICAgICAgICAgICAgICAgICAgICAgY29uc3QgaW1wb3J0ZWROZXN0aW5nID0gcGF0aC5sZW5ndGg7XG4gICAgICAgICAgICAgICAgICAgICAgY29uc3Qgcm91dGVOZXN0aW5nID0gcm91dGUuc3BsaXQoXCIvXCIpLmxlbmd0aDtcbiAgICAgICAgICAgICAgICAgICAgICBjb25zdCBmaWxlID0gcGF0aC5hdCgtMSk7XG4gICAgICAgICAgICAgICAgICAgICAgY29uc3QgY29ycmVjdGlvbiA9XG4gICAgICAgICAgICAgICAgICAgICAgICByb3V0ZU5lc3RpbmcgPiBpbXBvcnRlZE5lc3RpbmdcbiAgICAgICAgICAgICAgICAgICAgICAgICAgPyBuZXcgQXJyYXkocm91dGVOZXN0aW5nIC0gaW1wb3J0ZWROZXN0aW5nKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgLmZpbGwoXCIuLlwiKVxuICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgLmpvaW4oXCIvXCIpICsgXCIvXCJcbiAgICAgICAgICAgICAgICAgICAgICAgICAgOiBcIi4vXCI7XG5cbiAgICAgICAgICAgICAgICAgICAgICBpZiAoIWZpbGUpXG4gICAgICAgICAgICAgICAgICAgICAgICB0aGlzLmVycm9yKFxuICAgICAgICAgICAgICAgICAgICAgICAgICBgRmFpbGVkIHRvIGNvcnJlY3QgaW1wb3J0ICR7aWR9IGluIHJvdXRlICR7cm91dGV9YCxcbiAgICAgICAgICAgICAgICAgICAgICAgICk7XG5cbiAgICAgICAgICAgICAgICAgICAgICBjb25zdCBwYXR0ZXJuID0gbmV3IFJlZ0V4cChcbiAgICAgICAgICAgICAgICAgICAgICAgIGBpbXBvcnRcXFxccyooPzpcXFxcey4qXFxcXH1cXFxccypmcm9tXFxcXHMqKT8oXCIuLyR7ZmlsZX1cInwnLi8ke2ZpbGV9JylgLFxuICAgICAgICAgICAgICAgICAgICAgICk7XG4gICAgICAgICAgICAgICAgICAgICAgY2h1bmsuY29kZSA9IGNodW5rLmNvZGUucmVwbGFjZShcbiAgICAgICAgICAgICAgICAgICAgICAgIHBhdHRlcm4sXG4gICAgICAgICAgICAgICAgICAgICAgICAobWF0Y2gsIGlkOiBzdHJpbmcpID0+IHtcbiAgICAgICAgICAgICAgICAgICAgICAgICAgcmV0dXJuIG1hdGNoLnJlcGxhY2UoXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWQsXG4gICAgICAgICAgICAgICAgICAgICAgICAgICAgaWQucmVwbGFjZShcIi4vXCIsIGNvcnJlY3Rpb24pLFxuICAgICAgICAgICAgICAgICAgICAgICAgICApO1xuICAgICAgICAgICAgICAgICAgICAgICAgfSxcbiAgICAgICAgICAgICAgICAgICAgICApO1xuICAgICAgICAgICAgICAgICAgICB9XG4gICAgICAgICAgICAgICAgICB9KTtcblxuICAgICAgICAgICAgICAgICAgY2h1bmsuZmlsZU5hbWUgPSBgYXNzZXRzLyR7cm91dGV9LnN2ZWx0ZS5qc2A7XG4gICAgICAgICAgICAgICAgfVxuICAgICAgICAgICAgICB9XG4gICAgICAgICAgICB9KTtcbiAgICAgICAgICB9XG4gICAgICAgIH0sXG4gICAgICB9O1xuICAgIH0pKCksXG4gIF0sXG59KTtcbiJdLAogICJtYXBwaW5ncyI6ICI7QUFBaVUsU0FBUyxZQUFZLG9CQUFvQjtBQUMxVyxTQUFTLFVBQVUsU0FBUyxlQUFlO0FBQzNDLE9BQU8sZ0JBQWdCLHdCQUF3QjtBQUMvQyxTQUFTLFNBQVMsZUFBZSxrQkFBa0I7QUFDbkQsU0FBUyxvQkFBb0I7QUFFN0IsSUFBTyxzQkFBUSxhQUFhO0FBQUEsRUFDMUIsU0FBUztBQUFBLElBQ1AsWUFBWSxDQUFDLFFBQVE7QUFBQSxFQUN2QjtBQUFBLEVBQ0EsU0FBUztBQUFBLElBQ04sdUJBQU07QUFDTCxVQUFJLFFBQVE7QUFFWixhQUFPO0FBQUEsUUFDTCxNQUFNO0FBQUEsUUFDTixTQUFTO0FBQUEsUUFDVCxNQUFNLFVBQVUsUUFBUSxVQUFVLFNBQVM7QUFDekMsY0FBSTtBQUNKLGNBQUksV0FBVyxrQkFBa0I7QUFDL0IsZ0JBQUksQ0FBQyxTQUFVLE9BQU0sSUFBSSxNQUFNLDhCQUE4QjtBQUU3RCxpQkFBSyxVQUFVLDBCQUEwQixLQUFLLFFBQVEsT0FBTztBQUMzRCxvQkFBTSxJQUFJLE1BQU0sMENBQTBDO0FBRTVELG1CQUFPLHFCQUFxQixRQUFRLENBQUMsQ0FBQztBQUFBLFVBQ3hDO0FBRUEsY0FBSSxX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
