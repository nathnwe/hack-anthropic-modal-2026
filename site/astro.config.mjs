import { defineConfig } from "astro/config";

export default defineConfig({
  site: "https://nathnwe.github.io",
  base: "/hack-anthropic-modal-2026",
  output: "static",
  trailingSlash: "always",
  server: { port: 4322 },
  devToolbar: { enabled: false },
});
