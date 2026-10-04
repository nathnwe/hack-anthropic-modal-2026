import { defineConfig } from "astro/config";

export default defineConfig({
  // Publishing can use another account or host without changing local URLs.
  site: process.env.REWIRE_SITE_URL || "https://nathnwe.github.io",
  base: process.env.REWIRE_BASE_PATH || "/hack-anthropic-modal-2026",
  output: "static",
  trailingSlash: "always",
  server: { port: 4322 },
  devToolbar: { enabled: false },
});
