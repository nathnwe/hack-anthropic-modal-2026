// @ts-check
import { defineConfig } from 'astro/config';

import tailwindcss from '@tailwindcss/vite';

// https://astro.build/config
// `site` + `base` are required for GitHub Pages project sites: the repo is served
// from https://<owner>.github.io/hack-anthropic-modal-2026/, not from the domain root.
export default defineConfig({
  site: 'https://nathnwe.github.io',
  base: '/hack-anthropic-modal-2026',
  vite: {
    plugins: [tailwindcss()]
  }
});
