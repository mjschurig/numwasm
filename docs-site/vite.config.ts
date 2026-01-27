import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/vite';

export default defineConfig({
  plugins: [react(), tailwindcss()],
  base: '/',
  build: {
    outDir: 'dist',
  },
  ssr: {
    // Bundle these CJS packages into the SSR output to avoid ESM import issues
    noExternal: ['react-helmet-async'],
  },
});
