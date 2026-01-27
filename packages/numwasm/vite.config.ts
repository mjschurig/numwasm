import { defineConfig } from 'vite';
import dts from 'vite-plugin-dts';
import { resolve } from 'path';

export default defineConfig({
  plugins: [
    dts({
      include: ['src/ts/**/*.ts'],
      outDir: 'dist',
      rollupTypes: true,
    }),
  ],
  build: {
    lib: {
      entry: resolve(__dirname, 'src/ts/index.ts'),
      name: 'NumJS',
      formats: ['es', 'cjs'],
      fileName: (format) => `numjs.${format === 'es' ? 'mjs' : 'cjs'}`,
    },
    outDir: 'dist',
    emptyOutDir: false, // Don't delete WASM files
    rollupOptions: {
      external: ['module', 'path', 'url', 'fs'],
      output: {
        // Ensure WASM is referenced correctly
        assetFileNames: 'wasm/[name][extname]',
      },
    },
    target: ['es2020', 'node18'],
    minify: false, // Keep readable for debugging
    sourcemap: true,
  },
  // Handle WASM files
  assetsInclude: ['**/*.wasm'],
});
