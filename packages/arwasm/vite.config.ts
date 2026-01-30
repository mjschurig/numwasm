import { defineConfig } from 'vite';
import dts from 'vite-plugin-dts';
import { resolve } from 'path';

export default defineConfig({
  plugins: [
    dts({
      include: ['src/**/*.ts'],
      outDir: 'dist',
      rollupTypes: true,
    }),
  ],
  build: {
    lib: {
      entry: resolve(__dirname, 'src/index.ts'),
      name: 'ARWasm',
      formats: ['es', 'cjs'],
      fileName: (format) => `arwasm.${format === 'es' ? 'mjs' : 'cjs'}`,
    },
    outDir: 'dist',
    emptyOutDir: false,
    rollupOptions: {
      external: ['module', 'path', 'url', 'fs'],
    },
    target: ['es2020', 'node18'],
    minify: false,
    sourcemap: true,
  },
});
