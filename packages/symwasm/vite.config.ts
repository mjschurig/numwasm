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
      name: 'SymWasm',
      formats: ['es', 'cjs'],
      fileName: (format) => `symwasm.${format === 'es' ? 'mjs' : 'cjs'}`,
    },
    outDir: 'dist',
    rollupOptions: {
      external: ['module', 'path', 'url', 'fs'],
    },
    target: ['es2020', 'node18'],
    minify: false,
    sourcemap: true,
  },
});
