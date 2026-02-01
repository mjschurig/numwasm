import { defineConfig } from 'vitest/config';
import { resolve } from 'path';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/ts/**/*.test.ts'],
    testTimeout: 60000,
    hookTimeout: 30000,
  },
  resolve: {
    alias: {
      arwasm: resolve(__dirname, 'dist/arwasm.mjs'),
    },
  },
});
