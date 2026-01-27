import { defineConfig } from 'vitest/config';
import { resolve } from 'path';

export default defineConfig({
  test: {
    globals: true,
    environment: 'node',
    include: ['tests/ts/**/*.test.ts'],
    testTimeout: 30000,
    hookTimeout: 30000,
  },
  resolve: {
    alias: {
      // Map package name to dist for tests
      'numwasm': resolve(__dirname, 'dist/numjs.mjs'),
    },
  },
});
