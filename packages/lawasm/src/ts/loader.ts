/**
 * LAPACK WASM Module Loader
 *
 * Handles loading and initialization of the LAPACK WebAssembly module.
 * LAPACK (Linear Algebra PACKage) provides routines for solving linear systems,
 * eigenvalue problems, singular value decomposition, and matrix factorizations.
 *
 * This module bundles BLAS internally, so no separate BLAS module is needed.
 * Works in both Node.js and browser environments.
 */

import type { LAPACKModule, LAPACKModuleFactory } from './types.js';

let wasmModule: LAPACKModule | null = null;
let modulePromise: Promise<LAPACKModule> | null = null;

export interface LAPACKLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: LAPACKLoadConfig = {};

export function configureLAPACK(config: LAPACKLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure LAPACK WASM after module has started loading');
  }
  wasmConfig = { ...config };
}

const isNode =
  typeof process !== 'undefined' &&
  process.versions != null &&
  process.versions.node != null;

const isBrowser =
  typeof window !== 'undefined' ||
  (typeof self !== 'undefined' &&
    typeof (self as unknown as { importScripts?: unknown }).importScripts ===
      'function');

async function getModuleDir(): Promise<string> {
  const nodePath = await import('path');
  const nodeUrl = await import('url');
  try {
    if (import.meta && import.meta.url) {
      return nodePath.dirname(nodeUrl.fileURLToPath(import.meta.url));
    }
  } catch {
    // fall through
  }
  if (typeof __dirname !== 'undefined') {
    return __dirname;
  }
  return process.cwd();
}

async function loadWasmGlue(wasmPath: string): Promise<LAPACKModuleFactory> {
  try {
    if (import.meta && import.meta.url) {
      const nodeModule = await import('module');
      const req = nodeModule.createRequire(import.meta.url);
      return req(wasmPath);
    }
  } catch {
    // fall through
  }
  if (typeof require !== 'undefined') {
    return require(wasmPath);
  }
  throw new Error('Unable to load LAPACK WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [
    nodePath.join(moduleDir, 'wasm', 'lapack.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'lapack.cjs'));
      break;
    }
    const parent = nodePath.dirname(dir);
    if (parent === dir) break;
    dir = parent;
  }

  for (const candidate of candidates) {
    if (nodeFs.existsSync(candidate)) {
      return candidate;
    }
  }

  throw new Error(
    `Cannot find lapack.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'pnpm run build:wasm' to compile the LAPACK WASM module.`
  );
}

async function loadWasmNode(): Promise<LAPACKModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL(/* @vite-ignore */ './wasm/lapack.wasm', import.meta.url).href;
  const glueUrl = new URL(/* @vite-ignore */ './wasm/lapack.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<LAPACKModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('lapack.mjs', wasmUrlObj.href.replace(/lapack\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: LAPACKModuleFactory =
    glueModule.default || glueModule.createLAPACKModule;

  if (!createModule) {
    throw new Error('Failed to load LAPACK WASM glue code: no factory function found');
  }

  const module = await createModule({
    locateFile: (path: string) => {
      if (path.endsWith('.wasm')) {
        return wasmUrl!;
      }
      return path;
    },
  });

  return module;
}

export async function loadLAPACKModule(): Promise<LAPACKModule> {
  if (wasmModule) return wasmModule;
  if (modulePromise) return modulePromise;

  modulePromise = (async () => {
    try {
      if (isNode) {
        wasmModule = await loadWasmNode();
      } else if (isBrowser) {
        wasmModule = await loadWasmBrowser();
      } else {
        throw new Error('Unknown environment: neither Node.js nor browser detected');
      }
      return wasmModule;
    } catch (error) {
      modulePromise = null;
      throw new Error(
        `Failed to load LAPACK WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getLAPACKModule(): LAPACKModule {
  if (!wasmModule) {
    throw new Error(
      'LAPACK WASM module not loaded. Call loadLAPACKModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isLAPACKLoaded(): boolean {
  return wasmModule !== null;
}

export function resetLAPACKModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
