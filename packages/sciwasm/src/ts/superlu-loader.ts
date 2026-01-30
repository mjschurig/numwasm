/**
 * SuperLU WASM Module Loader
 *
 * Handles loading and initialization of the SuperLU WebAssembly module.
 * SuperLU is a library for solving sparse systems of linear equations.
 * Works in both Node.js and browser environments.
 */

import type { SuperLUModule, SuperLUModuleFactory } from './superlu-types.js';

let wasmModule: SuperLUModule | null = null;
let modulePromise: Promise<SuperLUModule> | null = null;

export interface SuperLULoadConfig {
  wasmUrl?: string;
}

let wasmConfig: SuperLULoadConfig = {};

export function configureSuperLU(config: SuperLULoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure SuperLU WASM after module has started loading');
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

async function loadWasmGlue(wasmPath: string): Promise<SuperLUModuleFactory> {
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
  throw new Error('Unable to load SuperLU WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [
    nodePath.join(moduleDir, 'wasm', 'superlu.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'superlu.cjs'));
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
    `Cannot find superlu.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'bash scripts/build-superlu.sh' to compile the SuperLU WASM module.`
  );
}

async function loadWasmNode(): Promise<SuperLUModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('./wasm/superlu.wasm', import.meta.url).href;
  const glueUrl = new URL('./wasm/superlu.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<SuperLUModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('superlu.mjs', wasmUrlObj.href.replace(/superlu\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: SuperLUModuleFactory =
    glueModule.default || glueModule.createSuperLUModule;

  if (!createModule) {
    throw new Error('Failed to load SuperLU WASM glue code: no factory function found');
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

export async function loadSuperLUModule(): Promise<SuperLUModule> {
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
        `Failed to load SuperLU WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getSuperLUModule(): SuperLUModule {
  if (!wasmModule) {
    throw new Error(
      'SuperLU WASM module not loaded. Call loadSuperLUModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isSuperLULoaded(): boolean {
  return wasmModule !== null;
}

export function resetSuperLUModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
