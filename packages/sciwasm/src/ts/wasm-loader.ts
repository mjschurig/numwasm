/**
 * WASM Module Loader for sciwasm
 *
 * Handles loading and initialization of the SciWASM WebAssembly module.
 * Works in both Node.js and browser environments.
 * Pattern follows numwasm's wasm-loader.ts.
 */

import type { SciWasmModule, WasmModuleFactory } from './wasm-types.js';

let wasmModule: SciWasmModule | null = null;
let modulePromise: Promise<SciWasmModule> | null = null;

export interface WasmLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: WasmLoadConfig = {};

export function configureWasm(config: WasmLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure WASM after module has started loading');
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

async function loadWasmGlue(wasmPath: string): Promise<WasmModuleFactory> {
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
  throw new Error('Unable to load WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  // Try multiple candidate paths:
  // 1. <moduleDir>/wasm/sciwasm.cjs  (when running from dist/)
  // 2. <packageRoot>/dist/wasm/sciwasm.cjs  (when running from src/ts/ during dev/test)
  const candidates = [
    nodePath.join(moduleDir, 'wasm', 'sciwasm.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'sciwasm.cjs'));
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
    `Cannot find sciwasm.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'bash scripts/build-wasm.sh' to compile the WASM module.`
  );
}

async function loadWasmNode(): Promise<SciWasmModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('./wasm/sciwasm.wasm', import.meta.url).href;
  const glueUrl = new URL('./wasm/sciwasm.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<SciWasmModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('sciwasm.mjs', wasmUrlObj.href.replace(/sciwasm\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: WasmModuleFactory =
    glueModule.default || glueModule.createSciWASMModule;

  if (!createModule) {
    throw new Error('Failed to load WASM glue code: no factory function found');
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

export async function loadWasmModule(): Promise<SciWasmModule> {
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
        `Failed to load WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getWasmModule(): SciWasmModule {
  if (!wasmModule) {
    throw new Error(
      'WASM module not loaded. Call loadWasmModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isWasmLoaded(): boolean {
  return wasmModule !== null;
}

export function resetWasmModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
