/**
 * LINPACK WASM Module Loader
 *
 * Handles loading and initialization of the LINPACK WebAssembly module.
 * LINPACK is a classic library for solving linear algebra problems.
 * Works in both Node.js and browser environments.
 */

import type { LINPACKModule, LINPACKModuleFactory, LINPACKModuleConfig } from './types.js';

let wasmModule: LINPACKModule | null = null;
let modulePromise: Promise<LINPACKModule> | null = null;

export interface LINPACKLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: LINPACKLoadConfig = {};

export function configureLINPACK(config: LINPACKLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure LINPACK WASM after module has started loading');
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

async function loadWasmGlue(wasmPath: string): Promise<LINPACKModuleFactory> {
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
  throw new Error('Unable to load LINPACK WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [
    nodePath.join(moduleDir, 'wasm', 'linpack.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'linpack.cjs'));
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
    `Cannot find linpack.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'pnpm run build:wasm' to compile the LINPACK WASM module.`
  );
}

async function loadWasmNode(): Promise<LINPACKModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('./wasm/linpack.wasm', import.meta.url).href;
  const glueUrl = new URL('./wasm/linpack.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<LINPACKModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('linpack.mjs', wasmUrlObj.href.replace(/linpack\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: LINPACKModuleFactory =
    glueModule.default || glueModule.createLINPACKModule;

  if (!createModule) {
    throw new Error('Failed to load LINPACK WASM glue code: no factory function found');
  }

  const config: LINPACKModuleConfig = {
    locateFile: (path: string) => {
      if (path.endsWith('.wasm')) {
        return wasmUrl!;
      }
      return path;
    },
  };
  const module = await createModule(config);

  return module;
}

export async function loadLINPACKModule(): Promise<LINPACKModule> {
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
        `Failed to load LINPACK WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getLINPACKModule(): LINPACKModule {
  if (!wasmModule) {
    throw new Error(
      'LINPACK WASM module not loaded. Call loadLINPACKModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isLINPACKLoaded(): boolean {
  return wasmModule !== null;
}

export function resetLINPACKModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
