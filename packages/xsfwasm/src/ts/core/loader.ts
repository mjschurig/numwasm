/**
 * XSF (Special Functions) WASM Module Loader
 *
 * Handles loading and initialization of the XSF WebAssembly module.
 * Works in both Node.js and browser environments.
 */

import type { XSFModule, XSFModuleFactory } from './types.js';

let wasmModule: XSFModule | null = null;
let modulePromise: Promise<XSFModule> | null = null;

export interface XSFLoadConfig {
  wasmUrl?: string;
}

let wasmConfig: XSFLoadConfig = {};

export function configureXSF(config: XSFLoadConfig): void {
  if (wasmModule || modulePromise) {
    throw new Error('Cannot configure XSF WASM after module has started loading');
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

async function loadWasmGlue(wasmPath: string): Promise<XSFModuleFactory> {
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
  throw new Error('Unable to load XSF WASM module: no suitable require function available');
}

async function findWasmCjs(): Promise<string> {
  const nodePath = await import('path');
  const nodeFs = await import('fs');
  const moduleDir = await getModuleDir();

  const candidates = [
    nodePath.join(moduleDir, '..', 'wasm', 'xsf.cjs'),
  ];

  // Walk up to find package root (directory containing package.json)
  let dir = moduleDir;
  for (let i = 0; i < 10; i++) {
    const pkgPath = nodePath.join(dir, 'package.json');
    if (nodeFs.existsSync(pkgPath)) {
      candidates.push(nodePath.join(dir, 'dist', 'wasm', 'xsf.cjs'));
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
    `Cannot find xsf.cjs. Searched:\n${candidates.map(c => '  - ' + c).join('\n')}\nRun 'pnpm run build:wasm' to compile the XSF WASM module.`
  );
}

async function loadWasmNode(): Promise<XSFModule> {
  const wasmPath = await findWasmCjs();
  const createModule = await loadWasmGlue(wasmPath);
  return await createModule();
}

function getBrowserAssetUrls(): { wasmUrl: string; glueUrl: string } {
  const wasmUrl = new URL('../wasm/xsf.wasm', import.meta.url).href;
  const glueUrl = new URL('../wasm/xsf.mjs', import.meta.url).href;
  return { wasmUrl, glueUrl };
}

async function loadWasmBrowser(): Promise<XSFModule> {
  let wasmUrl: string;
  let glueUrl: string;

  if (wasmConfig.wasmUrl) {
    wasmUrl = wasmConfig.wasmUrl;
    const wasmUrlObj = new URL(wasmUrl);
    glueUrl = new URL('xsf.mjs', wasmUrlObj.href.replace(/xsf\.wasm$/, '')).href;
  } else {
    const urls = getBrowserAssetUrls();
    wasmUrl = urls.wasmUrl;
    glueUrl = urls.glueUrl;
  }

  const glueModule = await import(/* @vite-ignore */ glueUrl);
  const createModule: XSFModuleFactory =
    glueModule.default || glueModule.createXSFModule;

  if (!createModule) {
    throw new Error('Failed to load XSF WASM glue code: no factory function found');
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

export async function loadXSFModule(): Promise<XSFModule> {
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
        `Failed to load XSF WASM module: ${error instanceof Error ? error.message : String(error)}`
      );
    }
  })();

  return modulePromise;
}

export function getXSFModule(): XSFModule {
  if (!wasmModule) {
    throw new Error(
      'XSF WASM module not loaded. Call loadXSFModule() first and await the result.'
    );
  }
  return wasmModule;
}

export function isXSFLoaded(): boolean {
  return wasmModule !== null;
}

export function resetXSFModule(): void {
  wasmModule = null;
  modulePromise = null;
  wasmConfig = {};
}
