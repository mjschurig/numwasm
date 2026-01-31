/**
 * WASM Module Loader
 *
 * Re-exports from _core/wasm-loader.js for backwards compatibility.
 */

export {
  loadWasmModule,
  isWasmLoaded,
  getWasmModule,
  configureWasm,
  resetWasmModule,
} from "./_core/wasm-loader.js";
export type { WasmLoadConfig } from "./_core/wasm-loader.js";
