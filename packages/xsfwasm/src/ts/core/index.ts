/**
 * Core module - loader, types, and utilities
 */

export {
  loadXSFModule,
  getXSFModule,
  isXSFLoaded,
  resetXSFModule,
  configureXSF,
  type XSFLoadConfig,
} from './loader.js';

export type { XSFModule, XSFModuleFactory, XSFModuleOptions } from './types.js';

export { ensureLoaded, toFloat64Array, type ArrayInput } from './utils.js';
