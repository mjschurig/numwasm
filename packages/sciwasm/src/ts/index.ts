/**
 * sciwasm - SciPy-inspired scientific computing in TypeScript/WebAssembly.
 *
 * Built on top of numwasm, sciwasm provides advanced scientific computing
 * functions including optimization, integration, interpolation, signal
 * processing, linear algebra extensions, statistics, and more.
 */

export * as optimize from './optimize/index.js';
export * as integrate from './integrate/index.js';
export * as interpolate from './interpolate/index.js';
export * as stats from './stats/index.js';
export * as signal from './signal/index.js';
export * as spatial from './spatial/index.js';
export * as special from './special/index.js';
export * as sparse from './sparse/index.js';
export * as cluster from './cluster/index.js';
export * as io from './io/index.js';
export * as ndimage from './ndimage/index.js';
export * as constants from './constants/index.js';
export { NotImplementedError } from './errors.js';
