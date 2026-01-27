/**
 * Input/output for MATLAB .mat files and other formats.
 * @module io
 */

import { NotImplementedError } from '../errors.js';

/**
 * Load a MATLAB .mat file.
 * Mirrors scipy.io.loadmat.
 */
export function loadmat(
  _fileOrBuffer: string | ArrayBuffer,
  _options?: { squeeze_me?: boolean; struct_as_record?: boolean }
): Record<string, unknown> {
  throw new NotImplementedError('sciwasm.io.loadmat');
}

/**
 * Save data to a MATLAB .mat file.
 * Mirrors scipy.io.savemat.
 */
export function savemat(
  _filename: string,
  _mdict: Record<string, unknown>,
  _options?: { format?: '4' | '5' }
): void {
  throw new NotImplementedError('sciwasm.io.savemat');
}
