/**
 * Factory registry to break circular dependencies between sparse matrix formats.
 * All sparse matrix modules register their constructors here, and other modules
 * use this to create instances without direct imports.
 */
import type { SparseConstructorArrays, COOConstructorArrays } from './types.js';
import type { CompressedSparseMatrix } from './compressed.js';
import type { SparseMatrix } from './base.js';

type CompressedSparseFactory = (
  arrays: SparseConstructorArrays,
  options: { shape: [number, number] }
) => CompressedSparseMatrix;

type COOSparseFactory = (
  arrays: COOConstructorArrays,
  options: { shape: [number, number] }
) => SparseMatrix;

let _csrFactory: CompressedSparseFactory | null = null;
let _cscFactory: CompressedSparseFactory | null = null;
let _cooFactory: COOSparseFactory | null = null;

export function registerCSRFactory(factory: CompressedSparseFactory): void {
  _csrFactory = factory;
}

export function registerCSCFactory(factory: CompressedSparseFactory): void {
  _cscFactory = factory;
}

export function registerCOOFactory(factory: COOSparseFactory): void {
  _cooFactory = factory;
}

export function createCSR(
  arrays: SparseConstructorArrays,
  shape: [number, number]
): CompressedSparseMatrix {
  if (!_csrFactory) throw new Error('CSRMatrix not registered');
  return _csrFactory(arrays, { shape });
}

export function createCSC(
  arrays: SparseConstructorArrays,
  shape: [number, number]
): CompressedSparseMatrix {
  if (!_cscFactory) throw new Error('CSCMatrix not registered');
  return _cscFactory(arrays, { shape });
}

export function createCOO(
  arrays: COOConstructorArrays,
  shape: [number, number]
): SparseMatrix {
  if (!_cooFactory) throw new Error('COOMatrix not registered');
  return _cooFactory(arrays, { shape });
}
