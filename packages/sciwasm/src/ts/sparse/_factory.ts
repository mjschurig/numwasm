/**
 * Factory registry to break circular dependencies between CSR and CSC.
 * Both modules register their constructors here, and other modules
 * use this to create instances without direct imports.
 */
import type { SparseConstructorArrays } from './types.js';
import type { CompressedSparseMatrix } from './compressed.js';

type SparseFactory = (
  arrays: SparseConstructorArrays,
  options: { shape: [number, number] }
) => CompressedSparseMatrix;

let _csrFactory: SparseFactory | null = null;
let _cscFactory: SparseFactory | null = null;

export function registerCSRFactory(factory: SparseFactory): void {
  _csrFactory = factory;
}

export function registerCSCFactory(factory: SparseFactory): void {
  _cscFactory = factory;
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
