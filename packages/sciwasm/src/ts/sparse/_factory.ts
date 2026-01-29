/**
 * Factory registry to break circular dependencies between sparse matrix formats.
 * All sparse matrix modules register their constructors here, and other modules
 * use this to create instances without direct imports.
 */
import type {
  SparseConstructorArrays,
  COOConstructorArrays,
  LILConstructorArrays,
  DIAConstructorArrays,
  BSRConstructorArrays,
} from './types.js';
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

type LILSparseFactory = (
  arrays: LILConstructorArrays,
  options: { shape: [number, number] }
) => SparseMatrix;

type DOKSparseFactory = (
  dict: Map<string, number>,
  options: { shape: [number, number] }
) => SparseMatrix;

type DIASparseFactory = (
  arrays: DIAConstructorArrays,
  options: { shape: [number, number] }
) => SparseMatrix;

type BSRSparseFactory = (
  arrays: BSRConstructorArrays,
  options: { shape: [number, number] }
) => SparseMatrix;

let _csrFactory: CompressedSparseFactory | null = null;
let _cscFactory: CompressedSparseFactory | null = null;
let _cooFactory: COOSparseFactory | null = null;
let _lilFactory: LILSparseFactory | null = null;
let _dokFactory: DOKSparseFactory | null = null;
let _diaFactory: DIASparseFactory | null = null;
let _bsrFactory: BSRSparseFactory | null = null;

export function registerCSRFactory(factory: CompressedSparseFactory): void {
  _csrFactory = factory;
}

export function registerCSCFactory(factory: CompressedSparseFactory): void {
  _cscFactory = factory;
}

export function registerCOOFactory(factory: COOSparseFactory): void {
  _cooFactory = factory;
}

export function registerLILFactory(factory: LILSparseFactory): void {
  _lilFactory = factory;
}

export function registerDOKFactory(factory: DOKSparseFactory): void {
  _dokFactory = factory;
}

export function registerDIAFactory(factory: DIASparseFactory): void {
  _diaFactory = factory;
}

export function registerBSRFactory(factory: BSRSparseFactory): void {
  _bsrFactory = factory;
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

export function createLIL(
  arrays: LILConstructorArrays,
  shape: [number, number]
): SparseMatrix {
  if (!_lilFactory) throw new Error('LILMatrix not registered');
  return _lilFactory(arrays, { shape });
}

export function createDOK(
  dict: Map<string, number>,
  shape: [number, number]
): SparseMatrix {
  if (!_dokFactory) throw new Error('DOKMatrix not registered');
  return _dokFactory(dict, { shape });
}

export function createDIA(
  arrays: DIAConstructorArrays,
  shape: [number, number]
): SparseMatrix {
  if (!_diaFactory) throw new Error('DIAMatrix not registered');
  return _diaFactory(arrays, { shape });
}

export function createBSR(
  arrays: BSRConstructorArrays,
  shape: [number, number]
): SparseMatrix {
  if (!_bsrFactory) throw new Error('BSRMatrix not registered');
  return _bsrFactory(arrays, { shape });
}
