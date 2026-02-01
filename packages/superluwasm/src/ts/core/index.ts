/**
 * Core module - Module loading, configuration, memory management, and error handling
 *
 * @module core
 */

export {
  loadSuperLUModule,
  getSuperLUModule,
  isSuperLULoaded,
  isSuperLULoading,
  resetSuperLUModule,
  configureSuperLU,
  getSuperLUConfig,
  type SuperLULoadConfig,
} from './loader.js';

export {
  // Memory allocation
  allocateDoubles,
  allocateFloats,
  allocateInts,
  allocateValues,

  // Memory reading
  readDoubles,
  readFloats,
  readInts,
  readInt,
  readDouble,
  readByte,
  readValues,

  // Memory writing
  writeDoubles,
  writeFloats,
  writeInts,
  writeInt,
  writeDouble,
  writeByte,

  // Memory management
  freeAll,
  withAllocations,

  // Conversions
  toFloat64Array,
  toFloat32Array,
  toInt32Array,
  getElementSize,
  isComplexType,

  // Error handling
  getSuperLUErrorMessage,
  checkSuperLUError,
  SUPERLU_ERRORS,

  // Constants
  DOUBLE_SIZE,
  FLOAT_SIZE,
  INT_SIZE,
  OPTIONS_STRUCT_SIZE,
  STAT_STRUCT_SIZE,
  SUPERMATRIX_SIZE,
  GLOBALLU_SIZE,
  MEM_USAGE_SIZE,
} from './helpers.js';
