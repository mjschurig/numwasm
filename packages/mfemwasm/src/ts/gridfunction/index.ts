// GridFunction class
export { GridFunction } from './GridFunction.js';

// Creation (standalone factory functions)
export { fromCoefficient } from './fromCoefficient.js';
export { fromFunction } from './fromFunction.js';

// Data Access (operations not on the class)
export { getValue, IntegrationPoint } from './getValue.js';
export { getGradient } from './getGradient.js';
export { getCurl } from './getCurl.js';
export { getDivergence } from './getDivergence.js';

// Projection (operations not on the class)
export { projectCoefficient } from './projectCoefficient.js';
export { projectFunction } from './projectFunction.js';
export { projectBdrCoefficient } from './projectBdrCoefficient.js';

// Norms
export { normL2 } from './normL2.js';
export { normLinf } from './normLinf.js';
export { normH1 } from './normH1.js';
export { computeError } from './computeError.js';

// Output
export { save } from './save.js';
export { toVTK, VTKExportOptions } from './toVTK.js';
