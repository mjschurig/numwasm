// FiniteElementSpace class
export { FiniteElementSpace } from './FiniteElementSpace.js';

// H1 (Continuous) Spaces
export { createH1 } from './createH1.js';
export { createH1Positive } from './createH1Positive.js';

// L2 (Discontinuous) Spaces
export { createL2 } from './createL2.js';
export { createDG } from './createDG.js';

// H(curl) Spaces (Nedelec)
export { createHcurl } from './createHcurl.js';
export { createND } from './createND.js';

// H(div) Spaces (Raviart-Thomas)
export { createHdiv } from './createHdiv.js';
export { createRT } from './createRT.js';

// NURBS Spaces
export { createNURBS } from './createNURBS.js';

// Space Property Functions (operations not on the class)
export { order } from './order.js';
export { getDofMap } from './getDofMap.js';
export { getBoundaryDofs } from './getBoundaryDofs.js';
export { getEssentialDofs } from './getEssentialDofs.js';
