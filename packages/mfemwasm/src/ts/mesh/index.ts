/**
 * Mesh module for mfemwasm.
 *
 * This module provides comprehensive mesh operations for finite element analysis,
 * including mesh creation, refinement, transformation, and I/O operations.
 *
 * ## Mesh Creation
 *
 * Create structured meshes on standard domains:
 * - {@link makeCartesian1D} - 1D line segment mesh
 * - {@link makeCartesian2D} - 2D rectangular mesh (quads or triangles)
 * - {@link makeCartesian3D} - 3D box mesh (hexes, tets, or wedges)
 *
 * Create meshes on curved domains:
 * - {@link makeCircle} - 2D circular disk mesh
 * - {@link makeDisk} - 2D disk mesh (elliptic mapping)
 * - {@link makeSphere} - 3D solid ball mesh
 * - {@link makeCylinder} - 3D solid cylinder mesh
 *
 * Load meshes from files:
 * - {@link fromGmsh} - Load Gmsh MSH format
 * - {@link fromVTK} - Load VTK legacy format
 *
 * ## Mesh Properties
 *
 * Query mesh information:
 * - {@link getElementVolumes} - Get volumes/areas/lengths of elements
 * - {@link getElementCentroids} - Get geometric centers of elements
 *
 * ## Mesh Refinement
 *
 * Refine meshes for improved accuracy:
 * - {@link refineUniform} - Globally refine all elements
 * - {@link refineLocal} - Refine selected elements (adaptive)
 * - {@link refineNonconforming} - Non-conforming local refinement
 *
 * Non-conforming mesh operations:
 * - {@link enableNCMesh} - Enable NCMesh support for derefinement
 * - {@link isNCMesh} - Check if mesh has NCMesh support
 * - {@link derefine} - Coarsen previously refined elements
 * - {@link getDerefTableSize} - Get number of possible derefinements
 *
 * ## Mesh Transformation
 *
 * Transform mesh coordinates:
 * - {@link transform} - Apply arbitrary coordinate transformation
 * - {@link scale} - Scale mesh uniformly or per-axis
 * - {@link translate} - Translate mesh by offset vector
 * - {@link rotate} - Rotate mesh around axis
 * - {@link setCurvature} - Set high-order geometry representation
 *
 * ## Mesh I/O
 *
 * Export meshes:
 * - {@link toMFEM} - Export to MFEM native format
 * - {@link toVTK} - Export to VTK format (for visualization)
 * - {@link toGmsh} - Export to Gmsh MSH format
 *
 * @example Basic workflow
 * ```typescript
 * import { makeCartesian2D, refineUniform, toVTK } from 'mfemwasm/mesh';
 *
 * // Create a mesh
 * const mesh = await makeCartesian2D(10, 10, 'tri');
 *
 * // Refine it
 * refineUniform(mesh, 2);
 *
 * // Export for visualization
 * const vtkStr = toVTK(mesh);
 *
 * // Clean up
 * mesh.destroy();
 * ```
 *
 * @example Adaptive refinement
 * ```typescript
 * import {
 *   makeCartesian2D,
 *   getElementCentroids,
 *   refineLocal,
 * } from 'mfemwasm/mesh';
 *
 * const mesh = await makeCartesian2D(8, 8);
 *
 * // Refine elements near a point
 * const centroids = getElementCentroids(mesh);
 * const nearOrigin = centroids
 *   .map((c, i) => ({ i, dist: Math.hypot(c[0] - 0.5, c[1] - 0.5) }))
 *   .filter(e => e.dist < 0.2)
 *   .map(e => e.i);
 *
 * refineLocal(mesh, nearOrigin);
 * mesh.destroy();
 * ```
 *
 * @module mesh
 */

// Mesh class
export { Mesh, BoundingBox } from './Mesh.js';

// Mesh Creation (standalone factory functions)
export { fromGmsh } from './fromGmsh.js';
export { fromVTK } from './fromVTK.js';
export { makeCartesian1D } from './makeCartesian1D.js';
export { makeCartesian2D, ElementType2D } from './makeCartesian2D.js';
export { makeCartesian3D, ElementType3D } from './makeCartesian3D.js';
export { makeCircle } from './makeCircle.js';
export { makeDisk } from './makeDisk.js';
export { makeSphere } from './makeSphere.js';
export { makeCylinder } from './makeCylinder.js';

// Mesh Property Functions (operations not on the class)
export { getElementVolumes } from './getElementVolumes.js';
export { getElementCentroids } from './getElementCentroids.js';

// Mesh Refinement
export { refineUniform } from './refineUniform.js';
export { refineLocal } from './refineLocal.js';
export { refineNonconforming } from './refineNonconforming.js';

// Non-conforming Mesh (NCMesh) and Derefinement
export { derefine, enableNCMesh, isNCMesh, getDerefTableSize } from './derefine.js';

// Mesh Transformation
export { transform } from './transform.js';
export { scale } from './scale.js';
export { translate } from './translate.js';
export { rotate } from './rotate.js';
export { setCurvature } from './setCurvature.js';

// Mesh Output
export { toMFEM } from './toMFEM.js';
export { toVTK } from './toVTK.js';
export { toGmsh } from './toGmsh.js';
