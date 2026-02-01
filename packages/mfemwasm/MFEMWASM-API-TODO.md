# MFEMWASM High-Level TypeScript API TODO

This document tracks the implementation status of all high-level TypeScript API functions for mfemwasm.

## Implementation Status Legend
- [ ] Not started
- [x] Implemented
- [~] Partially implemented (stub or incomplete)
- [D] Documented with TypeDoc

---

## 1. Mesh Operations (~35 functions)

### Mesh Class (mesh/Mesh.ts) [D]
- [x][D] `Mesh.create()` - Empty mesh
- [x][D] `Mesh.fromString(data)` - Load from MFEM format
- [x][D] `mesh.dimension` - Spatial dimension (1, 2, 3)
- [x][D] `mesh.spaceDimension` - Embedding dimension
- [x][D] `mesh.numElements` - Element count
- [x][D] `mesh.numVertices` - Vertex count
- [x][D] `mesh.numEdges` - Edge count
- [x][D] `mesh.numFaces` - Face count
- [x][D] `mesh.numBoundaryElements` - Boundary element count
- [x][D] `mesh.getAttribute(elemIdx)` - Element attribute
- [x][D] `mesh.getBdrAttribute(bdrIdx)` - Boundary attribute
- [x][D] `mesh.getBoundingBox()` - {min, max} coordinates
- [x][D] `mesh.refineUniform()` - Single uniform refinement
- [x][D] `mesh.isNURBS` - Whether mesh is NURBS
- [x][D] `mesh.nurbsOrder` - NURBS order (-1 if not NURBS)
- [x][D] `mesh.isNCMesh` - Whether NCMesh is enabled
- [x][D] `mesh.destroy()` - Free resources

### Mesh Creation Functions (mesh/*.ts)
- [x] `fromGmsh(data)` - Load from Gmsh format
- [x] `fromVTK(data)` - Load from VTK format
- [x] `makeCartesian1D(n, sx?)` - 1D Cartesian mesh
- [x][D] `makeCartesian2D(nx, ny, type?, sx?, sy?)` - 2D Cartesian mesh
- [x] `makeCartesian3D(nx, ny, nz, type?, ...)` - 3D Cartesian mesh
- [x] `makeCircle(npts, radius?)` - Circle mesh
- [x] `makeDisk(npts, radius?)` - Disk mesh
- [x] `makeSphere(npts, radius?)` - Sphere mesh
- [x] `makeCylinder(npts, length?, radius?)` - Cylinder mesh

### Mesh Property Functions (mesh/*.ts)
- [x] `getElementVolumes(mesh)` - Float64Array of volumes
- [x] `getElementCentroids(mesh)` - Array of centroid coordinates

### Mesh Refinement Functions (mesh/*.ts)
- [x][D] `refineUniform(mesh, times?)` - Multiple uniform refinements
- [x] `refineLocal(mesh, elements)` - Local adaptive refinement
- [x] `refineNonconforming(mesh, elements)` - Non-conforming refinement
- [x][D] `derefine(mesh, derefinements)` - Derefinement (coarsening)
- [x][D] `enableNCMesh(mesh)` - Enable non-conforming mesh support
- [x][D] `isNCMesh(mesh)` - Check NCMesh support
- [x][D] `getDerefTableSize(mesh)` - Get number of possible derefinements

### Mesh Transformation Functions (mesh/*.ts)
- [x][D] `transform(mesh, fn)` - Apply coordinate transformation
- [x] `scale(mesh, factor)` - Scale mesh
- [x] `translate(mesh, offset)` - Translate mesh
- [x] `rotate(mesh, angle, axis?)` - Rotate mesh
- [x] `setCurvature(mesh, order)` - Set high-order curvature

### Mesh Output Functions (mesh/*.ts)
- [x] `toMFEM(mesh)` - Export to MFEM format
- [x] `toVTK(mesh)` - Export to VTK format
- [x] `toGmsh(mesh)` - Export to Gmsh format

---

## 2. Finite Element Spaces (~15 functions)

### FiniteElementSpace Class (fespace/FiniteElementSpace.ts) [D]
- [x][D] `FiniteElementSpace.createH1(mesh, order, vdim?)` - H1 space
- [x][D] `fespace.ndofs` - Number of degrees of freedom
- [x][D] `fespace.vdim` - Vector dimension
- [x][D] `fespace.destroy()` - Free resources

### Space Creation Functions (fespace/*.ts) [D]
- [x][D] `createH1(mesh, order, vdim?)` - H1 space
- [x][D] `createH1Positive(mesh, order, vdim?)` - Bernstein basis
- [x][D] `createL2(mesh, order, vdim?)` - L2 space
- [x][D] `createDG(mesh, order, vdim?)` - Alias for L2
- [x][D] `createHcurl(mesh, order)` - H(curl) space
- [x][D] `createND(mesh, order)` - Nedelec alias
- [x][D] `createHdiv(mesh, order)` - H(div) space
- [x][D] `createRT(mesh, order)` - Raviart-Thomas alias
- [x][D] `createNURBS(mesh, order, vdim?)` - NURBS space (requires NURBS mesh)

### Space Property Functions (fespace/*.ts) [D]
- [x][D] `order(fespace)` - Polynomial order
- [x][D] `getDofMap(fespace, elemIdx)` - Element DoF indices
- [x][D] `getBoundaryDofs(fespace, bdrAttr)` - Boundary DoFs
- [x][D] `getEssentialDofs(fespace, bdrAttrs)` - Essential BC DoFs

---

## 3. Grid Functions (~20 functions)

### GridFunction Class (gridfunction/GridFunction.ts) [D]
- [x][D] `GridFunction.create(fespace)` - Create grid function
- [x][D] `gf.size` - Number of DOFs
- [x][D] `gf.getData()` - Get DOF values
- [x][D] `gf.setData(data)` - Set DOF values
- [x][D] `gf.projectConstant(value)` - Project constant
- [x][D] `gf.destroy()` - Free resources

### GridFunction Creation Functions (gridfunction/*.ts) [D]
- [x][D] `fromCoefficient(fespace, value)` - From constant coefficient
- [x][D] `fromFunction(fespace, fn)` - From JS function (limited)

### Data Access Functions (gridfunction/*.ts) [D]
- [x][D] `getValue(gf, elemIdx, ip)` - Evaluate at point
- [x][D] `getGradient(gf, elemIdx, ip)` - Gradient at point
- [x][D] `getCurl(gf, elemIdx, ip)` - Curl at point
- [x][D] `getDivergence(gf, elemIdx, ip)` - Divergence at point

### Projection Functions (gridfunction/*.ts) [D]
- [x][D] `projectCoefficient(gf, value)` - Project constant coefficient
- [x][D] `projectFunction(gf, fn)` - Project JS function (limited)
- [x][D] `projectBdrCoefficient(gf, value, bdrAttrs)` - Project on boundary

### Norm Functions (gridfunction/*.ts) [D]
- [x][D] `normL2(gf)` - L2 norm
- [x][D] `normLinf(gf)` - L∞ norm
- [x][D] `normH1(gf)` - H1 seminorm
- [x][D] `computeError(gf, exactValue)` - Error vs constant solution

### Output Functions (gridfunction/*.ts) [D]
- [x][D] `save(gf)` - Export to MFEM format string
- [x][D] `toVTK(gf, options?)` - Export to VTK format

---

## 4. Coefficients (~15 functions) - NOT YET IMPLEMENTED

### Constant Coefficients
- [ ] `ConstantCoefficient(value)`
- [ ] `VectorConstantCoefficient(values)`
- [ ] `MatrixConstantCoefficient(values)`

### Function Coefficients
- [ ] `FunctionCoefficient(fn)`
- [ ] `VectorFunctionCoefficient(fn)`
- [ ] `MatrixFunctionCoefficient(fn)`

### GridFunction Coefficients
- [ ] `GridFunctionCoefficient(gf)`
- [ ] `GradientCoefficient(gf)`
- [ ] `CurlCoefficient(gf)`
- [ ] `DivergenceCoefficient(gf)`

### Special Coefficients
- [ ] `PWConstCoefficient(attrs, values)`
- [ ] `VectorArrayCoefficient(coeffs)`
- [ ] `SumCoefficient(c1, c2)`
- [ ] `ProductCoefficient(c1, c2)`
- [ ] `PowerCoefficient(c, p)`

---

## 5-20. Additional Categories - NOT YET IMPLEMENTED

See original design document for complete list of:
- Bilinear Forms (~25 functions)
- Linear Forms (~15 functions)
- Boundary Conditions (~10 functions)
- Linear Solvers (~20 functions)
- Eigenvalue Solvers (~10 functions)
- Time Integration (~15 functions)
- Nonlinear Solvers (~10 functions)
- Parallel/MPI (~15 functions)
- Error Estimation (~10 functions)
- Visualization (~10 functions)
- Utility Functions (~15 functions)
- High-Level Problem Classes (~15 functions)
- Advanced Features (~15 functions)
- I/O Operations (~15 functions)
- Debugging/Diagnostics (~15 functions)

---

## Documentation Status

TypeDoc is configured in `typedoc.json`. Run `pnpm docs` to generate documentation.

### Fully Documented [D]
- `Mesh` class - Complete TypeDoc with examples
- `FiniteElementSpace` class - Complete TypeDoc with examples
- `GridFunction` class - Complete TypeDoc with examples
- `makeCartesian2D()` - Complete TypeDoc with examples
- `transform()` - Complete TypeDoc with examples
- `refineUniform()` - Complete TypeDoc with examples
- All fespace creation functions (createH1, createL2, createHcurl, etc.)
- All fespace property functions (order, getDofMap, etc.)
- All gridfunction functions (getValue, norms, save, toVTK, etc.)

### Needs Documentation
- Other mesh functions (makeCartesian1D, makeCartesian3D, etc.)

---

## Summary

| Category | Total | Implemented | Documented |
|----------|-------|-------------|------------|
| Mesh Class | 17 | 17 | 17 |
| Mesh Functions | 24 | 24 | 7 |
| FESpace Class | 4 | 4 | 4 |
| FESpace Functions | 13 | 13 | 13 |
| GridFunction Class | 6 | 6 | 6 |
| GridFunction Functions | 14 | 14 | 14 |
| **Total Core** | **78** | **78** | **61** |
