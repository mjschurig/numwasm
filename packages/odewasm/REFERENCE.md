# odewasm Function Reference & Phase Plan

High-quality ODE solvers compiled to WebAssembly, including Ernst Hairer's collection and classic Netlib solvers.

---

## Available Reference Solvers

This package includes reference implementations from multiple sources. The table below shows all available solvers and their implementation status.

### Hairer Non-Stiff Solvers (`reference/nonstiff/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **DOPRI5** | Dormand-Prince 5(4) explicit RK | General-purpose, non-stiff | ✅ Implemented |
| **DOP853** | Dormand-Prince 8(5,3) explicit RK | High accuracy, smooth problems | ✅ Implemented |
| **ODEX** | GBS extrapolation with midpoint rule | Variable order, non-stiff | ⬚ Not started |
| **ODEX2** | GBS extrapolation (variant) | Variable order, non-stiff | ⬚ Not started |
| **RETARD** | DOPRI5 for delay equations | Delay differential equations (DDEs) | ⬚ Not started |

### Hairer Stiff Solvers (`reference/stiff/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **RADAU5** | Implicit Radau IIA order 5 | Stiff systems | ✅ Implemented |
| **RADAU** | Implicit Radau IIA variable order (5,9,13) | Stiff systems, high accuracy | ⬚ Not started |
| **RADAUP** | RADAU with preconditioning | Large stiff systems | ⬚ Not started |
| **RODAS** | Rosenbrock order 3(4) | Stiff systems, no Newton iteration | ⬚ Not started |
| **SEULEX** | Extrapolation with implicit Euler | Stiff systems, variable order | ⬚ Not started |

### Delay Differential Equations (`reference/delay/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **RADAR5** | Radau IIA for stiff DDEs | Stiff delay differential equations | ⬚ Not started |

### Method of Lines / Parabolic PDEs (`reference/mol/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **ROCK2** | Stabilized explicit RK order 2 | Parabolic PDEs, mildly stiff | ⬚ Not started |
| **ROCK4** | Stabilized explicit RK order 4 | Parabolic PDEs, mildly stiff | ⬚ Not started |

### Mechanical Systems (`reference/mechanical/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **HEM5** | Half-explicit method order 5 | Constrained mechanical systems (index-2 DAE) | ⬚ Not started |
| **PHEM56** | Partitioned half-explicit | Mechanical systems | ⬚ Not started |

### Geometric Integration (`reference/gnicodes/`)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **GNI_IRK2** | Symplectic Gauss IRK | Hamiltonian systems, energy preservation | ⬚ Not started |

---

## Netlib Solvers (`reference/netlib/`)

Classic ODE/DAE solvers from the Netlib repository.

### Initial Value Problem Solvers

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **RKF45** | Runge-Kutta-Fehlberg 4(5) | Non-stiff, moderate accuracy | ⬚ Not started |
| **DVERK** | Verner 6(5) RK | Non-stiff, higher accuracy | ⬚ Not started |
| **ODE** | Adams-Bashforth-Moulton | Non-stiff, multistep | ⬚ Not started |
| **VODE** | Variable-coefficient ODE solver (BDF/Adams) | Stiff or non-stiff, ODEPACK family | ⬚ Not started |
| **ZVODE** | VODE for complex-valued ODEs | Complex stiff/non-stiff | ⬚ Not started |
| **VODPK** | VODE with Krylov methods | Large stiff systems | ⬚ Not started |
| **EPSODE** | Episode package | Stiff systems | ⬚ Not started |
| **RKSUITE** | RK suite (orders 2,3 / 4,5 / 7,8) | Non-stiff with error assessment | ⬚ Not started |

### Stabilized Explicit Methods

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **RKC** | Runge-Kutta-Chebyshev order 2 | Mildly stiff, parabolic PDEs | ⬚ Not started |
| **IRKC** | Implicit RKC (F90) | Stiff parabolic PDEs | ⬚ Not started |

### Differential-Algebraic Equations (DAEs)

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **DDASSL** | BDF for DAEs | DAE index 0-1, G(t,y,y')=0 | ⬚ Not started |
| **DDASRT** | DDASSL with root finding | DAEs with event detection | ⬚ Not started |
| **MEBDFI** | Modified extended BDF | DAE index up to 3 | ⬚ Not started |
| **MEBDFSO** | MEBDF with sparse output | Large DAE systems | ⬚ Not started |
| **MEBDFDAE** | MEBDF for DAEs | DAE systems | ⬚ Not started |

### Boundary Value Problem Solvers

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **COLNEW** | Collocation with B-splines | BVP for mixed-order ODEs | ⬚ Not started |
| **COLSYS** | Collocation system | BVP for ODEs | ⬚ Not started |
| **COLDAE** | Collocation for DAEs | BVP for DAE systems | ⬚ Not started |
| **COLMOD** | Modified collocation | BVP with singular points | ⬚ Not started |
| **TWPBVP** | Two-point BVP solver | General two-point BVPs | ⬚ Not started |
| **MIRKDC** | Mono-implicit RK deferred correction | BVP with defect control | ⬚ Not started |
| **ACDC** | Automatic continuation + DC | BVPs with parameter continuation | ⬚ Not started |

### Special Purpose

| Solver | Method | Use Case | Status |
|--------|--------|----------|--------|
| **MUS** | Multiple shooting | BVPs via shooting method | ⬚ Not started |
| **DRESOL** | Direct residual | Large sparse systems | ⬚ Not started |
| **YALE** | Yale sparse solver integration | Systems with Yale sparse format | ⬚ Not started |

### Archived Packages (tarballs)

| Package | Contents | Status |
|---------|----------|--------|
| **CVODE** | LLNL C-language VODE | ⬚ Not started |
| **DASPK** | Improved DDASSL with Krylov | ⬚ Not started |
| **DASKR** | DASPK with root finding | ⬚ Not started |
| **DGELDA** | DAE solver with sensitivity | ⬚ Not started |
| **PARSODES** | Parallel stiff ODE solvers | ⬚ Not started |

---

## Implementation Status

### Phase 1: Core Infrastructure (Foundation)

**Dependencies:** None

- [x] `loadODEModule` — TS — Load WASM module
- [x] `getODEModule` — TS — Get loaded module
- [x] `isODELoaded` — TS — Check load status
- [x] `resetODEModule` — TS — Reset module
- [x] `configureODE` — TS — Pre-load configuration
- [x] `_malloc` — WASM — Memory allocation (Emscripten)
- [x] `_free` — WASM — Memory deallocation (Emscripten)
- [x] `_wasm_malloc` — WASM — WASM-specific malloc
- [x] `_wasm_free` — WASM — WASM-specific free
- [x] `_wasm_set_fcn_callback` — WASM — Register ODE function
- [x] `_wasm_dopri5_work_size` — WASM — DOPRI5 work array size
- [x] `_wasm_dopri5_iwork_size` — WASM — DOPRI5 iwork array size

### Phase 2: Explicit Solvers (Non-stiff)

**Dependencies:** Phase 1

- [x] `_wasm_dopri5` — WASM — DOPRI5 solver
- [x] `_wasm_dop853` — WASM — DOP853 solver
- [x] `_wasm_dop853_work_size` — WASM — DOP853 work array size
- [x] `_wasm_dop853_iwork_size` — WASM — DOP853 iwork array size
- [x] `_wasm_set_solout_callback` — WASM — Solution output callback
- [x] `solveExplicit` — TS — High-level explicit solver wrapper
- [x] `ExplicitMethod` — TS — Method enum (RK45, DOP853)

### Phase 3: Dense Output (Interpolation)

**Dependencies:** Phase 2

- [x] `_wasm_contd5` — WASM — DOPRI5 interpolation
- [x] `_wasm_contd8` — WASM — DOP853 interpolation
- [x] `DenseOutput` — TS — Interpolation interface (type defined)
- [ ] Dense output in `solveExplicit` — TS — Actual interpolation in TS wrapper (partial: t_eval uses endpoint approximation)

### Phase 4: Implicit Solvers (Stiff)

**Dependencies:** Phase 1

- [x] `_wasm_radau5` — WASM — RADAU5 solver
- [x] `_wasm_set_jac_callback` — WASM — Jacobian callback
- [x] `_wasm_radau5_work_size` — WASM — RADAU5 work array size
- [x] `_wasm_radau5_iwork_size` — WASM — RADAU5 iwork array size
- [x] `_wasm_contr5` — WASM — RADAU5 interpolation
- [x] `solveRadau` — TS — High-level Radau wrapper
- [x] `ImplicitMethod` — TS — Method enum (Radau)
- [x] `JacobianFunction` — TS — Jacobian type

### Phase 5: Unified High-Level API

**Dependencies:** Phase 2, Phase 3, Phase 4

- [x] `solve_ivp` — TS — SciPy-compatible unified interface
- [x] `ODEFunction` — TS — ODE function type
- [x] `SolveIVPOptions` — TS — Options interface
- [x] `SolveIVPResult` — TS — Result interface
- [x] `STATUS_MESSAGES` — TS — Status code messages
- [x] `ODEMethod` — TS — Union of all methods
- [ ] `dense_output` in result — TS — Return DenseOutput interpolator from solve_ivp
- [ ] `t_eval` interpolation — TS — Proper dense output interpolation for t_eval points

---

## Future Phases (Not Yet Implemented)

### Phase 6: ODEX (Extrapolation Method)

**Dependencies:** Phase 1
**Source:** `reference/nonstiff/nonstiff/odex.f`

- [ ] `_wasm_odex` — WASM — ODEX solver (GBS extrapolation)
- [ ] `_wasm_odex_work_size` — WASM — Work array size
- [ ] `_wasm_odex_iwork_size` — WASM — Integer work array size
- [ ] `_wasm_contex` — WASM — Dense output interpolation
- [ ] `solveOdex` — TS — High-level wrapper
- [ ] Add `ODEX` to `ExplicitMethod` enum

### Phase 7: RODAS (Rosenbrock Method)

**Dependencies:** Phase 1
**Source:** `reference/stiff/stiff/rodas.f`

- [ ] `_wasm_rodas` — WASM — RODAS solver
- [ ] `_wasm_rodas_work_size` — WASM — Work array size
- [ ] `_wasm_rodas_iwork_size` — WASM — Integer work array size
- [ ] `_wasm_contro` — WASM — Dense output interpolation
- [ ] `_wasm_set_dfx_callback` — WASM — df/dx callback
- [ ] `solveRodas` — TS — High-level wrapper
- [ ] Add `Rodas` to `ImplicitMethod` enum

### Phase 8: SEULEX (Extrapolation for Stiff)

**Dependencies:** Phase 1
**Source:** `reference/stiff/stiff/seulex.f`

- [ ] `_wasm_seulex` — WASM — SEULEX solver
- [ ] `_wasm_seulex_work_size` — WASM — Work array size
- [ ] `_wasm_seulex_iwork_size` — WASM — Integer work array size
- [ ] `_wasm_contse` — WASM — Dense output interpolation
- [ ] `solveSeulex` — TS — High-level wrapper
- [ ] Add `Seulex` to `ImplicitMethod` enum

### Phase 9: RADAU (Variable Order Stiff)

**Dependencies:** Phase 1
**Source:** `reference/stiff/stiff/radau.f`

- [ ] `_wasm_radau` — WASM — RADAU solver (variable order 5,9,13)
- [ ] `_wasm_radau_work_size` — WASM — Work array size
- [ ] `_wasm_radau_iwork_size` — WASM — Integer work array size
- [ ] `_wasm_contra` — WASM — Dense output interpolation
- [ ] `solveRadauVO` — TS — High-level wrapper (variable order variant)

### Phase 10: RETARD (Delay Differential Equations)

**Dependencies:** Phase 2
**Source:** `reference/nonstiff/nonstiff/retard.f`

- [ ] `_wasm_retard` — WASM — RETARD solver for DDEs
- [ ] `_wasm_retard_work_size` — WASM — Work array size
- [ ] `_wasm_set_phi_callback` — WASM — Initial history function
- [ ] `_wasm_ylag` — WASM — Get lagged solution values
- [ ] `DDEFunction` — TS — DDE function type (accesses lagged values)
- [ ] `HistoryFunction` — TS — Initial history φ(t) for t < t0
- [ ] `solveDDE` — TS — High-level DDE wrapper

### Phase 11: RADAR5 (Stiff DDEs)

**Dependencies:** Phase 4, Phase 10
**Source:** `reference/delay/RADAR5-V2.2/radar5.f`

- [ ] `_wasm_radar5` — WASM — RADAR5 solver for stiff DDEs
- [ ] `_wasm_radar5_work_size` — WASM — Work array size
- [ ] `solveStiffDDE` — TS — High-level stiff DDE wrapper

### Phase 12: ROCK Methods (Stabilized Explicit)

**Dependencies:** Phase 1
**Source:** `reference/mol/rock/rock2.f`, `rock4.f`

- [ ] `_wasm_rock2` — WASM — ROCK2 solver
- [ ] `_wasm_rock4` — WASM — ROCK4 solver
- [ ] `_wasm_rock2_work_size` — WASM — Work array size
- [ ] `_wasm_rock4_work_size` — WASM — Work array size
- [ ] `solveRock` — TS — High-level ROCK wrapper
- [ ] Add `ROCK2`, `ROCK4` to method enum

### Phase 13: HEM5 (Constrained Mechanical Systems)

**Dependencies:** Phase 1
**Source:** `reference/mechanical/mechanic/hem5.f`

- [ ] `_wasm_hem5` — WASM — HEM5 solver
- [ ] `_wasm_hem5_work_size` — WASM — Work array size
- [ ] `MechanicalProblem` — TS — Interface for q', Mv', constraints
- [ ] `solveMechanical` — TS — High-level wrapper for DAE index-2

### Phase 14: GNI_IRK (Symplectic Methods)

**Dependencies:** Phase 1
**Source:** `reference/gnicodes/gni_irk2.c`

- [ ] `_wasm_gni_irk2` — WASM — Symplectic Gauss method
- [ ] `_wasm_gni_irk4` — WASM — 4-stage Gauss method
- [ ] `_wasm_gni_irk6` — WASM — 6-stage Gauss method
- [ ] `HamiltonianSystem` — TS — Interface for q'' = f(t, q)
- [ ] `solveHamiltonian` — TS — High-level symplectic solver

### Phase 15: VODE/ZVODE (Netlib ODEPACK)

**Dependencies:** Phase 1
**Source:** `reference/netlib/vode.f`, `zvode.f`

- [ ] `_wasm_dvode` — WASM — VODE solver (BDF/Adams)
- [ ] `_wasm_zvode` — WASM — ZVODE for complex ODEs
- [ ] `_wasm_vode_work_size` — WASM — Work array size
- [ ] `solveVode` — TS — High-level wrapper
- [ ] Add `VODE`, `BDF`, `Adams` to method enum

### Phase 16: RKF45 / RKSUITE (Classic Explicit)

**Dependencies:** Phase 1
**Source:** `reference/netlib/rkf45.f`, `rksuite.f`

- [ ] `_wasm_rkf45` — WASM — RKF45 solver
- [ ] `_wasm_rksuite_setup` — WASM — RKSUITE initialization
- [ ] `_wasm_rksuite_ut` — WASM — RKSUITE integration
- [ ] `solveRKF45` — TS — High-level wrapper
- [ ] `solveRKSuite` — TS — High-level wrapper

### Phase 17: RKC (Stabilized Chebyshev)

**Dependencies:** Phase 1
**Source:** `reference/netlib/rkc.f`

- [ ] `_wasm_rkc` — WASM — RKC solver
- [ ] `_wasm_rkcint` — WASM — RKC interpolation
- [ ] `_wasm_set_spcrad_callback` — WASM — Spectral radius callback
- [ ] `solveRKC` — TS — High-level wrapper
- [ ] Add `RKC` to method enum

### Phase 18: DDASSL/DDASRT (DAE Solvers)

**Dependencies:** Phase 1
**Source:** `reference/netlib/ddassl.f`, `ddasrt.f`

- [ ] `_wasm_ddassl` — WASM — DDASSL solver
- [ ] `_wasm_ddasrt` — WASM — DDASRT with root finding
- [ ] `_wasm_set_res_callback` — WASM — Residual callback G(t,y,y')
- [ ] `_wasm_set_root_callback` — WASM — Root function callback
- [ ] `DAEFunction` — TS — DAE residual type
- [ ] `solveDAE` — TS — High-level DAE wrapper

### Phase 19: COLNEW (Boundary Value Problems)

**Dependencies:** Phase 1
**Source:** `reference/netlib/colnew.f`

- [ ] `_wasm_colnew` — WASM — COLNEW solver
- [ ] `_wasm_set_fsub_callback` — WASM — ODE RHS callback
- [ ] `_wasm_set_gsub_callback` — WASM — Boundary condition callback
- [ ] `_wasm_set_dfsub_callback` — WASM — Jacobian of f
- [ ] `_wasm_set_dgsub_callback` — WASM — Jacobian of g
- [ ] `BVPProblem` — TS — BVP problem interface
- [ ] `solveBVP` — TS — High-level BVP wrapper

---

## Summary

| Category | Solvers | Implemented | Status |
|----------|---------|-------------|--------|
| **Hairer Non-Stiff** | 5 | 2 | 40% |
| **Hairer Stiff** | 5 | 1 | 20% |
| **Hairer Delay** | 1 | 0 | 0% |
| **Hairer MOL** | 2 | 0 | 0% |
| **Hairer Mechanical** | 2 | 0 | 0% |
| **Hairer Geometric** | 1 | 0 | 0% |
| **Netlib IVP** | 8 | 0 | 0% |
| **Netlib Stabilized** | 2 | 0 | 0% |
| **Netlib DAE** | 5 | 0 | 0% |
| **Netlib BVP** | 7 | 0 | 0% |
| **Netlib Special** | 3 | 0 | 0% |
| **Total Solvers** | **41** | **3** | **7%** |

### Core Functions Status

| Category | Complete | Partial | Not Started | Total |
|----------|----------|---------|-------------|-------|
| **Phase 1-5 (Core)** | 36 | 1 | 2 | 39 |
| **Phase 6-14 (Hairer)** | 0 | 0 | 51 | 51 |
| **Phase 15-19 (Netlib)** | 0 | 0 | ~40 | ~40 |
| **Total** | **36** | **1** | **~93** | **~130** |

**Core Implementation: 95% complete (36/39)**
**Overall Reference Coverage: ~28% complete (37/130)**

### Remaining Core Work

1. **Dense output interpolation in `solveExplicit`** — Currently t_eval only captures endpoints; should use `_wasm_contd5`/`_wasm_contd8` to interpolate at arbitrary points within steps
2. **`dense_output` in `SolveIVPResult`** — The `sol` field should be populated with a `DenseOutput` function when `dense_output=true`
3. **`t_eval` proper interpolation** — Use the dense output functions to evaluate at exact t_eval points

---

## Currently Implemented WASM Functions (19 total)

### Core Solvers

| Function | Description | Status |
|----------|-------------|--------|
| `_wasm_dopri5` | Dormand-Prince 5(4) explicit RK solver | ✅ |
| `_wasm_dop853` | Dormand-Prince 8(5,3) high-order explicit RK | ✅ |
| `_wasm_radau5` | Implicit Radau IIA order 5 for stiff systems | ✅ |

### Dense Output (Interpolation)

| Function | Description | Status |
|----------|-------------|--------|
| `_wasm_contd5` | Continuous output interpolation for DOPRI5 | ✅ |
| `_wasm_contd8` | Continuous output interpolation for DOP853 | ✅ |
| `_wasm_contr5` | Continuous output interpolation for RADAU5 | ✅ |

### Callback Setters

| Function | Description | Status |
|----------|-------------|--------|
| `_wasm_set_fcn_callback` | Set ODE function f(t,y) callback | ✅ |
| `_wasm_set_solout_callback` | Set solution output callback | ✅ |
| `_wasm_set_jac_callback` | Set Jacobian callback (for RADAU5) | ✅ |

### Work Array Sizing

| Function | Description | Status |
|----------|-------------|--------|
| `_wasm_dopri5_work_size` | Get work array size for DOPRI5 | ✅ |
| `_wasm_dopri5_iwork_size` | Get iwork array size for DOPRI5 | ✅ |
| `_wasm_dop853_work_size` | Get work array size for DOP853 | ✅ |
| `_wasm_dop853_iwork_size` | Get iwork array size for DOP853 | ✅ |
| `_wasm_radau5_work_size` | Get work array size for RADAU5 | ✅ |
| `_wasm_radau5_iwork_size` | Get iwork array size for RADAU5 | ✅ |

### Memory Management

| Function | Description | Status |
|----------|-------------|--------|
| `_wasm_malloc` | Allocate WASM memory | ✅ |
| `_wasm_free` | Free WASM memory | ✅ |
| `_malloc` | Emscripten malloc | ✅ |
| `_free` | Emscripten free | ✅ |

---

## Currently Implemented TypeScript Exports (15 total)

### High-Level API

| Export | Type | Description | Status |
|--------|------|-------------|--------|
| `solve_ivp` | function | SciPy-compatible IVP solver | ✅ |
| `solveExplicit` | function | Direct wrapper for DOPRI5/DOP853 | ✅ |
| `solveRadau` | function | Direct wrapper for RADAU5 | ✅ |

### Module Loading

| Export | Type | Description | Status |
|--------|------|-------------|--------|
| `loadODEModule` | function | Load and initialize WASM module | ✅ |
| `getODEModule` | function | Get loaded module (throws if not loaded) | ✅ |
| `isODELoaded` | function | Check if module is loaded | ✅ |
| `resetODEModule` | function | Reset module state | ✅ |
| `configureODE` | function | Configure before loading | ✅ |

### Types

| Export | Type | Description | Status |
|--------|------|-------------|--------|
| `ODEModule` | interface | WASM module interface | ✅ |
| `ODEFunction` | type | `(t: number, y: number[]) => number[]` | ✅ |
| `JacobianFunction` | type | `(t: number, y: number[]) => number[][]` | ✅ |
| `SolveIVPOptions` | interface | Options for solve_ivp | ✅ |
| `SolveIVPResult` | interface | Result from solve_ivp | ✅ |
| `DenseOutput` | interface | Interpolation function | ✅ (type only) |

### Enums & Constants

| Export | Type | Description | Status |
|--------|------|-------------|--------|
| `ExplicitMethod` | enum | `RK45`, `DOP853` | ✅ |
| `ImplicitMethod` | enum | `Radau` | ✅ |
| `STATUS_MESSAGES` | const | Status code → message mapping | ✅ |

---

## Reference File Locations

```
reference/
├── nonstiff/           # Hairer explicit methods for non-stiff problems
│   ├── nonstiff/
│   │   ├── dopri5.f    # ✅ Implemented
│   │   ├── dop853.f    # ✅ Implemented
│   │   ├── odex.f      # Variable-order extrapolation
│   │   ├── odex2.f     # ODEX variant
│   │   └── retard.f    # Delay differential equations
│   └── cprog/          # C versions of DOPRI5/DOP853
│
├── stiff/              # Hairer implicit methods for stiff problems
│   └── stiff/
│       ├── radau5.f    # ✅ Implemented
│       ├── radau.f     # Variable-order Radau
│       ├── radaup.f    # Preconditioned Radau
│       ├── rodas.f     # Rosenbrock method
│       ├── seulex.f    # Extrapolation with implicit Euler
│       ├── decsol.f    # LU decomposition routines
│       └── lapack.f    # LAPACK subset
│
├── delay/              # Delay differential equations
│   └── RADAR5-V2.2/
│       └── radar5.f    # Stiff DDEs
│
├── mol/                # Method of lines (parabolic PDEs)
│   └── rock/
│       ├── rock2.f     # 2nd order stabilized explicit
│       └── rock4.f     # 4th order stabilized explicit
│
├── mechanical/         # Constrained mechanical systems
│   └── mechanic/
│       ├── hem5.f      # Half-explicit method
│       └── phem56.f    # Partitioned half-explicit
│
├── gnicodes/           # Geometric numerical integration
│   └── gni_irk2.c      # Symplectic Gauss methods
│
└── netlib/             # Classic Netlib solvers
    ├── rkf45.f         # Runge-Kutta-Fehlberg 4(5)
    ├── dverk.f         # Verner 6(5)
    ├── ode.f           # Adams-Bashforth-Moulton
    ├── vode.f          # ODEPACK variable-coefficient
    ├── zvode.f         # Complex VODE
    ├── vodpk.f         # VODE with Krylov
    ├── rkc.f           # RK-Chebyshev
    ├── irkc.f90        # Implicit RKC
    ├── rksuite.f       # RK Suite (2,3)/(4,5)/(7,8)
    ├── ddassl.f        # DAE solver (BDF)
    ├── ddasrt.f        # DDASSL + root finding
    ├── mebdfi.f        # Modified extended BDF
    ├── colnew.f        # BVP collocation
    ├── colsys.f        # BVP collocation system
    ├── coldae.f        # BVP for DAEs
    ├── twpbvp.f        # Two-point BVP
    └── ...             # Many more
```

---

## Dependency Graph

```
Phase 1 (Foundation)
    │
    ├───────────────────────┬──────────────────────┬──────────────────────┐
    ↓                       ↓                      ↓                      ↓
Phase 2 (Explicit)    Phase 4 (Implicit)    Phase 6 (ODEX)    Phase 15-19 (Netlib)
    │                       │                                        │
    ↓                       ├──────────────────────┤                  │
Phase 3 (Dense)             ↓                      ↓                  │
    │                  Phase 7 (RODAS)       Phase 8 (SEULEX)         │
    │                       │                                         │
    │                  Phase 9 (RADAU VO)                             │
    │                       │                                         │
    └───────────┬───────────┘                                         │
                ↓                                                      │
         Phase 5 (Unified API) ←──────────────────────────────────────┘
                │
     ┌──────────┼──────────────────────────┐
     ↓          ↓                          ↓
Phase 10    Phase 12 (ROCK)          Phase 13 (HEM5)
(RETARD)        │                          │
     │          │                          │
     ↓          │                          │
Phase 11        │                          │
(RADAR5)        │                          │
                │                          │
                └──────────────────────────┘
                            │
                            ↓
                     Phase 14 (GNI_IRK)
```

---

## Usage Example

```typescript
import { solve_ivp } from 'odewasm';

// Solve dy/dt = -y (exponential decay)
const result = await solve_ivp(
  (t, y) => [-y[0]],  // ODE function
  [0, 5],              // time span
  [1],                 // initial condition
  { method: 'RK45', rtol: 1e-6 }
);

console.log(result.t);  // time points
console.log(result.y);  // solution values
```

---

## References

### Hairer Solvers
- [Hairer's ODE Solvers](http://www.unige.ch/~hairer/software.html)
- E. Hairer, S.P. Nørsett, G. Wanner: *Solving Ordinary Differential Equations I: Nonstiff Problems*
- E. Hairer, G. Wanner: *Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems*
- E. Hairer, C. Lubich, G. Wanner: *Geometric Numerical Integration*

### Netlib Solvers
- [Netlib ODE Repository](http://www.netlib.org/ode/)
- P.N. Brown, G.D. Byrne, A.C. Hindmarsh: *VODE: A Variable Coefficient ODE Solver*, SIAM J. Sci. Stat. Comput. 10 (1989)
- L.R. Petzold: *A Description of DASSL*, Sandia Report SAND82-8637
- U. Ascher, J. Christiansen, R.D. Russell: *A Collocation Solver for Mixed Order Systems of Boundary Value Problems*, Math. Comp. 33 (1979)

### SciPy Reference
- [SciPy solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
