#ifndef SCIWASM_NELDER_MEAD_H
#define SCIWASM_NELDER_MEAD_H

/**
 * Result struct for Nelder-Mead optimization.
 * Layout must match the TypeScript binding's struct reading code.
 */
typedef struct {
    int nfev;      /* number of function evaluations */
    int nit;       /* number of iterations */
    int status;    /* 0=success, 1=maxfev reached, 2=maxiter reached */
    double fval;   /* final function value */
} NelderMeadResult;

/**
 * Nelder-Mead simplex optimization.
 *
 * @param n          Number of variables
 * @param x0         Initial guess (n doubles), not modified
 * @param x_out      Solution output (n doubles)
 * @param sim_out    Final simplex output ((n+1)*n doubles, row-major), can be NULL
 * @param fsim_out   Final simplex function values ((n+1) doubles), can be NULL
 * @param result     Result struct output
 * @param xatol      Absolute tolerance in x for convergence
 * @param fatol      Absolute tolerance in f for convergence
 * @param maxiter    Max iterations (0 = n*200)
 * @param maxfev     Max function evaluations (0 = n*200)
 * @param adaptive   1 for adaptive parameters, 0 for standard
 * @param has_bounds 1 if bounds provided
 * @param lower      Lower bounds (n doubles, ignored if !has_bounds)
 * @param upper      Upper bounds (n doubles, ignored if !has_bounds)
 * @param initial_simplex  Optional initial simplex ((n+1)*n doubles, row-major), NULL to auto-generate
 * @param func_ptr   Function pointer for objective: double f(const double* x, int n)
 * @return 0 on success
 */
int nelder_mead_minimize(
    int n,
    const double* x0,
    double* x_out,
    double* sim_out,
    double* fsim_out,
    NelderMeadResult* result,
    double xatol,
    double fatol,
    int maxiter,
    int maxfev,
    int adaptive,
    int has_bounds,
    const double* lower,
    const double* upper,
    const double* initial_simplex,
    double (*func_ptr)(const double* x, int n)
);

#endif
