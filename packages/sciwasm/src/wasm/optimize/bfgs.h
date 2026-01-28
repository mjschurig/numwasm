/**
 * BFGS quasi-Newton optimization algorithm.
 *
 * C translation of scipy.optimize._optimize._minimize_bfgs
 * with Moré-Thuente line search (MINPACK-2 DCSRCH/DCSTEP).
 */

#ifndef BFGS_H
#define BFGS_H

/**
 * Result struct for BFGS optimization.
 * Layout must match TypeScript expectations.
 */
typedef struct {
    int nfev;    /* Number of function evaluations */
    int njev;    /* Number of gradient evaluations */
    int nit;     /* Number of iterations */
    int status;  /* 0=success, 1=maxiter, 2=line_search_failed, 3=nan */
    double fval; /* Final function value */
} BFGSResult;

/**
 * Minimize a scalar function using the BFGS algorithm.
 *
 * @param n          Dimension of the problem
 * @param x0         Initial point (length n)
 * @param x_out      Output: solution vector (length n)
 * @param hess_inv_out Output: inverse Hessian approximation (n*n), row-major
 * @param jac_out    Output: final gradient (length n)
 * @param result     Output: result struct
 * @param gtol       Gradient tolerance for convergence
 * @param maxiter    Maximum number of iterations (0 = n*200)
 * @param c1         Armijo condition parameter
 * @param c2         Curvature condition parameter
 * @param has_jac    1 if grad_ptr computes analytical gradient, 0 for finite diff
 * @param eps        Step size for finite-difference gradient (if has_jac=0)
 * @param func_ptr   Objective function: f(x, n) -> double
 * @param grad_ptr   Gradient function: grad(x, grad_out, n) -> void
 *                   If has_jac=0, this is ignored and finite differences are used.
 * @return 0 on success
 */
int bfgs_minimize(
    int n,
    const double* x0,
    double* x_out,
    double* hess_inv_out,
    double* jac_out,
    BFGSResult* result,
    double gtol,
    int maxiter,
    double c1,
    double c2,
    int has_jac,
    double eps,
    double (*func_ptr)(const double* x, int n),
    void (*grad_ptr)(const double* x, double* grad_out, int n)
);

#endif /* BFGS_H */
