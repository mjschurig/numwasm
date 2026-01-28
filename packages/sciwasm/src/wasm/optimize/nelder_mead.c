/**
 * Nelder-Mead simplex optimization algorithm.
 *
 * Direct C translation of scipy.optimize._optimize._minimize_neldermead
 * (scipy/optimize/_optimize.py, lines 700-968).
 *
 * References:
 *   Gao, F. and Han, L. "Implementing the Nelder-Mead simplex algorithm
 *   with adaptive parameters." Computational Optimization and Applications,
 *   51:1, pp. 259-277, 2012.
 */

#include "nelder_mead.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Clip value to [lo, hi] */
static double clip(double val, double lo, double hi) {
    if (val < lo) return lo;
    if (val > hi) return hi;
    return val;
}

/* Simple insertion sort for argsort of small arrays */
static void argsort(const double* vals, int* indices, int n) {
    for (int i = 0; i < n; i++) indices[i] = i;
    for (int i = 1; i < n; i++) {
        int key_idx = indices[i];
        double key_val = vals[key_idx];
        int j = i - 1;
        while (j >= 0 && vals[indices[j]] > key_val) {
            indices[j + 1] = indices[j];
            j--;
        }
        indices[j + 1] = key_idx;
    }
}

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
) {
    const int np1 = n + 1;

    /* Adaptive or standard parameters */
    double rho, chi, psi, sigma;
    if (adaptive) {
        double dim = (double)n;
        rho = 1.0;
        chi = 1.0 + 2.0 / dim;
        psi = 0.75 - 1.0 / (2.0 * dim);
        sigma = 1.0 - 1.0 / dim;
    } else {
        rho = 1.0;
        chi = 2.0;
        psi = 0.5;
        sigma = 0.5;
    }

    const double nonzdelt = 0.05;
    const double zdelt = 0.00025;

    /* Default maxiter/maxfev */
    if (maxiter <= 0 && maxfev <= 0) {
        maxiter = n * 200;
        maxfev = n * 200;
    } else if (maxiter <= 0) {
        maxiter = (int)2e9; /* effectively infinity */
    } else if (maxfev <= 0) {
        maxfev = (int)2e9;
    }

    /* Allocate simplex: (n+1) vertices, each of dimension n */
    double* sim = (double*)malloc(np1 * n * sizeof(double));
    double* fsim = (double*)malloc(np1 * sizeof(double));
    int* ind = (int*)malloc(np1 * sizeof(int));
    double* xbar = (double*)malloc(n * sizeof(double));
    double* xr = (double*)malloc(n * sizeof(double));
    double* xe = (double*)malloc(n * sizeof(double));
    double* xc = (double*)malloc(n * sizeof(double));
    double* x0_clipped = (double*)malloc(n * sizeof(double));

    if (!sim || !fsim || !ind || !xbar || !xr || !xe || !xc || !x0_clipped) {
        free(sim); free(fsim); free(ind); free(xbar);
        free(xr); free(xe); free(xc); free(x0_clipped);
        result->status = -1;
        result->nfev = 0;
        result->nit = 0;
        result->fval = INFINITY;
        return -1;
    }

    /* Clip x0 to bounds if needed */
    memcpy(x0_clipped, x0, n * sizeof(double));
    if (has_bounds) {
        for (int i = 0; i < n; i++) {
            x0_clipped[i] = clip(x0_clipped[i], lower[i], upper[i]);
        }
    }

    /* Initialize simplex */
    if (initial_simplex) {
        memcpy(sim, initial_simplex, np1 * n * sizeof(double));
    } else {
        /* Row 0 = x0 */
        memcpy(sim, x0_clipped, n * sizeof(double));
        /* Rows 1..n: perturb each dimension */
        for (int k = 0; k < n; k++) {
            memcpy(sim + (k + 1) * n, x0_clipped, n * sizeof(double));
            double val = x0_clipped[k];
            if (val != 0.0) {
                sim[(k + 1) * n + k] = (1.0 + nonzdelt) * val;
            } else {
                sim[(k + 1) * n + k] = zdelt;
            }
        }
    }

    /* Clip simplex to bounds and handle reflection for degenerate cases */
    if (has_bounds) {
        for (int i = 0; i < np1; i++) {
            for (int j = 0; j < n; j++) {
                double v = sim[i * n + j];
                if (v > upper[j]) {
                    /* Reflect into interior */
                    v = 2.0 * upper[j] - v;
                }
                /* Clip to bounds */
                sim[i * n + j] = clip(v, lower[j], upper[j]);
            }
        }
    }

    /* Evaluate initial simplex */
    int nfev = 0;
    for (int k = 0; k < np1; k++) {
        fsim[k] = INFINITY;
    }
    for (int k = 0; k < np1; k++) {
        if (nfev >= maxfev) break;
        fsim[k] = func_ptr(sim + k * n, n);
        nfev++;
    }

    /* Sort simplex by function value */
    argsort(fsim, ind, np1);
    {
        double* sim_tmp = (double*)malloc(np1 * n * sizeof(double));
        double* fsim_tmp = (double*)malloc(np1 * sizeof(double));
        if (sim_tmp && fsim_tmp) {
            for (int i = 0; i < np1; i++) {
                memcpy(sim_tmp + i * n, sim + ind[i] * n, n * sizeof(double));
                fsim_tmp[i] = fsim[ind[i]];
            }
            memcpy(sim, sim_tmp, np1 * n * sizeof(double));
            memcpy(fsim, fsim_tmp, np1 * sizeof(double));
        }
        free(sim_tmp);
        free(fsim_tmp);
    }

    int iterations = 1;

    while (nfev < maxfev && iterations < maxiter) {
        /* Check convergence */
        double max_x_diff = 0.0;
        for (int i = 1; i < np1; i++) {
            for (int j = 0; j < n; j++) {
                double d = fabs(sim[i * n + j] - sim[j]);
                if (d > max_x_diff) max_x_diff = d;
            }
        }
        double max_f_diff = 0.0;
        for (int i = 1; i < np1; i++) {
            double d = fabs(fsim[0] - fsim[i]);
            if (d > max_f_diff) max_f_diff = d;
        }
        if (max_x_diff <= xatol && max_f_diff <= fatol) {
            break;
        }

        /* Compute centroid of all vertices except the worst */
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) { /* n vertices (excluding last = worst) */
                s += sim[i * n + j];
            }
            xbar[j] = s / n;
        }

        /* Reflection */
        for (int j = 0; j < n; j++) {
            xr[j] = (1.0 + rho) * xbar[j] - rho * sim[(np1 - 1) * n + j];
        }
        if (has_bounds) {
            for (int j = 0; j < n; j++) {
                xr[j] = clip(xr[j], lower[j], upper[j]);
            }
        }
        if (nfev >= maxfev) break;
        double fxr = func_ptr(xr, n);
        nfev++;

        int doshrink = 0;

        if (fxr < fsim[0]) {
            /* Expansion */
            for (int j = 0; j < n; j++) {
                xe[j] = (1.0 + rho * chi) * xbar[j] - rho * chi * sim[(np1 - 1) * n + j];
            }
            if (has_bounds) {
                for (int j = 0; j < n; j++) {
                    xe[j] = clip(xe[j], lower[j], upper[j]);
                }
            }
            if (nfev >= maxfev) {
                /* Accept reflection before breaking */
                memcpy(sim + (np1 - 1) * n, xr, n * sizeof(double));
                fsim[np1 - 1] = fxr;
                goto sort_and_continue;
            }
            double fxe = func_ptr(xe, n);
            nfev++;

            if (fxe < fxr) {
                memcpy(sim + (np1 - 1) * n, xe, n * sizeof(double));
                fsim[np1 - 1] = fxe;
            } else {
                memcpy(sim + (np1 - 1) * n, xr, n * sizeof(double));
                fsim[np1 - 1] = fxr;
            }
        } else {
            /* fxr >= fsim[0] */
            if (fxr < fsim[np1 - 2]) {
                /* Accept reflection */
                memcpy(sim + (np1 - 1) * n, xr, n * sizeof(double));
                fsim[np1 - 1] = fxr;
            } else {
                /* Contraction */
                if (fxr < fsim[np1 - 1]) {
                    /* Outside contraction */
                    for (int j = 0; j < n; j++) {
                        xc[j] = (1.0 + psi * rho) * xbar[j] - psi * rho * sim[(np1 - 1) * n + j];
                    }
                    if (has_bounds) {
                        for (int j = 0; j < n; j++) {
                            xc[j] = clip(xc[j], lower[j], upper[j]);
                        }
                    }
                    if (nfev >= maxfev) goto sort_and_continue;
                    double fxc = func_ptr(xc, n);
                    nfev++;

                    if (fxc <= fxr) {
                        memcpy(sim + (np1 - 1) * n, xc, n * sizeof(double));
                        fsim[np1 - 1] = fxc;
                    } else {
                        doshrink = 1;
                    }
                } else {
                    /* Inside contraction */
                    for (int j = 0; j < n; j++) {
                        xc[j] = (1.0 - psi) * xbar[j] + psi * sim[(np1 - 1) * n + j];
                    }
                    if (has_bounds) {
                        for (int j = 0; j < n; j++) {
                            xc[j] = clip(xc[j], lower[j], upper[j]);
                        }
                    }
                    if (nfev >= maxfev) goto sort_and_continue;
                    double fxcc = func_ptr(xc, n);
                    nfev++;

                    if (fxcc < fsim[np1 - 1]) {
                        memcpy(sim + (np1 - 1) * n, xc, n * sizeof(double));
                        fsim[np1 - 1] = fxcc;
                    } else {
                        doshrink = 1;
                    }
                }

                if (doshrink) {
                    /* Shrink all vertices toward the best */
                    for (int i = 1; i < np1; i++) {
                        for (int j = 0; j < n; j++) {
                            sim[i * n + j] = sim[j] + sigma * (sim[i * n + j] - sim[j]);
                        }
                        if (has_bounds) {
                            for (int j = 0; j < n; j++) {
                                sim[i * n + j] = clip(sim[i * n + j], lower[j], upper[j]);
                            }
                        }
                        if (nfev >= maxfev) break;
                        fsim[i] = func_ptr(sim + i * n, n);
                        nfev++;
                    }
                }
            }
        }

sort_and_continue:
        iterations++;

        /* Sort simplex by function value */
        argsort(fsim, ind, np1);
        {
            double* sim_tmp = (double*)malloc(np1 * n * sizeof(double));
            double* fsim_tmp = (double*)malloc(np1 * sizeof(double));
            if (sim_tmp && fsim_tmp) {
                for (int i = 0; i < np1; i++) {
                    memcpy(sim_tmp + i * n, sim + ind[i] * n, n * sizeof(double));
                    fsim_tmp[i] = fsim[ind[i]];
                }
                memcpy(sim, sim_tmp, np1 * n * sizeof(double));
                memcpy(fsim, fsim_tmp, np1 * sizeof(double));
            }
            free(sim_tmp);
            free(fsim_tmp);
        }
    }

    /* Write output */
    memcpy(x_out, sim, n * sizeof(double)); /* best vertex = sim[0] */
    result->fval = fsim[0];
    result->nfev = nfev;
    result->nit = iterations;

    if (nfev >= maxfev) {
        result->status = 1;
    } else if (iterations >= maxiter) {
        result->status = 2;
    } else {
        result->status = 0;
    }

    /* Copy final simplex if requested */
    if (sim_out) {
        memcpy(sim_out, sim, np1 * n * sizeof(double));
    }
    if (fsim_out) {
        memcpy(fsim_out, fsim, np1 * sizeof(double));
    }

    free(sim);
    free(fsim);
    free(ind);
    free(xbar);
    free(xr);
    free(xe);
    free(xc);
    free(x0_clipped);

    return 0;
}
