/**
 * BFGS quasi-Newton optimization algorithm.
 *
 * Direct C translation of:
 *  - scipy.optimize._optimize._minimize_bfgs
 *  - scipy.optimize._dcsrch.DCSRCH (Moré-Thuente line search)
 *  - scipy.optimize._dcsrch.dcstep
 *
 * References:
 *   Nocedal, J., and Wright, S. J. "Numerical Optimization." Springer, 2006.
 *   Moré, J. J. and Thuente, D. J. "Line search algorithms with guaranteed
 *   sufficient decrease." ACM TOMS 20(3):286-307, 1994.
 */

#include "bfgs.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* ============================== */
/* Helper: vector operations      */
/* ============================== */

static double dot(const double* a, const double* b, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

static double vec_max_abs(const double* v, int n) {
    double m = 0.0;
    for (int i = 0; i < n; i++) {
        double a = fabs(v[i]);
        if (a > m) m = a;
    }
    return m;
}

static double vec_norm2(const double* v, int n) {
    return sqrt(dot(v, v, n));
}

static inline double dmax(double a, double b) { return a > b ? a : b; }
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dclip(double v, double lo, double hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}
static inline double dsign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

/* ============================== */
/* DCSTEP: safeguarded step       */
/* ============================== */

typedef struct {
    double stx, fx, dx;
    double sty, fy, dy;
    double stp;
    int brackt;
} DcstepState;

static void dcstep(
    double *stx, double *fx, double *dx,
    double *sty, double *fy, double *dy,
    double *stp, double fp, double dp,
    int *brackt, double stpmin, double stpmax
) {
    double sgnd = dsign(dp) * dsign(*dx);
    double stpf, stpc, stpq;
    double theta, s, gamma, p, q, r;

    if (fp > *fx) {
        /* Case 1: Higher function value. Minimum is bracketed. */
        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = dmax(fabs(theta), dmax(fabs(*dx), fabs(dp)));
        gamma = s * sqrt((theta / s) * (theta / s) - (*dx / s) * (dp / s));
        if (*stp < *stx) gamma = -gamma;
        p = (gamma - *dx) + theta;
        q = ((gamma - *dx) + gamma) + dp;
        r = p / q;
        stpc = *stx + r * (*stp - *stx);
        stpq = *stx + ((*dx / ((*fx - fp) / (*stp - *stx) + *dx)) / 2.0) * (*stp - *stx);
        if (fabs(stpc - *stx) <= fabs(stpq - *stx)) {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.0;
        }
        *brackt = 1;

    } else if (sgnd < 0.0) {
        /* Case 2: Lower value, opposite sign derivatives. Bracketed. */
        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = dmax(fabs(theta), dmax(fabs(*dx), fabs(dp)));
        gamma = s * sqrt((theta / s) * (theta / s) - (*dx / s) * (dp / s));
        if (*stp > *stx) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + *dx;
        r = p / q;
        stpc = *stp + r * (*stx - *stp);
        stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
        if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = 1;

    } else if (fabs(dp) < fabs(*dx)) {
        /* Case 3: Lower value, same sign, magnitude decreases. */
        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = dmax(fabs(theta), dmax(fabs(*dx), fabs(dp)));
        double disc = (theta / s) * (theta / s) - (*dx / s) * (dp / s);
        if (disc < 0.0) disc = 0.0;
        gamma = s * sqrt(disc);
        if (*stp > *stx) gamma = -gamma;
        p = (gamma - dp) + theta;
        q = (gamma + (*dx - dp)) + gamma;
        r = p / q;
        if (r < 0.0 && gamma != 0.0) {
            stpc = *stp + r * (*stx - *stp);
        } else if (*stp > *stx) {
            stpc = stpmax;
        } else {
            stpc = stpmin;
        }
        stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);

        if (*brackt) {
            if (fabs(stpc - *stp) < fabs(stpq - *stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            if (*stp > *stx) {
                stpf = dmin(*stp + 0.66 * (*sty - *stp), stpf);
            } else {
                stpf = dmax(*stp + 0.66 * (*sty - *stp), stpf);
            }
        } else {
            if (fabs(stpc - *stp) > fabs(stpq - *stp)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            stpf = dclip(stpf, stpmin, stpmax);
        }

    } else {
        /* Case 4: Lower value, same sign, magnitude does not decrease. */
        if (*brackt) {
            theta = 3.0 * (fp - *fy) / (*sty - *stp) + *dy + dp;
            s = dmax(fabs(theta), dmax(fabs(*dy), fabs(dp)));
            gamma = s * sqrt((theta / s) * (theta / s) - (*dy / s) * (dp / s));
            if (*stp > *sty) gamma = -gamma;
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + *dy;
            r = p / q;
            stpc = *stp + r * (*sty - *stp);
            stpf = stpc;
        } else if (*stp > *stx) {
            stpf = stpmax;
        } else {
            stpf = stpmin;
        }
    }

    /* Update the interval which contains a minimizer. */
    if (fp > *fx) {
        *sty = *stp;
        *fy = fp;
        *dy = dp;
    } else {
        if (sgnd < 0.0) {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = fp;
        *dx = dp;
    }

    *stp = stpf;
}


/* ============================== */
/* DCSRCH: Moré-Thuente search    */
/* ============================== */

/* Return codes for dcsrch */
#define DCSRCH_FG    0  /* Need function/gradient evaluation */
#define DCSRCH_CONV  1  /* Converged */
#define DCSRCH_WARN  2  /* Warning */
#define DCSRCH_ERROR 3  /* Error */

typedef struct {
    double ftol, gtol, xtol;
    double stpmin, stpmax;
    /* State */
    int stage;
    int brackt;
    double finit, ginit, gtest;
    double stx, fx, gx;
    double sty, fy, gy;
    double stmin, stmax;
    double width, width1;
} DcsrchState;

static void dcsrch_init(DcsrchState* s, double ftol, double gtol, double xtol,
                         double stpmin, double stpmax) {
    s->ftol = ftol;
    s->gtol = gtol;
    s->xtol = xtol;
    s->stpmin = stpmin;
    s->stpmax = stpmax;
    s->stage = 0; /* not yet started */
}

static int dcsrch_iterate(DcsrchState* s, double* stp, double f, double g) {
    const double p5 = 0.5;
    const double p66 = 0.66;
    const double xtrapl = 1.1;
    const double xtrapu = 4.0;

    if (s->stage == 0) {
        /* START */
        if (*stp < s->stpmin) return DCSRCH_ERROR;
        if (*stp > s->stpmax) return DCSRCH_ERROR;
        if (g >= 0.0) return DCSRCH_ERROR;
        if (s->ftol < 0.0) return DCSRCH_ERROR;
        if (s->gtol < 0.0) return DCSRCH_ERROR;
        if (s->xtol < 0.0) return DCSRCH_ERROR;
        if (s->stpmin < 0.0) return DCSRCH_ERROR;
        if (s->stpmax < s->stpmin) return DCSRCH_ERROR;

        s->brackt = 0;
        s->stage = 1;
        s->finit = f;
        s->ginit = g;
        s->gtest = s->ftol * s->ginit;
        s->width = s->stpmax - s->stpmin;
        s->width1 = s->width / p5;

        s->stx = 0.0;
        s->fx = s->finit;
        s->gx = s->ginit;
        s->sty = 0.0;
        s->fy = s->finit;
        s->gy = s->ginit;
        s->stmin = 0.0;
        s->stmax = *stp + xtrapu * (*stp);

        return DCSRCH_FG;
    }

    /* Non-start iteration */
    double ftest = s->finit + (*stp) * s->gtest;

    /* Stage transition */
    if (s->stage == 1 && f <= ftest && g >= 0.0) {
        s->stage = 2;
    }

    /* Test for warnings */
    if (s->brackt && (*stp <= s->stmin || *stp >= s->stmax))
        return DCSRCH_WARN;
    if (s->brackt && s->stmax - s->stmin <= s->xtol * s->stmax)
        return DCSRCH_WARN;
    if (*stp == s->stpmax && f <= ftest && g <= s->gtest)
        return DCSRCH_WARN;
    if (*stp == s->stpmin && (f > ftest || g >= s->gtest))
        return DCSRCH_WARN;

    /* Test for convergence */
    if (f <= ftest && fabs(g) <= s->gtol * (-s->ginit))
        return DCSRCH_CONV;

    /* Modified function for stage 1 */
    if (s->stage == 1 && f <= s->fx && f > ftest) {
        double fm = f - (*stp) * s->gtest;
        double fxm = s->fx - s->stx * s->gtest;
        double fym = s->fy - s->sty * s->gtest;
        double gm = g - s->gtest;
        double gxm = s->gx - s->gtest;
        double gym = s->gy - s->gtest;

        dcstep(&s->stx, &fxm, &gxm, &s->sty, &fym, &gym,
               stp, fm, gm, &s->brackt, s->stmin, s->stmax);

        s->fx = fxm + s->stx * s->gtest;
        s->fy = fym + s->sty * s->gtest;
        s->gx = gxm + s->gtest;
        s->gy = gym + s->gtest;
    } else {
        dcstep(&s->stx, &s->fx, &s->gx, &s->sty, &s->fy, &s->gy,
               stp, f, g, &s->brackt, s->stmin, s->stmax);
    }

    /* Bisection check */
    if (s->brackt) {
        if (fabs(s->sty - s->stx) >= p66 * s->width1) {
            *stp = s->stx + p5 * (s->sty - s->stx);
        }
        s->width1 = s->width;
        s->width = fabs(s->sty - s->stx);
    }

    /* Update step bounds */
    if (s->brackt) {
        s->stmin = dmin(s->stx, s->sty);
        s->stmax = dmax(s->stx, s->sty);
    } else {
        s->stmin = *stp + xtrapl * (*stp - s->stx);
        s->stmax = *stp + xtrapu * (*stp - s->stx);
    }

    /* Force step within global bounds */
    *stp = dclip(*stp, s->stpmin, s->stpmax);

    /* If further progress is not possible, use best point */
    if ((s->brackt && (*stp <= s->stmin || *stp >= s->stmax)) ||
        (s->brackt && s->stmax - s->stmin <= s->xtol * s->stmax)) {
        *stp = s->stx;
    }

    return DCSRCH_FG;
}


/* ============================== */
/* Line search with Wolfe conds   */
/* ============================== */

/**
 * Perform Wolfe line search.
 *
 * Given f, g, search direction pk:
 *   phi(alpha) = f(xk + alpha*pk)
 *   derphi(alpha) = g(xk + alpha*pk) · pk
 *
 * Returns step size alpha satisfying strong Wolfe conditions,
 * or -1.0 if line search failed.
 */
static double line_search_wolfe(
    int n,
    const double* xk,
    const double* pk,
    const double* gfk,
    double fk,
    double old_fval_prev, /* f_{k-1}, for initial step guess */
    double c1,
    double c2,
    double (*func_ptr)(const double* x, int n),
    void (*grad_fn)(const double* x, double* grad_out, int n),
    double* x_work,   /* scratch: length n */
    double* g_work,    /* scratch: length n */
    double* fval_out,  /* output: f(xk + alpha*pk) */
    double* gfkp1_out, /* output: gradient at new point, length n */
    int* nfev,
    int* njev
) {
    const int maxiter = 100;

    /* Directional derivative at alpha=0 */
    double derphi0 = dot(gfk, pk, n);
    if (derphi0 >= 0.0) return -1.0; /* Not a descent direction */

    /* Initial step guess (scipy heuristic) */
    double alpha1;
    if (old_fval_prev > fk) {
        alpha1 = dmin(1.0, 1.01 * 2.0 * (fk - old_fval_prev) / derphi0);
        if (alpha1 <= 0.0) alpha1 = 1.0;
    } else {
        alpha1 = 1.0;
    }

    /* Initialize DCSRCH */
    DcsrchState ds;
    dcsrch_init(&ds, c1, c2, 1e-14, 1e-100, 1e100);

    double stp = alpha1;
    double phi_val = fk;
    double dphi_val = derphi0;

    /* First call: START */
    int task = dcsrch_iterate(&ds, &stp, phi_val, dphi_val);
    if (task != DCSRCH_FG) return -1.0;

    for (int iter = 0; iter < maxiter; iter++) {
        /* Compute x_work = xk + stp * pk */
        for (int i = 0; i < n; i++) {
            x_work[i] = xk[i] + stp * pk[i];
        }

        /* Evaluate phi(stp) = f(x_work) */
        phi_val = func_ptr(x_work, n);
        (*nfev)++;

        /* Evaluate derphi(stp) = grad(x_work) · pk */
        grad_fn(x_work, g_work, n);
        (*njev)++;
        dphi_val = dot(g_work, pk, n);

        task = dcsrch_iterate(&ds, &stp, phi_val, dphi_val);

        if (!isfinite(stp)) return -1.0;

        if (task == DCSRCH_CONV) {
            /* Success! Compute final point values */
            for (int i = 0; i < n; i++) {
                x_work[i] = xk[i] + stp * pk[i];
            }
            *fval_out = phi_val;
            /* g_work already has the gradient at the last evaluated point */
            /* But the step may have been adjusted, so recompute if needed */
            /* Actually, DCSRCH_CONV means the last evaluated point satisfies Wolfe */
            memcpy(gfkp1_out, g_work, n * sizeof(double));
            return stp;
        }

        if (task == DCSRCH_WARN || task == DCSRCH_ERROR) {
            return -1.0;
        }

        /* task == DCSRCH_FG: continue with new stp */
    }

    return -1.0; /* maxiter reached */
}


/* ============================== */
/* Finite-difference gradient     */
/* ============================== */

static void finite_diff_gradient(
    const double* x, double* grad_out, int n,
    double eps,
    double (*func_ptr)(const double* x, int n),
    double* x_work,
    int* nfev
) {
    double f0 = func_ptr(x, n);
    (*nfev)++;

    memcpy(x_work, x, n * sizeof(double));
    for (int i = 0; i < n; i++) {
        double xi_orig = x_work[i];
        x_work[i] = xi_orig + eps;
        double fp = func_ptr(x_work, n);
        (*nfev)++;
        grad_out[i] = (fp - f0) / eps;
        x_work[i] = xi_orig;
    }
}


/* ============================== */
/* BFGS main algorithm            */
/* ============================== */

/* Context for gradient evaluation (handles both analytical and finite-diff) */
typedef struct {
    int has_jac;
    double eps;
    double (*func_ptr)(const double* x, int n);
    void (*grad_ptr)(const double* x, double* grad_out, int n);
    double* fd_work; /* scratch for finite diff */
    int* nfev_ptr;   /* pointer to nfev counter (for finite diff) */
    int n;
} GradContext;

static void eval_gradient(const double* x, double* grad_out, int n, GradContext* ctx) {
    if (ctx->has_jac) {
        ctx->grad_ptr(x, grad_out, n);
    } else {
        finite_diff_gradient(x, grad_out, n, ctx->eps, ctx->func_ptr,
                             ctx->fd_work, ctx->nfev_ptr);
    }
}

/* Wrapper for line search that uses GradContext */
static void grad_fn_wrapper_ctx(const double* x, double* grad_out, int n);
static GradContext* g_grad_ctx = NULL; /* Thread-unsafe but WASM is single-threaded */

static void grad_fn_wrapper_ctx(const double* x, double* grad_out, int n) {
    eval_gradient(x, grad_out, n, g_grad_ctx);
}


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
) {
    if (maxiter <= 0) maxiter = n * 200;
    if (eps <= 0.0) eps = 1.4901161193847656e-8; /* sqrt(machine epsilon) */

    int nfev = 0;
    int njev = 0;

    /* Allocate working memory */
    double* xk = (double*)malloc(n * sizeof(double));
    double* gfk = (double*)malloc(n * sizeof(double));
    double* pk = (double*)malloc(n * sizeof(double));
    double* Hk = (double*)malloc(n * n * sizeof(double));
    double* sk = (double*)malloc(n * sizeof(double));
    double* yk = (double*)malloc(n * sizeof(double));
    double* gfkp1 = (double*)malloc(n * sizeof(double));
    double* x_work = (double*)malloc(n * sizeof(double));
    double* g_work = (double*)malloc(n * sizeof(double));
    double* fd_work = (double*)malloc(n * sizeof(double));
    /* Temporaries for Hessian update */
    double* Hy = (double*)malloc(n * sizeof(double));
    double* A1Hk = (double*)malloc(n * n * sizeof(double));

    if (!xk || !gfk || !pk || !Hk || !sk || !yk || !gfkp1 ||
        !x_work || !g_work || !fd_work || !Hy || !A1Hk) {
        free(xk); free(gfk); free(pk); free(Hk); free(sk); free(yk);
        free(gfkp1); free(x_work); free(g_work); free(fd_work);
        free(Hy); free(A1Hk);
        result->status = -1;
        result->nfev = 0; result->njev = 0; result->nit = 0;
        result->fval = INFINITY;
        return -1;
    }

    /* Set up gradient context */
    GradContext grad_ctx;
    grad_ctx.has_jac = has_jac;
    grad_ctx.eps = eps;
    grad_ctx.func_ptr = func_ptr;
    grad_ctx.grad_ptr = grad_ptr;
    grad_ctx.fd_work = fd_work;
    grad_ctx.nfev_ptr = &nfev;
    grad_ctx.n = n;
    g_grad_ctx = &grad_ctx;

    /* Initialize x */
    memcpy(xk, x0, n * sizeof(double));

    /* Initial function evaluation */
    double old_fval = func_ptr(xk, n);
    nfev++;

    /* Initial gradient */
    eval_gradient(xk, gfk, n, &grad_ctx);
    njev++;

    /* Initialize Hk = Identity */
    memset(Hk, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) Hk[i * n + i] = 1.0;

    /* Initial step guess: old_old_fval = old_fval + ||gfk|| / 2 */
    double gnorm = vec_max_abs(gfk, n);
    double old_old_fval = old_fval + vec_norm2(gfk, n) / 2.0;

    int k = 0;
    int warnflag = 0;

    while (gnorm > gtol && k < maxiter) {
        /* Search direction: pk = -Hk @ gfk */
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                s += Hk[i * n + j] * gfk[j];
            }
            pk[i] = -s;
        }

        /* Line search */
        double new_fval;
        double alpha_k = line_search_wolfe(
            n, xk, pk, gfk, old_fval, old_old_fval,
            c1, c2, func_ptr, grad_fn_wrapper_ctx,
            x_work, g_work, &new_fval, gfkp1,
            &nfev, &njev
        );

        if (alpha_k < 0.0) {
            /* Line search failed */
            warnflag = 2;
            break;
        }

        /* sk = alpha_k * pk */
        for (int i = 0; i < n; i++) {
            sk[i] = alpha_k * pk[i];
        }

        /* Update x: xk = xk + sk */
        for (int i = 0; i < n; i++) {
            xk[i] += sk[i];
        }

        old_old_fval = old_fval;
        old_fval = new_fval;

        /* yk = gfkp1 - gfk */
        for (int i = 0; i < n; i++) {
            yk[i] = gfkp1[i] - gfk[i];
        }
        memcpy(gfk, gfkp1, n * sizeof(double));

        k++;
        gnorm = vec_max_abs(gfk, n);
        if (gnorm <= gtol) break;

        /* Check for non-finite values */
        if (!isfinite(old_fval)) {
            warnflag = 2;
            break;
        }

        /* BFGS Hessian inverse update:
         * rhok = 1 / (yk · sk)
         * A1 = I - rhok * sk * yk^T
         * A2 = I - rhok * yk * sk^T
         * Hk = A1 @ Hk @ A2 + rhok * sk * sk^T
         */
        double rhok_inv = dot(yk, sk, n);
        double rhok;
        if (rhok_inv == 0.0) {
            rhok = 1000.0;
        } else {
            rhok = 1.0 / rhok_inv;
        }

        /* Compute A1 @ Hk:
         * A1 = I - rhok * sk * yk^T
         * (A1 @ Hk)[i][j] = Hk[i][j] - rhok * sk[i] * (yk · Hk_col_j)
         *
         * Actually, let's compute it directly:
         * A1Hk[i][j] = sum_m A1[i][m] * Hk[m][j]
         *            = sum_m (delta[i][m] - rhok*sk[i]*yk[m]) * Hk[m][j]
         *            = Hk[i][j] - rhok*sk[i] * sum_m yk[m]*Hk[m][j]
         *
         * Let Hy[j] = sum_m yk[m] * Hk[m][j]  (yk^T @ Hk, row j of result)
         */
        /* Compute Hy = Hk^T @ yk  (since Hk is symmetric, this is Hk @ yk) */
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int m = 0; m < n; m++) {
                s += yk[m] * Hk[m * n + j];
            }
            Hy[j] = s;
        }

        /* A1Hk[i][j] = Hk[i][j] - rhok * sk[i] * Hy[j] */
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A1Hk[i * n + j] = Hk[i * n + j] - rhok * sk[i] * Hy[j];
            }
        }

        /* Hk_new = A1Hk @ A2 + rhok * sk * sk^T
         * A2 = I - rhok * yk * sk^T
         * (A1Hk @ A2)[i][j] = sum_m A1Hk[i][m] * A2[m][j]
         *                    = sum_m A1Hk[i][m] * (delta[m][j] - rhok*yk[m]*sk[j])
         *                    = A1Hk[i][j] - rhok * sk[j] * sum_m A1Hk[i][m] * yk[m]
         */
        /* Hk_new[i][j] = (A1Hk @ A2)[i][j] + rhok * sk[i] * sk[j]
         * where (A1Hk @ A2)[i][j] = A1Hk[i][j] - rhok * (A1Hk[i] · yk) * sk[j]
         */
        for (int i = 0; i < n; i++) {
            double ti = 0.0;
            for (int m = 0; m < n; m++) {
                ti += A1Hk[i * n + m] * yk[m];
            }
            for (int j = 0; j < n; j++) {
                Hk[i * n + j] = A1Hk[i * n + j] - rhok * ti * sk[j] + rhok * sk[i] * sk[j];
            }
        }
    }

    /* Write output */
    memcpy(x_out, xk, n * sizeof(double));
    memcpy(jac_out, gfk, n * sizeof(double));
    memcpy(hess_inv_out, Hk, n * n * sizeof(double));

    result->fval = old_fval;
    result->nfev = nfev;
    result->njev = njev;
    result->nit = k;

    if (warnflag == 2) {
        result->status = 2; /* line search or precision loss */
    } else if (k >= maxiter) {
        result->status = 1; /* maxiter */
    } else if (!isfinite(gnorm) || !isfinite(old_fval)) {
        result->status = 3; /* nan */
    } else {
        result->status = 0; /* success */
    }

    /* Cleanup */
    free(xk); free(gfk); free(pk); free(Hk); free(sk); free(yk);
    free(gfkp1); free(x_work); free(g_work); free(fd_work);
    free(Hy); free(A1Hk);
    g_grad_ctx = NULL;

    return 0;
}
