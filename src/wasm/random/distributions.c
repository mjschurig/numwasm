/**
 * NumJS Random - Distribution Algorithms Implementation
 *
 * Implements various probability distributions using the BitGenerator interface.
 */

#include "distributions.h"
#include "ziggurat_constants.h"
#include <math.h>
#include <float.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Constants ============ */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============ Uniform Distributions ============ */

EXPORT double random_standard_uniform(bitgen_t *state) {
    return state->next_double(state->state);
}

EXPORT float random_standard_uniform_f(bitgen_t *state) {
    return uint32_to_float(state->next_uint32(state->state));
}

EXPORT void random_uniform_fill(bitgen_t *state, int32_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = state->next_double(state->state);
    }
}

EXPORT void random_uniform_fill_f(bitgen_t *state, int32_t cnt, float *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = uint32_to_float(state->next_uint32(state->state));
    }
}

/* ============ Normal Distribution (Ziggurat) ============ */

/**
 * Handle the unlikely cases in the Ziggurat algorithm for normal distribution.
 * This is called when the sample falls in the tail or needs rejection sampling.
 */
static double standard_normal_unlikely(bitgen_t *state, int idx, int sign, uint64_t rabs) {
    double x, xx, yy;

    if (idx == 0) {
        /* Tail: sample from exponential distribution */
        for (;;) {
            xx = -ZIGGURAT_NOR_INV_R * log(1.0 - state->next_double(state->state));
            yy = -log(1.0 - state->next_double(state->state));
            if (yy + yy > xx * xx) {
                return sign ? -(ZIGGURAT_NOR_R + xx) : (ZIGGURAT_NOR_R + xx);
            }
        }
    } else {
        /* Rejection sampling within the layer */
        x = rabs * wi_double[idx];
        if (sign) x = -x;

        if (((fi_double[idx - 1] - fi_double[idx]) * state->next_double(state->state) +
             fi_double[idx]) < exp(-0.5 * x * x)) {
            return x;
        }
        /* Retry - recursively call main function */
        return random_standard_normal(state);
    }
}

EXPORT double random_standard_normal(bitgen_t *state) {
    uint64_t r;
    int sign;
    uint64_t rabs;
    int idx;
    double x;

    r = state->next_uint64(state->state);
    idx = r & 0xFF;           /* Use bottom 8 bits for layer index */
    r >>= 8;
    sign = r & 0x1;           /* Use next bit for sign */
    rabs = (r >> 1) & 0x000FFFFFFFFFFFFFULL;  /* Use remaining 53 bits */
    x = rabs * wi_double[idx];

    if (sign) x = -x;

    /* Fast path: ~98.9% of samples */
    if (rabs < ki_double[idx]) {
        return x;
    }

    /* Slow path: tail or rejection sampling */
    return standard_normal_unlikely(state, idx, sign, rabs);
}

/**
 * Single precision normal using float Ziggurat tables.
 */
static float standard_normal_unlikely_f(bitgen_t *state, int idx, int sign, uint32_t rabs) {
    float x, xx, yy;

    if (idx == 0) {
        /* Tail */
        for (;;) {
            xx = -ZIGGURAT_NOR_INV_R_F * logf(1.0f - uint32_to_float(state->next_uint32(state->state)));
            yy = -logf(1.0f - uint32_to_float(state->next_uint32(state->state)));
            if (yy + yy > xx * xx) {
                return sign ? -(ZIGGURAT_NOR_R_F + xx) : (ZIGGURAT_NOR_R_F + xx);
            }
        }
    } else {
        x = rabs * wi_float[idx];
        if (sign) x = -x;

        if (((fi_float[idx - 1] - fi_float[idx]) * uint32_to_float(state->next_uint32(state->state)) +
             fi_float[idx]) < expf(-0.5f * x * x)) {
            return x;
        }
        return random_standard_normal_f(state);
    }
}

EXPORT float random_standard_normal_f(bitgen_t *state) {
    uint32_t r;
    int sign;
    uint32_t rabs;
    int idx;
    float x;

    r = state->next_uint32(state->state);
    idx = r & 0xFF;
    r >>= 8;
    sign = r & 0x1;
    rabs = (r >> 1) & 0x7FFFFF;  /* 23 bits for mantissa */
    x = rabs * wi_float[idx];

    if (sign) x = -x;

    if (rabs < ki_float[idx]) {
        return x;
    }

    return standard_normal_unlikely_f(state, idx, sign, rabs);
}

EXPORT void random_standard_normal_fill(bitgen_t *state, int32_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_normal(state);
    }
}

EXPORT void random_standard_normal_fill_f(bitgen_t *state, int32_t cnt, float *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_normal_f(state);
    }
}

EXPORT double random_normal(bitgen_t *state, double loc, double scale) {
    return loc + scale * random_standard_normal(state);
}

/* ============ Exponential Distribution (Ziggurat) ============ */

/**
 * Handle unlikely cases for exponential Ziggurat.
 */
static double standard_exponential_unlikely(bitgen_t *state, int idx, double x) {
    if (idx == 0) {
        /* Tail: add exponential to the tail boundary */
        return ZIGGURAT_EXP_R - log(1.0 - state->next_double(state->state));
    } else if ((fe_double[idx - 1] - fe_double[idx]) * state->next_double(state->state) +
               fe_double[idx] < exp(-x)) {
        return x;
    } else {
        /* Retry */
        return random_standard_exponential(state);
    }
}

EXPORT double random_standard_exponential(bitgen_t *state) {
    uint64_t ri;
    int idx;
    double x;

    ri = state->next_uint64(state->state);
    ri >>= 3;                    /* Discard bottom 3 bits */
    idx = ri & 0xFF;             /* Use bottom 8 bits for layer */
    ri >>= 8;                    /* Use remaining bits for position */
    x = ri * we_double[idx];

    /* Fast path: ~98.9% of samples */
    if (ri < ke_double[idx]) {
        return x;
    }

    return standard_exponential_unlikely(state, idx, x);
}

static float standard_exponential_unlikely_f(bitgen_t *state, int idx, float x) {
    if (idx == 0) {
        return ZIGGURAT_EXP_R_F - logf(1.0f - uint32_to_float(state->next_uint32(state->state)));
    } else if ((fe_float[idx - 1] - fe_float[idx]) * uint32_to_float(state->next_uint32(state->state)) +
               fe_float[idx] < expf(-x)) {
        return x;
    } else {
        return random_standard_exponential_f(state);
    }
}

EXPORT float random_standard_exponential_f(bitgen_t *state) {
    uint32_t ri;
    int idx;
    float x;

    ri = state->next_uint32(state->state);
    idx = ri & 0xFF;
    ri >>= 8;
    x = ri * we_float[idx];

    if (ri < ke_float[idx]) {
        return x;
    }

    return standard_exponential_unlikely_f(state, idx, x);
}

EXPORT void random_standard_exponential_fill(bitgen_t *state, int32_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_exponential(state);
    }
}

EXPORT void random_standard_exponential_fill_f(bitgen_t *state, int32_t cnt, float *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_exponential_f(state);
    }
}

EXPORT void random_standard_exponential_inv_fill(bitgen_t *state, int32_t cnt, double *out) {
    /* Inverse CDF method: -log(1 - U) */
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = -log(1.0 - state->next_double(state->state));
    }
}

EXPORT double random_exponential(bitgen_t *state, double scale) {
    return scale * random_standard_exponential(state);
}

/* ============ Gamma Distribution (Marsaglia-Tsang) ============ */

EXPORT double random_standard_gamma(bitgen_t *state, double shape) {
    double b, c, U, V, X, Y;

    if (shape == 1.0) {
        /* Gamma(1) = Exp(1) */
        return random_standard_exponential(state);
    } else if (shape == 0.0) {
        return 0.0;
    } else if (shape < 1.0) {
        /* Ahrens-Dieter method for shape < 1 */
        for (;;) {
            U = state->next_double(state->state);
            V = random_standard_exponential(state);
            if (U <= 1.0 - shape) {
                X = pow(U, 1.0 / shape);
                if (X <= V) {
                    return X;
                }
            } else {
                Y = -log((1.0 - U) / shape);
                X = pow(1.0 - shape + shape * Y, 1.0 / shape);
                if (X <= V + Y) {
                    return X;
                }
            }
        }
    } else {
        /* Marsaglia-Tsang method for shape > 1 */
        b = shape - 1.0 / 3.0;
        c = 1.0 / sqrt(9.0 * b);
        for (;;) {
            do {
                X = random_standard_normal(state);
                V = 1.0 + c * X;
            } while (V <= 0.0);

            V = V * V * V;
            U = state->next_double(state->state);

            /* Squeeze step: accept with high probability */
            if (U < 1.0 - 0.0331 * (X * X) * (X * X)) {
                return b * V;
            }

            /* Slower check using log */
            if (log(U) < 0.5 * X * X + b * (1.0 - V + log(V))) {
                return b * V;
            }
        }
    }
}

EXPORT float random_standard_gamma_f(bitgen_t *state, float shape) {
    float b, c, U, V, X, Y;

    if (shape == 1.0f) {
        return random_standard_exponential_f(state);
    } else if (shape == 0.0f) {
        return 0.0f;
    } else if (shape < 1.0f) {
        for (;;) {
            U = uint32_to_float(state->next_uint32(state->state));
            V = random_standard_exponential_f(state);
            if (U <= 1.0f - shape) {
                X = powf(U, 1.0f / shape);
                if (X <= V) {
                    return X;
                }
            } else {
                Y = -logf((1.0f - U) / shape);
                X = powf(1.0f - shape + shape * Y, 1.0f / shape);
                if (X <= V + Y) {
                    return X;
                }
            }
        }
    } else {
        b = shape - 1.0f / 3.0f;
        c = 1.0f / sqrtf(9.0f * b);
        for (;;) {
            do {
                X = random_standard_normal_f(state);
                V = 1.0f + c * X;
            } while (V <= 0.0f);

            V = V * V * V;
            U = uint32_to_float(state->next_uint32(state->state));

            if (U < 1.0f - 0.0331f * (X * X) * (X * X)) {
                return b * V;
            }

            if (logf(U) < 0.5f * X * X + b * (1.0f - V + logf(V))) {
                return b * V;
            }
        }
    }
}

EXPORT void random_standard_gamma_fill(bitgen_t *state, int32_t cnt, double shape, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_gamma(state, shape);
    }
}

EXPORT double random_gamma(bitgen_t *state, double shape, double scale) {
    return scale * random_standard_gamma(state, shape);
}

/* ============ Beta Distribution ============ */

EXPORT double random_beta(bitgen_t *state, double a, double b) {
    double Ga, Gb;

    if (a <= 1.0 && b <= 1.0) {
        /* Johnk's algorithm for small parameters */
        double U, V, X, Y;
        for (;;) {
            U = state->next_double(state->state);
            V = state->next_double(state->state);
            X = pow(U, 1.0 / a);
            Y = pow(V, 1.0 / b);

            if (X + Y <= 1.0) {
                if (X + Y > 0) {
                    return X / (X + Y);
                } else {
                    /* Handle underflow for very small a, b */
                    double logX = log(U) / a;
                    double logY = log(V) / b;
                    double logM = logX > logY ? logX : logY;
                    logX -= logM;
                    logY -= logM;
                    return exp(logX - log(exp(logX) + exp(logY)));
                }
            }
        }
    } else {
        /* Use ratio of independent gamma variates */
        Ga = random_standard_gamma(state, a);
        Gb = random_standard_gamma(state, b);
        return Ga / (Ga + Gb);
    }
}

EXPORT void random_beta_fill(bitgen_t *state, int32_t cnt, double a, double b, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_beta(state, a, b);
    }
}

/* ============ Chi-Square Distribution ============ */

EXPORT double random_chisquare(bitgen_t *state, double df) {
    return 2.0 * random_standard_gamma(state, df / 2.0);
}

EXPORT double random_noncentral_chisquare(bitgen_t *state, double df, double nonc) {
    if (nonc == 0.0) {
        return random_chisquare(state, df);
    }
    if (df > 1.0) {
        double Chi2 = random_chisquare(state, df - 1.0);
        double n = random_standard_normal(state) + sqrt(nonc);
        return Chi2 + n * n;
    } else {
        int64_t i = random_poisson(state, nonc / 2.0);
        return random_chisquare(state, df + 2.0 * i);
    }
}

/* ============ F Distribution ============ */

EXPORT double random_f(bitgen_t *state, double dfnum, double dfden) {
    return (random_chisquare(state, dfnum) * dfden) /
           (random_chisquare(state, dfden) * dfnum);
}

EXPORT double random_noncentral_f(bitgen_t *state, double dfnum, double dfden, double nonc) {
    double t = random_noncentral_chisquare(state, dfnum, nonc) * dfden;
    return t / (random_chisquare(state, dfden) * dfnum);
}

/* ============ Student's t Distribution ============ */

EXPORT double random_standard_t(bitgen_t *state, double df) {
    double N = random_standard_normal(state);
    double G = random_standard_gamma(state, df / 2.0);
    return sqrt(df / 2.0) * N / sqrt(G);
}

/* ============ Other Continuous Distributions ============ */

EXPORT double random_standard_cauchy(bitgen_t *state) {
    /* Ratio of two standard normals gives Cauchy */
    return random_standard_normal(state) / random_standard_normal(state);
}

EXPORT double random_pareto(bitgen_t *state, double a) {
    return exp(random_standard_exponential(state) / a) - 1.0;
}

EXPORT double random_weibull(bitgen_t *state, double a) {
    if (a == 0.0) {
        return 0.0;
    }
    return pow(random_standard_exponential(state), 1.0 / a);
}

EXPORT double random_power(bitgen_t *state, double a) {
    return pow(1.0 - exp(-random_standard_exponential(state)), 1.0 / a);
}

EXPORT double random_laplace(bitgen_t *state, double loc, double scale) {
    double U = state->next_double(state->state);
    if (U < 0.5) {
        return loc + scale * log(2.0 * U);
    } else {
        return loc - scale * log(2.0 * (1.0 - U));
    }
}

EXPORT double random_gumbel(bitgen_t *state, double loc, double scale) {
    double U = 1.0 - state->next_double(state->state);
    return loc - scale * log(-log(U));
}

EXPORT double random_logistic(bitgen_t *state, double loc, double scale) {
    double U = state->next_double(state->state);
    return loc + scale * log(U / (1.0 - U));
}

EXPORT double random_lognormal(bitgen_t *state, double mean, double sigma) {
    return exp(random_normal(state, mean, sigma));
}

EXPORT double random_rayleigh(bitgen_t *state, double scale) {
    return scale * sqrt(2.0 * random_standard_exponential(state));
}

EXPORT double random_wald(bitgen_t *state, double mean, double scale) {
    double mu_2l = mean / (2.0 * scale);
    double Y = random_standard_normal(state);
    Y = mean * Y * Y;
    double X = mean + mu_2l * (Y - sqrt(4.0 * scale * Y + Y * Y));
    double U = state->next_double(state->state);
    if (U <= mean / (mean + X)) {
        return X;
    } else {
        return mean * mean / X;
    }
}

EXPORT double random_triangular(bitgen_t *state, double left, double mode, double right) {
    double base = right - left;
    double leftbase = mode - left;
    double ratio = leftbase / base;
    double U = state->next_double(state->state);
    if (U <= ratio) {
        return left + sqrt(U * base * leftbase);
    } else {
        return right - sqrt((1.0 - U) * base * (right - mode));
    }
}

EXPORT double random_vonmises(bitgen_t *state, double mu, double kappa) {
    double s;
    double U, V, W, Y, Z;
    double result, mod;
    int neg;

    if (kappa < 1e-8) {
        /* Uniform for very small kappa */
        return M_PI * (2.0 * state->next_double(state->state) - 1.0);
    }

    /* Best-Fisher algorithm */
    s = 0.5 / kappa;
    double r = s + sqrt(1.0 + s * s);

    for (;;) {
        U = state->next_double(state->state);
        Z = cos(M_PI * U);
        W = (1.0 + r * Z) / (r + Z);
        Y = kappa * (r - W);
        V = state->next_double(state->state);

        if (Y * (2.0 - Y) - V >= 0.0) {
            break;
        }
        if (log(Y / V) + 1.0 - Y >= 0.0) {
            break;
        }
    }

    U = state->next_double(state->state);
    result = acos(W);
    if (U < 0.5) {
        result = -result;
    }
    result += mu;

    /* Normalize to [-pi, pi) */
    neg = (result < 0);
    mod = fabs(result);
    mod = fmod(mod + M_PI, 2.0 * M_PI) - M_PI;
    if (neg) {
        mod *= -1;
    }

    return mod;
}

/* ============ Discrete Distributions ============ */

/**
 * Binomial distribution using BTPE algorithm.
 * BTPE = BinomialTriangleParetoExponential
 * Efficient for large n*min(p, 1-p).
 */
EXPORT int64_t random_binomial_btpe(bitgen_t *state, binomial_t *binomial, int64_t n, double p) {
    double r, q, fm, p1, xm, xl, xr, c, laml, lamr, p2, p3, p4;
    double a, U, V, s, F, rho, t, A, nrq, x1, x2, f1, f2, z, z2, w, w2;
    int64_t m, y, k;

    if (!(binomial->initialized) || (binomial->nsave != n) || (binomial->psave != p)) {
        /* Initialize cached values */
        binomial->nsave = n;
        binomial->psave = p;
        binomial->initialized = 1;

        r = p < 0.5 ? p : 1.0 - p;
        q = 1.0 - r;
        fm = n * r + r;
        m = (int64_t)fm;
        p1 = floor(2.195 * sqrt(n * r * q) - 4.6 * q) + 0.5;
        xm = m + 0.5;
        xl = xm - p1;
        xr = xm + p1;
        c = 0.134 + 20.5 / (15.3 + m);
        a = (fm - xl) / (fm - xl * r);
        laml = a * (1.0 + a / 2.0);
        a = (xr - fm) / (xr * q);
        lamr = a * (1.0 + a / 2.0);
        p2 = p1 * (1.0 + 2.0 * c);
        p3 = p2 + c / laml;
        p4 = p3 + c / lamr;

        binomial->r = r;
        binomial->q = q;
        binomial->fm = fm;
        binomial->m = m;
        binomial->p1 = p1;
        binomial->xm = xm;
        binomial->xl = xl;
        binomial->xr = xr;
        binomial->c = c;
        binomial->laml = laml;
        binomial->lamr = lamr;
        binomial->p2 = p2;
        binomial->p3 = p3;
        binomial->p4 = p4;
    } else {
        r = binomial->r;
        q = binomial->q;
        fm = binomial->fm;
        m = binomial->m;
        p1 = binomial->p1;
        xm = binomial->xm;
        xl = binomial->xl;
        xr = binomial->xr;
        c = binomial->c;
        laml = binomial->laml;
        lamr = binomial->lamr;
        p2 = binomial->p2;
        p3 = binomial->p3;
        p4 = binomial->p4;
    }

step10:
    nrq = n * r * q;
    U = state->next_double(state->state) * p4;
    V = state->next_double(state->state);

    if (U <= p1) {
        /* Triangular region */
        y = (int64_t)(xm - p1 * V + U);
        goto step50;
    }

    if (U <= p2) {
        /* Parallelogram region */
        double x = xl + (U - p1) / c;
        V = V * c + 1.0 - fabs(m - x + 0.5) / p1;
        if (V > 1.0) goto step10;
        y = (int64_t)x;
        goto step50;
    }

    if (U <= p3) {
        /* Left exponential tail */
        y = (int64_t)(xl + log(V) / laml);
        if (y < 0) goto step10;
        V *= (U - p2) * laml;
    } else {
        /* Right exponential tail */
        y = (int64_t)(xr - log(V) / lamr);
        if (y > n) goto step10;
        V *= (U - p3) * lamr;
    }

    /* Step 40: acceptance/rejection comparison */
    k = llabs(y - m);
    if ((k > 20) && (k < (nrq / 2.0 - 1))) {
        /* Squeeze step */
        s = r / q;
        a = s * (n + 1);
        F = 1.0;
        if (m < y) {
            for (int64_t i = m + 1; i <= y; i++) {
                F *= a / i - s;
            }
        } else if (m > y) {
            for (int64_t i = y + 1; i <= m; i++) {
                F /= a / i - s;
            }
        }
        if (V > F) goto step10;
        goto step50;
    }

    /* Step 50: use Stirling approximation */
step50:
    if (p > 0.5) {
        y = n - y;
    }
    return y;
}

/**
 * Binomial using inversion method.
 * Efficient for small n*min(p, 1-p).
 */
EXPORT int64_t random_binomial_inversion(bitgen_t *state, int64_t n, double p) {
    double q, qn, r, U;
    int64_t x;

    if (p > 0.5) {
        p = 1.0 - p;
        x = n;
    } else {
        x = 0;
    }

    q = 1.0 - p;
    qn = exp(n * log(q));
    r = p / q;
    U = state->next_double(state->state);

    while (U > qn) {
        U -= qn;
        qn *= r * (n - x) / (x + 1);
        x++;
    }

    return x;
}

EXPORT int64_t random_binomial(bitgen_t *state, double p, int64_t n) {
    double q;

    if (p <= 0.0 || n == 0) {
        return 0;
    }
    if (p >= 1.0) {
        return n;
    }

    q = 1.0 - p;
    if (p <= 0.5) {
        if (n * p < 30.0) {
            return random_binomial_inversion(state, n, p);
        } else {
            binomial_t binomial = {0};
            return random_binomial_btpe(state, &binomial, n, p);
        }
    } else {
        if (n * q < 30.0) {
            return n - random_binomial_inversion(state, n, q);
        } else {
            binomial_t binomial = {0};
            return random_binomial_btpe(state, &binomial, n, p);
        }
    }
}

/* 32-bit wrapper for JavaScript compatibility */
EXPORT int32_t random_binomial32(bitgen_t *state, double p, int32_t n) {
    return (int32_t)random_binomial(state, p, (int64_t)n);
}

EXPORT int64_t random_negative_binomial(bitgen_t *state, double n, double p) {
    double Y = random_standard_gamma(state, n);
    return random_poisson(state, Y * (1.0 - p) / p);
}

/* 32-bit wrapper for JavaScript compatibility */
EXPORT int32_t random_negative_binomial32(bitgen_t *state, double n, double p) {
    return (int32_t)random_negative_binomial(state, n, p);
}

EXPORT int64_t random_poisson(bitgen_t *state, double lam) {
    if (lam >= 10.0) {
        /* PTRS algorithm for large lambda */
        double slam = sqrt(lam);
        double loglam = log(lam);
        double b = 0.931 + 2.53 * slam;
        double a = -0.059 + 0.02483 * b;
        double invalpha = 1.1239 + 1.1328 / (b - 3.4);
        double vr = 0.9277 - 3.6224 / (b - 2.0);

        for (;;) {
            double U, V, us, k;
            U = state->next_double(state->state) - 0.5;
            V = state->next_double(state->state);
            us = 0.5 - fabs(U);
            k = floor((2.0 * a / us + b) * U + lam + 0.43);

            if (us >= 0.07 && V <= vr) {
                return (int64_t)k;
            }
            if (k < 0 || (us < 0.013 && V > us)) {
                continue;
            }
            if (log(V) + log(invalpha) - log(a / (us * us) + b) <=
                -lam + k * loglam - lgamma(k + 1)) {
                return (int64_t)k;
            }
        }
    } else if (lam == 0.0) {
        return 0;
    } else {
        /* Direct method for small lambda */
        double enlam = exp(-lam);
        int64_t k = 0;
        double prod = 1.0;

        for (;;) {
            prod *= state->next_double(state->state);
            if (prod > enlam) {
                k++;
            } else {
                return k;
            }
        }
    }
}

EXPORT int64_t random_geometric(bitgen_t *state, double p) {
    if (p >= 1.0) {
        return 1;
    }
    if (p <= 0.0) {
        return 0;
    }
    return (int64_t)ceil(log(1.0 - state->next_double(state->state)) / log(1.0 - p));
}

/* 32-bit wrappers for JavaScript compatibility */
EXPORT int32_t random_poisson32(bitgen_t *state, double lam) {
    return (int32_t)random_poisson(state, lam);
}

EXPORT int32_t random_geometric32(bitgen_t *state, double p) {
    return (int32_t)random_geometric(state, p);
}

EXPORT int64_t random_hypergeometric(bitgen_t *state, int64_t ngood, int64_t nbad, int64_t nsample) {
    int64_t d1, K, Z;
    double d2;

    if (nsample > 10) {
        /* Hypergeometric random variate - ratio method */
        int64_t total = ngood + nbad;
        d1 = nsample;
        d2 = (double)ngood + (double)nbad - nsample;
        Z = 0;

        while (d1 > 0) {
            double Y = (double)ngood / (double)(ngood + nbad);
            K = random_binomial(state, Y, d1);
            Z += K;
            ngood -= K;
            d1 -= K;

            Y = (double)nbad / (double)(ngood + nbad);
            K = random_binomial(state, Y, d1);
            nbad -= K;
            d1 -= K;
        }

        return Z;
    } else {
        /* Direct method for small nsample */
        int64_t good = ngood;
        int64_t bad = nbad;
        int64_t selected = 0;
        int64_t total = ngood + nbad;

        for (int64_t i = 0; i < nsample; i++) {
            double p = (double)good / (double)total;
            if (state->next_double(state->state) < p) {
                selected++;
                good--;
            } else {
                bad--;
            }
            total--;
        }

        return selected;
    }
}

EXPORT int64_t random_logseries(bitgen_t *state, double p) {
    double q, r, U, V, result;

    r = log(1.0 - p);

    for (;;) {
        V = state->next_double(state->state);
        if (V >= p) {
            return 1;
        }
        U = state->next_double(state->state);
        q = 1.0 - exp(r * U);
        if (V <= q * q) {
            result = floor(1.0 + log(V) / log(q));
            if (result < 1.0) {
                continue;
            }
            return (int64_t)result;
        }
        if (V >= q) {
            return 1;
        }
        return 2;
    }
}

/* 32-bit wrappers for JavaScript compatibility */
EXPORT int32_t random_hypergeometric32(bitgen_t *state, int32_t ngood, int32_t nbad, int32_t nsample) {
    return (int32_t)random_hypergeometric(state, (int64_t)ngood, (int64_t)nbad, (int64_t)nsample);
}

EXPORT int32_t random_logseries32(bitgen_t *state, double p) {
    return (int32_t)random_logseries(state, p);
}

EXPORT int64_t random_zipf(bitgen_t *state, double a) {
    double am1, b;

    am1 = a - 1.0;
    b = pow(2.0, am1);

    for (;;) {
        double T, U, V, X;

        U = 1.0 - state->next_double(state->state);
        V = state->next_double(state->state);
        X = floor(pow(U, -1.0 / am1));

        if (X < 1.0) {
            continue;
        }

        T = pow(1.0 + 1.0 / X, am1);

        if (V * X * (T - 1.0) / (b - 1.0) <= T / b) {
            return (int64_t)X;
        }
    }
}

/* 32-bit wrapper for JavaScript compatibility */
EXPORT int32_t random_zipf32(bitgen_t *state, double a) {
    return (int32_t)random_zipf(state, a);
}

/* ============ Bounded Integers (Lemire's Algorithm) ============ */

EXPORT uint64_t random_bounded_uint64(bitgen_t *state, uint64_t off, uint64_t rng, uint64_t mask) {
    uint64_t val;

    if (rng == 0) {
        return off;
    } else if ((rng & (rng + 1)) == 0) {
        /* Power of 2: use simple masking */
        return off + (state->next_uint64(state->state) & rng);
    } else {
        /* Lemire's algorithm for unbiased bounded integers */
        /* threshold = 2^64 % (rng + 1) = (-rng - 1) % (rng + 1) */
        uint64_t threshold = (-rng - 1) % (rng + 1);
        for (;;) {
            val = state->next_uint64(state->state);
            if (val >= threshold) {
                return off + (val % (rng + 1));
            }
        }
    }
}

EXPORT uint32_t random_bounded_uint32(bitgen_t *state, uint32_t off, uint32_t rng, uint32_t mask) {
    uint32_t val;

    if (rng == 0) {
        return off;
    } else if ((rng & (rng + 1)) == 0) {
        return off + (state->next_uint32(state->state) & rng);
    } else {
        uint32_t threshold = (-rng - 1) % (rng + 1);
        for (;;) {
            val = state->next_uint32(state->state);
            if (val >= threshold) {
                return off + (val % (rng + 1));
            }
        }
    }
}

EXPORT int64_t random_integers(bitgen_t *state, int64_t low, int64_t high) {
    if (low == high) {
        return low;
    }
    uint64_t rng = (uint64_t)(high - low);
    return low + (int64_t)random_bounded_uint64(state, 0, rng, 0);
}

/* 32-bit wrapper for JavaScript compatibility */
EXPORT int32_t random_integers32(bitgen_t *state, int32_t low, int32_t high) {
    if (low == high) {
        return low;
    }
    uint32_t rng = (uint32_t)(high - low);
    return low + (int32_t)random_bounded_uint32(state, 0, rng, 0);
}

EXPORT void random_integers_fill(bitgen_t *state, int32_t cnt, int64_t low, int64_t high, int64_t *out) {
    uint64_t rng = (uint64_t)(high - low);
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = low + (int64_t)random_bounded_uint64(state, 0, rng, 0);
    }
}

/* 32-bit wrapper for JavaScript compatibility */
EXPORT void random_integers32_fill(bitgen_t *state, int32_t cnt, int32_t low, int32_t high, int32_t *out) {
    uint32_t rng = (uint32_t)(high - low);
    for (int32_t i = 0; i < cnt; i++) {
        out[i] = low + (int32_t)random_bounded_uint32(state, 0, rng, 0);
    }
}
