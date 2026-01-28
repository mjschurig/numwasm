# sciwasm Implementation Todos

Inventory of all modules/functions needed to replace stubs and fully implement the sciwasm package, based on the scipy reference implementation in `/scipy`.

Legend: ✅ = implemented, 🔲 = stubbed (exists but throws NotImplementedError), ⬜ = not yet created

---

## Workflow

For each function/module, follow these steps in order:

### Step 1: Identify the underlying C/Fortran code in `/packages/sciwasm/reference/scipy`
- Look at the reference scipy implementation for the function
- Trace through the Python code to find the underlying C, C++, or Fortran source (e.g., in `scipy/special/cephes/`, `scipy/linalg/src/`, LAPACK/BLAS wrappers, etc.)
- These are the computational kernels that do the real work

### Step 2: Compile C/C++/Fortran to WebAssembly
- Copy the identified C/C++/Fortran source files into the sciwasm build
- Compile them to `.wasm` using Emscripten (or similar toolchain)
- Expose the necessary entry points for calling from JS/TS

### Step 3: Write TypeScript interfaces only
- Write thin TypeScript wrappers that call into the compiled WASM module
- Do **not** reimplement numerical algorithms in TypeScript — the TypeScript layer is purely an interface/binding layer
- Handle type conversions, array marshalling, and ergonomic API surface in TS

### Step 4: Identify and copy reference tests from `/scipy`
- Find the corresponding test files in the scipy repo (e.g., `scipy/special/tests/test_basic.py`)
- Copy the test cases — these are the ground truth for correctness

### Step 5: Port tests to TypeScript
- Translate the copied Python tests into the sciwasm TypeScript test suite
- Only adjust what is necessary for the language difference (syntax, imports, assertion style)
- Do **not** weaken or remove test cases — preserve the same coverage and numerical tolerances
- Wire the tests into the existing test runner

---

## 1. `special` — Special Mathematical Functions

### Current Stubs
- ✅ `gamma(z)` — Gamma function
- ✅ `gammaln(x)` — Log-gamma function
- 🔲 `beta(a, b)` — Beta function
- 🔲 `erf(z)` — Error function
- 🔲 `erfc(x)` — Complementary error function
- 🔲 `j0(x)` — Bessel J₀
- 🔲 `j1(x)` — Bessel J₁
- ✅ `factorial(n, exact?)` — Factorial
- 🔲 `comb(N, k, exact?)` — Binomial coefficient
- 🔲 `perm(N, k, exact?)` — Permutations

### Priority Additions
- ⬜ `betaln(a, b)` — Log of beta function
- ⬜ `gammainc(a, x)` — Regularized lower incomplete gamma
- ⬜ `gammaincc(a, x)` — Regularized upper incomplete gamma
- ⬜ `digamma(x)` / `psi(x)` — Digamma function
- ⬜ `loggamma(z)` — Principal branch of log-gamma
- ⬜ `erfinv(y)` — Inverse error function
- ⬜ `erfcinv(y)` — Inverse complementary error function
- ⬜ `erfcx(x)` — Scaled complementary error function
- ⬜ `y0(x)` — Bessel Y₀
- ⬜ `y1(x)` — Bessel Y₁
- ⬜ `jv(v, z)` — Bessel Jν (arbitrary order)
- ⬜ `yv(v, z)` — Bessel Yν (arbitrary order)
- ⬜ `iv(v, z)` — Modified Bessel Iν
- ⬜ `kv(v, z)` — Modified Bessel Kν
- ⬜ `expit(x)` — Sigmoid / logistic function
- ⬜ `logit(x)` — Logit function
- ⬜ `softmax(x)` — Softmax
- ⬜ `log_softmax(x)` — Log-softmax
- ⬜ `logsumexp(a)` — Log of sum of exponentials
- ⬜ `xlogy(x, y)` — x * log(y), 0 if x == 0
- ⬜ `xlog1py(x, y)` — x * log1p(y), 0 if x == 0
- ⬜ `exprel(x)` — Relative error exponential
- ⬜ `zeta(x)` — Riemann zeta function
- ⬜ `lambertw(z)` — Lambert W function
- ✅ `factorial2(n)` — Double factorial
- ✅ `factorialk(n, k)` — Multifactorial
- ✅ `rgamma(x)` — Reciprocal gamma function
- ⬜ `stirling2(N, K)` — Stirling numbers of the second kind
- ⬜ `ndtr(x)` — Standard normal CDF
- ⬜ `ndtri(y)` — Inverse standard normal CDF
- ⬜ `boxcox(x, lmbda)` — Box-Cox transformation
- ⬜ `boxcox1p(x, lmbda)` — Box-Cox for 1+x
- ⬜ `inv_boxcox(y, lmbda)` — Inverse Box-Cox
- ⬜ `inv_boxcox1p(y, lmbda)` — Inverse Box-Cox for 1+x
- ⬜ `eval_legendre(n, x)` — Evaluate Legendre polynomial
- ⬜ `eval_chebyt(n, x)` — Evaluate Chebyshev T polynomial
- ⬜ `eval_hermite(n, x)` — Evaluate Hermite polynomial
- ⬜ `eval_laguerre(n, x)` — Evaluate Laguerre polynomial
- ⬜ `hyp2f1(a, b, c, z)` — Gauss hypergeometric 2F1
- ⬜ `hyp1f1(a, b, z)` — Confluent hypergeometric 1F1
- ⬜ `ellipk(m)` — Complete elliptic integral K
- ⬜ `ellipe(m)` — Complete elliptic integral E
- ⬜ `sinc(x)` — Sinc function
- ⬜ `cbrt(x)` — Cube root
- ⬜ `entr(x)` — Elementwise entropy
- ⬜ `rel_entr(x, y)` — Relative entropy
- ⬜ `kl_div(x, y)` — Kullback-Leibler divergence
- ⬜ `softplus(x)` — Softplus function

---

## 2. `stats` — Statistical Functions

### Current Stubs
- ✅ `describe(a)` — Descriptive statistics
- 🔲 `norm(loc?, scale?)` — Normal distribution
- 🔲 `t(df, loc?, scale?)` — Student's t distribution
- 🔲 `f(dfn, dfd, loc?, scale?)` — F distribution
- 🔲 `chi2(df, loc?, scale?)` — Chi-squared distribution
- 🔲 `pearsonr(x, y)` — Pearson correlation
- 🔲 `spearmanr(a, b?)` — Spearman rank correlation
- 🔲 `ttest_ind(a, b, options?)` — Independent samples t-test
- 🔲 `ttest_1samp(a, popmean)` — One-sample t-test
- 🔲 `kstest(rvs, cdf)` — Kolmogorov-Smirnov test

### Priority Additions — Distributions
- ⬜ `uniform(loc?, scale?)` — Uniform distribution
- ⬜ `expon(loc?, scale?)` — Exponential distribution
- ⬜ `gamma(a, loc?, scale?)` — Gamma distribution
- ⬜ `beta(a, b, loc?, scale?)` — Beta distribution
- ⬜ `lognorm(s, loc?, scale?)` — Log-normal distribution
- ⬜ `weibull_min(c, loc?, scale?)` — Weibull minimum distribution
- ⬜ `cauchy(loc?, scale?)` — Cauchy distribution
- ⬜ `laplace(loc?, scale?)` — Laplace distribution
- ⬜ `pareto(b, loc?, scale?)` — Pareto distribution
- ⬜ `poisson(mu, loc?)` — Poisson distribution (discrete)
- ⬜ `binom(n, p, loc?)` — Binomial distribution (discrete)
- ⬜ `bernoulli(p, loc?)` — Bernoulli distribution (discrete)

Each distribution needs: `pdf`/`pmf`, `cdf`, `ppf`, `rvs`, `mean`, `std`, `var`, `entropy`, `fit`

### Priority Additions — Statistical Tests
- ⬜ `ttest_rel(a, b)` — Paired t-test
- ⬜ `mannwhitneyu(x, y)` — Mann-Whitney U test
- ⬜ `wilcoxon(x, y?)` — Wilcoxon signed-rank test
- ⬜ `kruskal(*args)` — Kruskal-Wallis H test
- ⬜ `f_oneway(*args)` — One-way ANOVA
- ⬜ `chi2_contingency(observed)` — Chi-square contingency test
- ⬜ `fisher_exact(table)` — Fisher's exact test
- ⬜ `shapiro(x)` — Shapiro-Wilk normality test
- ⬜ `normaltest(a)` — D'Agostino-Pearson normality test
- ⬜ `anderson(x)` — Anderson-Darling test
- ⬜ `levene(*args)` — Levene's test for equal variances
- ⬜ `bartlett(*args)` — Bartlett's test for equal variances
- ⬜ `chisquare(f_obs, f_exp?)` — Chi-square goodness of fit
- ⬜ `power_divergence(f_obs, f_exp?)` — Power divergence statistic

### Priority Additions — Correlation & Regression
- ⬜ `kendalltau(x, y)` — Kendall's tau
- ⬜ `pointbiserialr(x, y)` — Point-biserial correlation
- ⬜ `linregress(x, y)` — Simple linear regression

### Priority Additions — Descriptive
- ⬜ `mode(a)` — Modal value
- ✅ `moment(a, moment?)` — Central moment
- ✅ `skew(a)` — Skewness
- ✅ `kurtosis(a)` — Kurtosis
- ⬜ `sem(a)` — Standard error of the mean
- ⬜ `zscore(a)` — Z-score standardization
- ⬜ `iqr(x)` — Interquartile range
- ⬜ `trim_mean(a, proportiontocut)` — Trimmed mean
- ⬜ `rankdata(a)` — Rank data
- ⬜ `percentileofscore(a, score)` — Percentile of score
- ⬜ `scoreatpercentile(a, per)` — Score at percentile

### Priority Additions — Other
- ⬜ `entropy(pk, qk?)` — Shannon entropy
- ⬜ `differential_entropy(values)` — Differential entropy
- ⬜ `gaussian_kde(dataset)` — Kernel density estimation
- ⬜ `rv_continuous` — Base class for continuous distributions
- ⬜ `rv_discrete` — Base class for discrete distributions
- ⬜ `bootstrap(data, statistic)` — Bootstrap confidence intervals
- ⬜ `permutation_test(data, statistic)` — Permutation test

---

## 3. `optimize` — Optimization and Root Finding

### Current Stubs
- ✅ `minimize(fun, x0, options?)` — Minimize scalar function of one or more variables (Nelder-Mead, BFGS, L-BFGS-B)
- 🔲 `least_squares(fun, x0)` — Nonlinear least-squares
- 🔲 `root_scalar(f, options?)` — Find root of scalar function
- 🔲 `linprog(c, options?)` — Linear programming
- 🔲 `curve_fit(f, xdata, ydata, p0?)` — Nonlinear curve fitting

### Priority Additions
- ⬜ `minimize_scalar(fun, options?)` — Minimize scalar function of one variable
- ⬜ `root(fun, x0, options?)` — Find root of vector function
- ⬜ `brentq(f, a, b)` — Brent's method root finding
- ⬜ `brenth(f, a, b)` — Brent's method (hyperbolic extrapolation)
- ⬜ `bisect(f, a, b)` — Bisection root finding
- ⬜ `newton(func, x0)` — Newton-Raphson root finding
- ⬜ `ridder(f, a, b)` — Ridder's method
- ⬜ `toms748(f, a, b)` — TOMS 748 root finding
- ⬜ `fixed_point(func, x0)` — Fixed-point iteration
- ⬜ `differential_evolution(func, bounds)` — Global optimization
- ⬜ `basinhopping(func, x0)` — Global optimization (basin hopping)
- ⬜ `dual_annealing(func, bounds)` — Global optimization (simulated annealing)
- ⬜ `linear_sum_assignment(cost_matrix)` — Hungarian algorithm
- ⬜ `milp(c, constraints?)` — Mixed-integer linear programming
- ⬜ `nnls(A, b)` — Non-negative least squares
- ⬜ `lsq_linear(A, b)` — Bounded linear least squares
- ⬜ `bracket(func)` — Bracket a minimum
- ⬜ `approx_fprime(xk, f)` — Finite-difference gradient approximation
- ⬜ `check_grad(func, grad, x0)` — Check gradient correctness
- ✅ `OptimizeResult` — Result class
- ✅ `Bounds` — Variable bounds
- ⬜ `LinearConstraint` — Linear constraint
- ⬜ `NonlinearConstraint` — Nonlinear constraint

---

## 4. `integrate` — Integration and ODEs

### Current Stubs
- ✅ `quad(func, a, b, options?)` — Adaptive quadrature
- 🔲 `dblquad(func, a, b, gfun, hfun)` — Double integration
- 🔲 `tplquad(func, a, b, gfun, hfun, qfun, rfun)` — Triple integration
- 🔲 `trapezoid(y, x?, dx?)` — Trapezoidal rule
- 🔲 `simpson(y, x?, dx?)` — Simpson's rule
- 🔲 `odeint(func, y0, t, options?)` — ODE solver (legacy)

### Priority Additions
- ⬜ `quad_vec(func, a, b)` — Vector-valued quadrature
- ⬜ `nquad(func, ranges)` — N-dimensional quadrature
- ⬜ `fixed_quad(func, a, b, n?)` — Fixed-order Gaussian quadrature
- ⬜ `cumulative_trapezoid(y, x?)` — Cumulative trapezoidal
- ⬜ `cumulative_simpson(y, x?)` — Cumulative Simpson's
- ⬜ `romb(y, dx?)` — Romberg integration
- ⬜ `solve_ivp(fun, t_span, y0, method?)` — Modern ODE solver (RK45, RK23, DOP853, Radau, BDF, LSODA)
- ⬜ `solve_bvp(fun, bc, x, y)` — Boundary value problem solver
- ⬜ `newton_cotes(rn, equal?)` — Newton-Cotes integration weights

---

## 5. `interpolate` — Interpolation

### Current Stubs
- 🔲 `interp1d(x, y, options?)` — 1-D interpolation (function)
- 🔲 `CubicSpline` — Cubic spline interpolation (class)
- 🔲 `PchipInterpolator` — PCHIP monotonic cubic (class)
- 🔲 `griddata(points, values, xi, method?)` — Unstructured N-D interpolation

### Priority Additions
- ⬜ `Akima1DInterpolator(x, y)` — Akima 1-D interpolation
- ⬜ `BarycentricInterpolator(xi, yi?)` — Barycentric interpolation
- ⬜ `KroghInterpolator(xi, yi)` — Krogh interpolation
- ⬜ `CubicHermiteSpline(x, y, dydx)` — Cubic Hermite spline
- ⬜ `BSpline(t, c, k)` — B-spline basis
- ⬜ `PPoly(c, x)` — Piecewise polynomial
- ⬜ `BPoly(c, x)` — Bernstein polynomial basis
- ⬜ `make_interp_spline(x, y, k?)` — Build interpolating B-spline
- ⬜ `make_lsq_spline(x, y, t, k?)` — Build least-squares B-spline
- ⬜ `make_smoothing_spline(x, y)` — Build smoothing spline
- ⬜ `RegularGridInterpolator(points, values)` — N-D regular grid interpolation
- ⬜ `LinearNDInterpolator(points, values)` — Piecewise linear N-D
- ⬜ `NearestNDInterpolator(points, values)` — Nearest-neighbor N-D
- ⬜ `CloughTocher2DInterpolator(points, values)` — Clough-Tocher 2-D
- ⬜ `RBFInterpolator(y, d)` — Radial basis function interpolation
- ⬜ `NdPPoly(c, x)` — N-D piecewise polynomial
- ⬜ `NdBSpline(t, c, k)` — N-D B-spline

---

## 6. `signal` — Signal Processing

### Current Stubs
- 🔲 `convolve(in1, in2, mode?)` — 1-D convolution
- 🔲 `fftconvolve(in1, in2, mode?)` — FFT-based convolution
- 🔲 `butter(N, Wn, options?)` — Butterworth filter design
- 🔲 `sosfilt(sos, x)` — Second-order sections filtering
- 🔲 `firwin(numtaps, cutoff, options?)` — FIR filter design (window method)
- 🔲 `welch(x, options?)` — Power spectral density (Welch)
- 🔲 `spectrogram(x, options?)` — Spectrogram

### Priority Additions — Convolution & Correlation
- ⬜ `correlate(in1, in2, mode?)` — Cross-correlation
- ⬜ `convolve2d(in1, in2, mode?)` — 2-D convolution
- ⬜ `correlate2d(in1, in2, mode?)` — 2-D cross-correlation
- ⬜ `oaconvolve(in1, in2, mode?)` — Overlap-add convolution
- ⬜ `correlation_lags(in1_len, in2_len, mode?)` — Lag indices for correlation

### Priority Additions — Filtering
- ⬜ `lfilter(b, a, x)` — IIR/FIR filter
- ⬜ `filtfilt(b, a, x)` — Zero-phase filtering
- ⬜ `sosfiltfilt(sos, x)` — Zero-phase SOS filtering
- ⬜ `medfilt(volume, kernel_size?)` — Median filter
- ⬜ `wiener(im, mysize?)` — Wiener filter
- ⬜ `savgol_filter(x, window_length, polyorder)` — Savitzky-Golay filter
- ⬜ `deconvolve(signal, divisor)` — Deconvolution
- ⬜ `hilbert(x)` — Analytic signal via Hilbert transform
- ⬜ `hilbert2(x)` — 2-D Hilbert transform
- ⬜ `envelope(z)` — Envelope of analytic signal
- ⬜ `detrend(data, type?)` — Remove trend from data
- ⬜ `decimate(x, q)` — Downsample after anti-alias filter
- ⬜ `resample(x, num)` — Resample using Fourier method
- ⬜ `resample_poly(x, up, down)` — Resample using polyphase filter
- ⬜ `upfirdn(h, x, up?, down?)` — Upsample, FIR filter, downsample

### Priority Additions — Filter Design
- ⬜ `cheby1(N, rp, Wn, options?)` — Chebyshev type I filter
- ⬜ `cheby2(N, rs, Wn, options?)` — Chebyshev type II filter
- ⬜ `ellip(N, rp, rs, Wn, options?)` — Elliptic (Cauer) filter
- ⬜ `bessel(N, Wn, options?)` — Bessel/Thomson filter
- ⬜ `iirfilter(N, Wn, options?)` — IIR digital/analog filter design
- ⬜ `iirdesign(wp, ws, gpass, gstop)` — IIR filter from specs
- ⬜ `firwin2(numtaps, freq, gain)` — FIR filter (frequency sampling)
- ⬜ `firls(numtaps, bands, desired)` — FIR filter (least-squares)
- ⬜ `remez(numtaps, bands, desired)` — FIR filter (Parks-McClellan)
- ⬜ `kaiserord(ripple, width)` — Kaiser window FIR order estimation
- ⬜ `freqz(b, a?, worN?)` — Frequency response of digital filter
- ⬜ `sosfreqz(sos, worN?)` — Frequency response of SOS filter
- ⬜ `freqs(b, a, worN)` — Frequency response of analog filter
- ⬜ `bilinear(b, a, fs)` — Bilinear transformation
- ⬜ `bilinear_zpk(z, p, k, fs)` — Bilinear for zpk
- ⬜ `savgol_coeffs(window_length, polyorder)` — Savitzky-Golay coefficients

### Priority Additions — Spectral Analysis
- ⬜ `periodogram(x, options?)` — Periodogram PSD estimate
- ⬜ `csd(x, y, options?)` — Cross spectral density
- ⬜ `coherence(x, y, options?)` — Magnitude squared coherence
- ⬜ `stft(x, options?)` — Short-time Fourier transform
- ⬜ `istft(Zxx, options?)` — Inverse STFT

### Priority Additions — Peak Finding
- ⬜ `find_peaks(x, options?)` — Find peaks in signal
- ⬜ `peak_prominences(x, peaks)` — Peak prominences
- ⬜ `peak_widths(x, peaks)` — Peak widths

### Priority Additions — Windows
- ⬜ `get_window(window, Nx)` — Get window function by name
- ⬜ `hann(M)` — Hann window
- ⬜ `hamming(M)` — Hamming window
- ⬜ `blackman(M)` — Blackman window
- ⬜ `kaiser(M, beta)` — Kaiser window
- ⬜ `tukey(M, alpha?)` — Tukey window
- ⬜ `gaussian(M, std)` — Gaussian window
- ⬜ `bartlett(M)` — Bartlett window

---

## 7. `spatial` — Spatial Algorithms and Distance

### Current Stubs
- 🔲 `KDTree` — kd-tree (`query`, `query_ball_point`)
- 🔲 `Delaunay` — Delaunay tessellation
- 🔲 `ConvexHull` — Convex hull
- 🔲 `Voronoi` — Voronoi diagram
- 🔲 `distance.euclidean(u, v)` — Euclidean distance
- 🔲 `distance.cosine(u, v)` — Cosine distance
- 🔲 `distance.cdist(XA, XB, metric?)` — Pairwise distance (all pairs)
- 🔲 `distance.pdist(X, metric?)` — Condensed pairwise distance

### Priority Additions — Spatial Structures
- ⬜ `SphericalVoronoi(points)` — Spherical Voronoi diagram
- ⬜ `HalfspaceIntersection(halfspaces, interior_point)` — Halfspace intersection

### Priority Additions — Functions
- ⬜ `distance_matrix(x, y)` — Full distance matrix
- ⬜ `minkowski_distance(u, v, p)` — Minkowski distance
- ⬜ `procrustes(data1, data2)` — Procrustes analysis

### Priority Additions — Distance Metrics
- ⬜ `distance.cityblock(u, v)` — Manhattan distance
- ⬜ `distance.chebyshev(u, v)` — Chebyshev distance
- ⬜ `distance.minkowski(u, v, p)` — Minkowski distance
- ⬜ `distance.mahalanobis(u, v, VI)` — Mahalanobis distance
- ⬜ `distance.correlation(u, v)` — Correlation distance
- ⬜ `distance.hamming(u, v)` — Hamming distance
- ⬜ `distance.jaccard(u, v)` — Jaccard distance
- ⬜ `distance.braycurtis(u, v)` — Bray-Curtis distance
- ⬜ `distance.canberra(u, v)` — Canberra distance
- ⬜ `distance.sqeuclidean(u, v)` — Squared Euclidean distance
- ⬜ `distance.seuclidean(u, v, V)` — Standardized Euclidean distance
- ⬜ `distance.directed_hausdorff(u, v)` — Directed Hausdorff distance
- ⬜ `distance.squareform(X)` — Convert condensed ↔ square distance matrix

---

## 8. `sparse` — Sparse Matrices

### Current Stubs
- ✅ `csr_matrix(data)` — Compressed Sparse Row
- ✅ `csc_matrix(data)` — Compressed Sparse Column
- ✅ `eye(m, n?, k?)` — Sparse identity matrix
- ✅ `diags(diagonals, offsets?, shape?)` — Diagonal sparse matrix

### Priority Additions — Formats
- ✅ `coo_matrix(data)` / `coo_array(data)` — Coordinate format
- ⬜ `lil_matrix(shape)` / `lil_array(shape)` — List of lists (construction)
- ⬜ `bsr_matrix(data)` / `bsr_array(data)` — Block sparse row
- ⬜ `dok_matrix(shape)` / `dok_array(shape)` — Dictionary of keys
- ⬜ `dia_matrix(data)` / `dia_array(data)` — Diagonal format

### Priority Additions — Construction
- ⬜ `random(m, n, density?)` — Random sparse matrix
- ⬜ `kron(A, B)` — Kronecker product
- ⬜ `kronsum(A, B)` — Kronecker sum
- ⬜ `block_diag(mats)` — Block diagonal matrix
- ⬜ `hstack(blocks)` — Horizontal stack
- ⬜ `vstack(blocks)` — Vertical stack
- ⬜ `tril(A, k?)` — Lower triangle
- ⬜ `triu(A, k?)` — Upper triangle
- ⬜ `issparse(x)` — Check if sparse

### Priority Additions — `sparse.linalg`
- ⬜ `linalg.spsolve(A, b)` — Solve sparse system
- ⬜ `linalg.eigs(A, k?)` — Eigenvalues (sparse, largest)
- ⬜ `linalg.eigsh(A, k?)` — Eigenvalues (sparse, symmetric)
- ⬜ `linalg.svds(A, k?)` — SVD (sparse, truncated)
- ⬜ `linalg.inv(A)` — Sparse inverse
- ⬜ `linalg.norm(x)` — Sparse norm
- ⬜ `linalg.expm(A)` — Sparse matrix exponential
- ⬜ `linalg.cg(A, b)` — Conjugate gradient solver
- ⬜ `linalg.gmres(A, b)` — GMRES solver
- ⬜ `linalg.bicgstab(A, b)` — BiCGSTAB solver
- ⬜ `linalg.splu(A)` — Sparse LU decomposition
- ⬜ `linalg.spilu(A)` — Sparse incomplete LU
- ⬜ `linalg.LinearOperator` — Abstract linear operator

### Priority Additions — `sparse.csgraph`
- ⬜ `csgraph.shortest_path(csgraph)` — Shortest path (all algorithms)
- ⬜ `csgraph.dijkstra(csgraph)` — Dijkstra's algorithm
- ⬜ `csgraph.floyd_warshall(csgraph)` — Floyd-Warshall
- ⬜ `csgraph.bellman_ford(csgraph)` — Bellman-Ford
- ⬜ `csgraph.connected_components(csgraph)` — Connected components
- ⬜ `csgraph.laplacian(csgraph)` — Graph Laplacian
- ⬜ `csgraph.minimum_spanning_tree(csgraph)` — Minimum spanning tree
- ⬜ `csgraph.breadth_first_order(csgraph, i_start)` — BFS ordering
- ⬜ `csgraph.depth_first_order(csgraph, i_start)` — DFS ordering

---

## 9. `ndimage` — N-Dimensional Image Processing

### Current Stubs
- 🔲 `convolve(input, weights, options?)` — Multi-dimensional convolution
- 🔲 `gaussian_filter(input, sigma, options?)` — Gaussian filter
- 🔲 `label(input, structure?)` — Label connected features
- 🔲 `binary_erosion(input, structure?, iterations?)` — Binary erosion
- 🔲 `binary_dilation(input, structure?, iterations?)` — Binary dilation

### Priority Additions — Filters
- ⬜ `correlate(input, weights, options?)` — Multi-dimensional correlation
- ⬜ `uniform_filter(input, size?)` — Uniform (box) filter
- ⬜ `median_filter(input, size?)` — Median filter
- ⬜ `maximum_filter(input, size?)` — Maximum filter
- ⬜ `minimum_filter(input, size?)` — Minimum filter
- ⬜ `sobel(input, axis?)` — Sobel edge detection
- ⬜ `prewitt(input, axis?)` — Prewitt edge detection
- ⬜ `laplace(input)` — Laplacian filter
- ⬜ `gaussian_laplace(input, sigma)` — Gaussian Laplacian (LoG)
- ⬜ `gaussian_gradient_magnitude(input, sigma)` — Gaussian gradient magnitude
- ⬜ `generic_filter(input, function, size)` — Generic filter with callback
- ⬜ `rank_filter(input, rank, size)` — Rank filter
- ⬜ `percentile_filter(input, percentile, size)` — Percentile filter

### Priority Additions — Morphology
- ⬜ `binary_opening(input, structure?)` — Binary opening
- ⬜ `binary_closing(input, structure?)` — Binary closing
- ⬜ `binary_fill_holes(input, structure?)` — Fill holes in binary objects
- ⬜ `generate_binary_structure(rank, connectivity)` — Generate structuring element
- ⬜ `grey_erosion(input, size?)` — Greyscale erosion
- ⬜ `grey_dilation(input, size?)` — Greyscale dilation
- ⬜ `grey_opening(input, size?)` — Greyscale opening
- ⬜ `grey_closing(input, size?)` — Greyscale closing

### Priority Additions — Interpolation / Geometric
- ⬜ `zoom(input, zoom)` — Zoom (resize) array
- ⬜ `rotate(input, angle)` — Rotate array
- ⬜ `shift(input, shift)` — Shift array
- ⬜ `affine_transform(input, matrix)` — Affine transformation
- ⬜ `map_coordinates(input, coordinates)` — Map coordinates interpolation

### Priority Additions — Measurements
- ⬜ `center_of_mass(input, labels?)` — Center of mass
- ⬜ `find_objects(input)` — Find objects (bounding boxes)
- ⬜ `sum_labels(input, labels?)` — Sum by label
- ⬜ `mean(input, labels?)` — Mean by label
- ⬜ `variance(input, labels?)` — Variance by label
- ⬜ `standard_deviation(input, labels?)` — Std dev by label
- ⬜ `minimum(input, labels?)` — Minimum by label
- ⬜ `maximum(input, labels?)` — Maximum by label
- ⬜ `extrema(input, labels?)` — Min and max by label
- ⬜ `histogram(input, min, max, bins)` — Histogram by label

---

## 10. `cluster` — Clustering

### Current Stubs
- 🔲 `kmeans(data, k, options?)` — K-means clustering
- 🔲 `hierarchy.linkage(y, method?)` — Hierarchical clustering linkage
- 🔲 `hierarchy.fcluster(Z, t, criterion?)` — Form flat clusters

### Priority Additions — Vector Quantization
- ⬜ `vq.vq(obs, code_book)` — Assign codes from code book
- ⬜ `vq.whiten(obs)` — Normalize observations by std dev

### Priority Additions — Hierarchy
- ⬜ `hierarchy.dendrogram(Z)` — Generate dendrogram data
- ⬜ `hierarchy.cut_tree(Z, n_clusters?)` — Cut dendrogram
- ⬜ `hierarchy.leaves_list(Z)` — Leaf order
- ⬜ `hierarchy.optimal_leaf_ordering(Z, y)` — Optimal leaf ordering
- ⬜ `hierarchy.cophenet(Z)` — Cophenetic distances
- ⬜ `hierarchy.inconsistent(Z, d?)` — Inconsistency statistics
- ⬜ `hierarchy.maxdists(Z)` — Maximum distances
- ⬜ `hierarchy.ward(y)` — Ward linkage
- ⬜ `hierarchy.single(y)` — Single linkage
- ⬜ `hierarchy.complete(y)` — Complete linkage
- ⬜ `hierarchy.average(y)` — Average linkage
- ⬜ `hierarchy.weighted(y)` — Weighted linkage
- ⬜ `hierarchy.centroid(y)` — Centroid linkage
- ⬜ `hierarchy.median(y)` — Median linkage
- ⬜ `hierarchy.is_valid_linkage(Z)` — Validate linkage matrix
- ⬜ `hierarchy.is_monotonic(Z)` — Check monotonicity

---

## 11. `io` — Input/Output

### Current Stubs
- 🔲 `loadmat(fileOrBuffer, options?)` — Load MATLAB .mat files
- 🔲 `savemat(filename, mdict, options?)` — Save to MATLAB .mat format

### Priority Additions
- ⬜ `whosmat(filename)` — List variables in .mat file
- ⬜ `mmread(source)` — Read Matrix Market file
- ⬜ `mmwrite(target, a)` — Write Matrix Market file
- ⬜ `mminfo(source)` — Matrix Market file info
- ⬜ `wavfile.read(filename)` — Read WAV audio file
- ⬜ `wavfile.write(filename, rate, data)` — Write WAV audio file
- ⬜ `FortranFile(filename, mode?)` — Read/write Fortran unformatted files
- ⬜ `hb_read(source)` — Read Harwell-Boeing sparse format
- ⬜ `hb_write(target, m)` — Write Harwell-Boeing sparse format
- ⬜ `arff.loadarff(filename)` — Read ARFF file

---

## 12. `constants` — Physical and Mathematical Constants

### Currently Implemented
- ✅ `c` / `speed_of_light` — 299792458
- ✅ `h` / `Planck` — 6.62607015e-34
- ✅ `hbar` — 1.054571817e-34
- ✅ `G` — 6.67430e-11
- ✅ `g` — 9.80665
- ✅ `e` — 1.602176634e-19
- ✅ `k` / `Boltzmann` — 1.380649e-23
- ✅ `N_A` — 6.02214076e23
- ✅ `R` — 8.314462618
- ✅ `sigma` — 5.670374419e-8

### Current Stubs
- 🔲 `physical_constants(name)` — Look up constant by name

### Priority Additions — Mathematical Constants
- ⬜ `pi` — 3.14159265358979...
- ⬜ `golden` / `golden_ratio` — 1.61803398874989...

### Priority Additions — Physical Constants
- ⬜ `mu_0` — Vacuum magnetic permeability
- ⬜ `epsilon_0` — Vacuum electric permittivity
- ⬜ `alpha` / `fine_structure` — Fine-structure constant
- ⬜ `Wien` — Wien displacement law constant
- ⬜ `Rydberg` — Rydberg constant
- ⬜ `m_e` / `electron_mass` — Electron mass
- ⬜ `m_p` / `proton_mass` — Proton mass
- ⬜ `m_n` / `neutron_mass` — Neutron mass
- ⬜ `u` / `atomic_mass` — Atomic mass constant
- ⬜ `eV` / `electron_volt` — Electron volt (in joules)

### Priority Additions — Database Functions
- ⬜ `value(name)` — Get constant value by name
- ⬜ `unit(name)` — Get constant unit by name
- ⬜ `precision(name)` — Get constant precision by name
- ⬜ `find(sub)` — Search constants by substring

### Priority Additions — Unit Conversions
- ⬜ `convert_temperature(val, old_scale, new_scale)` — Temperature conversion
- ⬜ SI prefixes: `yotta`, `zetta`, `exa`, `peta`, `tera`, `giga`, `mega`, `kilo`, `hecto`, `deka`, `deci`, `centi`, `milli`, `micro`, `nano`, `pico`, `femto`, `atto`, `zepto`, `yocto`
- ⬜ Length: `inch`, `foot`, `yard`, `mile`, `mil`, `pt`, `point`, `survey_foot`, `survey_mile`, `nautical_mile`, `fermi`, `angstrom`, `micron`, `au`, `astronomical_unit`, `light_year`, `parsec`
- ⬜ Mass: `gram`, `metric_ton`, `grain`, `lb`, `pound`, `blob`, `slinch`, `slug`, `oz`, `ounce`, `stone`, `long_ton`, `short_ton`, `troy_ounce`, `troy_pound`, `carat`
- ⬜ Time: `minute`, `hour`, `day`, `week`, `year`, `Julian_year`
- ⬜ Pressure: `atm`, `atmosphere`, `bar`, `torr`, `mmHg`, `psi`
- ⬜ Energy: `calorie`, `calorie_th`, `calorie_IT`, `erg`, `Btu`, `Btu_IT`, `Btu_th`, `ton_TNT`
- ⬜ Power: `hp`, `horsepower`
- ⬜ Temperature: `zero_Celsius`, `degree_Fahrenheit`

---

## Summary

| Module | ✅ Done | 🔲 Stubbed | ⬜ To Create | Total |
|--------|---------|-----------|-------------|-------|
| special | 0 | 10 | ~35 | ~45 |
| stats | 4 | 9 | ~37 | ~50 |
| optimize | 3 | 4 | ~18 | ~25 |
| integrate | 1 | 5 | ~9 | ~15 |
| interpolate | 0 | 4 | ~14 | ~18 |
| signal | 0 | 7 | ~35 | ~42 |
| spatial | 0 | 8 | ~15 | ~23 |
| sparse | 5 | 0 | ~29 | ~34 |
| ndimage | 0 | 5 | ~25 | ~30 |
| cluster | 0 | 3 | ~15 | ~18 |
| io | 0 | 2 | ~10 | ~12 |
| constants | 10 | 1 | ~50+ | ~61 |
| **Total** | **22** | **59** | **~292** | **~373** |
