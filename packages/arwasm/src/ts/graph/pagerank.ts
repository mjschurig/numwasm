/**
 * PAGERANK - PageRank Algorithm
 *
 * Computes PageRank scores for nodes in a directed graph using
 * the power iteration method.
 *
 * PageRank models a random surfer who:
 * - With probability d, follows a random outgoing link
 * - With probability (1-d), jumps to a random page
 *
 * The stationary distribution of this Markov chain gives the PageRank scores.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for PageRank computation.
 */
export interface PagerankOptions {
  /**
   * Damping factor (probability of following a link).
   * Standard value is 0.85.
   * @default 0.85
   */
  damping?: number;

  /**
   * Convergence tolerance.
   * @default 1e-8
   */
  tol?: number;

  /**
   * Maximum number of iterations.
   * @default 100
   */
  maxiter?: number;

  /**
   * Personalization vector (non-uniform teleportation probabilities).
   * If provided, must sum to 1.
   * @default uniform distribution
   */
  personalization?: Float64Array;
}

/**
 * Result from PageRank computation.
 */
export interface PagerankResult {
  /**
   * PageRank scores for each node (sums to 1).
   */
  scores: Float64Array;

  /**
   * Number of iterations performed.
   */
  niter: number;

  /**
   * Whether the algorithm converged.
   */
  converged: boolean;

  /**
   * Final residual norm.
   */
  residual: number;
}

/**
 * Compute PageRank scores for a directed graph.
 *
 * Uses power iteration to find the dominant eigenvector of the
 * Google matrix G = d*M + (1-d)*v*e^T, where:
 * - M is the column-stochastic transition matrix
 * - d is the damping factor
 * - v is the personalization vector
 * - e is the all-ones vector
 *
 * @param outlinksMatvec - Function computing y = M*x where M is column-stochastic
 *                         (columns sum to 1, or 0 for dangling nodes)
 * @param n - Number of nodes in the graph
 * @param options - Algorithm options
 * @returns PageRank scores
 *
 * @example
 * ```ts
 * import { pagerank } from 'arwasm';
 *
 * // For a directed graph with adjacency matrix
 * const n = 1000;
 *
 * // Compute out-degrees for normalization
 * const outDegrees = new Float64Array(n);
 * // ... fill with out-degrees
 *
 * // Create column-stochastic transition matrix matvec
 * const outlinksMatvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   // For each edge (i,j): y[j] += x[i] / outDegree[i]
 *   // ... implement sparse matvec
 *   return y;
 * };
 *
 * const result = await pagerank(outlinksMatvec, n);
 * console.log('Top PageRank node:', result.scores.indexOf(Math.max(...result.scores)));
 * ```
 */
export async function pagerank(
  outlinksMatvec: MatVecFunction,
  n: number,
  options?: PagerankOptions
): Promise<PagerankResult> {
  const {
    damping = 0.85,
    tol = 1e-8,
    maxiter = 100,
    personalization,
  } = options ?? {};

  // Validate damping factor
  if (damping <= 0 || damping >= 1) {
    throw new Error('Damping factor must be in (0, 1)');
  }

  // Initialize personalization vector (uniform if not provided)
  const teleport = personalization ?? new Float64Array(n).fill(1 / n);

  // Validate personalization sums to ~1
  const teleportSum = teleport.reduce((a, b) => a + b, 0);
  if (Math.abs(teleportSum - 1) > 1e-6) {
    throw new Error('Personalization vector must sum to 1');
  }

  // Initialize PageRank vector uniformly
  let pr = new Float64Array(n).fill(1 / n);
  let converged = false;
  let niter = 0;
  let residual = Infinity;

  // Power iteration
  for (niter = 1; niter <= maxiter; niter++) {
    // Compute M*pr (follow links)
    const Mpr = outlinksMatvec(pr);
    const MprArr = Mpr instanceof Float64Array ? Mpr : new Float64Array(Mpr);

    // Handle dangling nodes: their probability mass goes to teleport
    // Dangling mass = sum of pr[i] for dangling nodes
    // We detect this as 1 - sum(M*pr) since M is column-stochastic
    const linkMass = MprArr.reduce((a, b) => a + b, 0);
    const danglingMass = 1 - linkMass;

    // New PageRank: pr' = d * M * pr + d * dangling_mass * teleport + (1-d) * teleport
    //                   = d * M * pr + (d * dangling_mass + (1-d)) * teleport
    const teleportCoeff = damping * danglingMass + (1 - damping);

    const newPr = new Float64Array(n);
    residual = 0;

    for (let i = 0; i < n; i++) {
      newPr[i] = damping * MprArr[i] + teleportCoeff * teleport[i];
      residual += Math.abs(newPr[i] - pr[i]);
    }

    // Normalize (should already sum to 1, but ensure numerical stability)
    const sum = newPr.reduce((a, b) => a + b, 0);
    for (let i = 0; i < n; i++) {
      newPr[i] /= sum;
    }

    pr = newPr;

    // Check convergence
    if (residual < tol) {
      converged = true;
      break;
    }
  }

  return {
    scores: pr,
    niter,
    converged,
    residual,
  };
}

/**
 * Alternative PageRank using ARPACK eigenvalue solver.
 *
 * Finds the dominant eigenvector of the Google matrix directly.
 * This may be more efficient for very large graphs or when high
 * precision is needed.
 *
 * @param outlinksMatvec - Function computing y = M*x
 * @param n - Number of nodes
 * @param options - Algorithm options
 * @returns PageRank scores
 */
export async function pagerankEigs(
  outlinksMatvec: MatVecFunction,
  n: number,
  options?: PagerankOptions
): Promise<PagerankResult> {
  const {
    damping = 0.85,
    tol = 1e-10,
    maxiter = n * 10,
    personalization,
  } = options ?? {};

  const teleport = personalization ?? new Float64Array(n).fill(1 / n);

  // Create Google matrix matvec: G*x = d*M*x + (1-d)*v*(e^T*x)
  const googleMatvec: MatVecFunction = (x: Float64Array) => {
    // Compute M*x
    const Mx = outlinksMatvec(x);
    const MxArr = Mx instanceof Float64Array ? Mx : new Float64Array(Mx);

    // Compute e^T*x = sum(x)
    const xSum = x.reduce((a, b) => a + b, 0);

    // Handle dangling nodes
    const linkMass = MxArr.reduce((a, b) => a + b, 0);
    const danglingMass = xSum - linkMass;

    // G*x = d*M*x + (d*dangling + (1-d))*v*sum(x)
    const teleportCoeff = (damping * danglingMass + (1 - damping) * xSum);

    const result = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      result[i] = damping * MxArr[i] + teleportCoeff * teleport[i];
    }
    return result;
  };

  // Find dominant eigenvector (eigenvalue should be 1)
  const result = await eigs(googleMatvec, n, 1, {
    which: 'LM',
    tol,
    maxiter,
    return_eigenvectors: true,
  });

  if (!result.success || !result.eigenvectors || result.eigenvectors.length === 0) {
    return {
      scores: new Float64Array(n).fill(1 / n),
      niter: result.niter,
      converged: false,
      residual: Infinity,
    };
  }

  // Normalize eigenvector to sum to 1 (and ensure non-negative)
  const ev = result.eigenvectors[0];
  let sum = 0;
  for (let i = 0; i < n; i++) {
    ev[i] = Math.abs(ev[i]); // Eigenvector may have arbitrary sign
    sum += ev[i];
  }
  for (let i = 0; i < n; i++) {
    ev[i] /= sum;
  }

  return {
    scores: ev,
    niter: result.niter,
    converged: result.success,
    residual: 0, // Not computed for eigenvector method
  };
}
