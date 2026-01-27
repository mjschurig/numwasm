/**
 * Statistical distributions and tests.
 * @module stats
 */

import { NotImplementedError } from '../errors.js';

/** Descriptive statistics result. */
export interface DescribeResult {
  nobs: number;
  minmax: [number, number];
  mean: number;
  variance: number;
  skewness: number;
  kurtosis: number;
}

/** Result of a statistical test. */
export interface TestResult {
  statistic: number;
  pvalue: number;
}

/** A continuous probability distribution. */
export interface ContinuousDistribution {
  pdf(x: number): number;
  cdf(x: number): number;
  ppf(q: number): number;
  rvs(size?: number): number | number[];
  mean(): number;
  std(): number;
}

/**
 * Compute descriptive statistics.
 * Mirrors scipy.stats.describe.
 */
export function describe(_a: number[]): DescribeResult {
  throw new NotImplementedError('sciwasm.stats.describe');
}

/**
 * Normal (Gaussian) distribution.
 * Mirrors scipy.stats.norm.
 */
export function norm(_loc?: number, _scale?: number): ContinuousDistribution {
  throw new NotImplementedError('sciwasm.stats.norm');
}

/**
 * Student's t distribution.
 * Mirrors scipy.stats.t.
 */
export function t(_df: number, _loc?: number, _scale?: number): ContinuousDistribution {
  throw new NotImplementedError('sciwasm.stats.t');
}

/**
 * F distribution.
 * Mirrors scipy.stats.f.
 */
export function f(_dfn: number, _dfd: number, _loc?: number, _scale?: number): ContinuousDistribution {
  throw new NotImplementedError('sciwasm.stats.f');
}

/**
 * Chi-squared distribution.
 * Mirrors scipy.stats.chi2.
 */
export function chi2(_df: number, _loc?: number, _scale?: number): ContinuousDistribution {
  throw new NotImplementedError('sciwasm.stats.chi2');
}

/**
 * Pearson correlation coefficient.
 * Mirrors scipy.stats.pearsonr.
 */
export function pearsonr(_x: number[], _y: number[]): TestResult {
  throw new NotImplementedError('sciwasm.stats.pearsonr');
}

/**
 * Spearman rank-order correlation coefficient.
 * Mirrors scipy.stats.spearmanr.
 */
export function spearmanr(_a: number[], _b?: number[]): TestResult {
  throw new NotImplementedError('sciwasm.stats.spearmanr');
}

/**
 * Independent two-sample t-test.
 * Mirrors scipy.stats.ttest_ind.
 */
export function ttest_ind(_a: number[], _b: number[], _options?: { equal_var?: boolean }): TestResult {
  throw new NotImplementedError('sciwasm.stats.ttest_ind');
}

/**
 * One-sample t-test.
 * Mirrors scipy.stats.ttest_1samp.
 */
export function ttest_1samp(_a: number[], _popmean: number): TestResult {
  throw new NotImplementedError('sciwasm.stats.ttest_1samp');
}

/**
 * Kolmogorov-Smirnov test.
 * Mirrors scipy.stats.kstest.
 */
export function kstest(_rvs: number[], _cdf: string | ((x: number) => number)): TestResult {
  throw new NotImplementedError('sciwasm.stats.kstest');
}
