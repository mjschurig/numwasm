/**
 * Standard test functions for optimization algorithms.
 * These match the test functions used in scipy.optimize tests.
 */

/**
 * Rosenbrock function: f(x) = sum(100*(x[i+1] - x[i]^2)^2 + (1-x[i])^2)
 * Minimum at x = [1, 1, ..., 1], f(x*) = 0
 */
export function rosenbrock(x: number[]): number {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++) {
    sum += 100 * (x[i + 1] - x[i] * x[i]) ** 2 + (1 - x[i]) ** 2;
  }
  return sum;
}

/**
 * Beale function (2D): f(x,y) = (1.5-x+xy)^2 + (2.25-x+xy^2)^2 + (2.625-x+xy^3)^2
 * Minimum at (3, 0.5), f(x*) = 0
 */
export function beale(x: number[]): number {
  const [a, b] = x;
  return (
    (1.5 - a + a * b) ** 2 +
    (2.25 - a + a * b * b) ** 2 +
    (2.625 - a + a * b * b * b) ** 2
  );
}

/**
 * Booth function (2D): f(x,y) = (x + 2y - 7)^2 + (2x + y - 5)^2
 * Minimum at (1, 3), f(x*) = 0
 */
export function booth(x: number[]): number {
  return (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2;
}

/**
 * Simple quadratic: f(x) = sum(x[i]^2)
 * Minimum at x = [0, 0, ..., 0], f(x*) = 0
 */
export function quadratic(x: number[]): number {
  let sum = 0;
  for (let i = 0; i < x.length; i++) {
    sum += x[i] * x[i];
  }
  return sum;
}

/**
 * Himmelblau function (2D): f(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
 * Has 4 identical minima, one at (3, 2), f(x*) = 0
 */
export function himmelblau(x: number[]): number {
  return (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2;
}
