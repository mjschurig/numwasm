/**
 * Error thrown when calling a function that has not been implemented yet.
 */
export class NotImplementedError extends Error {
  constructor(name: string) {
    super(name + " is not yet implemented. Coming soon in sciwasm.");
    this.name = "NotImplementedError";
  }
}

/**
 * Warning issued when numerical integration encounters issues.
 */
export class IntegrationWarning extends Error {
  constructor(message: string) {
    super(message);
    this.name = "IntegrationWarning";
  }
}

/**
 * Error thrown when a matrix is singular and cannot be factorized or inverted.
 */
export class SingularMatrixError extends Error {
  /** The column index where zero pivot was detected (if available) */
  readonly column?: number;

  constructor(column?: number) {
    super(
      column !== undefined
        ? `Matrix is singular: zero pivot at column ${column}`
        : 'Matrix is singular'
    );
    this.name = 'SingularMatrixError';
    this.column = column;
  }
}

/**
 * Error thrown for dimension mismatches in linear algebra operations.
 */
export class DimensionMismatchError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'DimensionMismatchError';
  }
}

/**
 * Error thrown by ARPACK eigenvalue solvers.
 */
export class ARPACKError extends Error {
  /** ARPACK error code */
  readonly code: number;

  constructor(code: number, message: string) {
    super(`ARPACK error ${code}: ${message}`);
    this.name = 'ARPACKError';
    this.code = code;
  }
}

/**
 * Error thrown when ARPACK does not converge.
 */
export class ARPACKNoConvergence extends Error {
  /** Number of eigenvalues that did converge */
  readonly nconv: number;

  constructor(nconv: number, k: number) {
    super(`ARPACK did not converge: only ${nconv} of ${k} eigenvalues converged`);
    this.name = 'ARPACKNoConvergence';
    this.nconv = nconv;
  }
}
