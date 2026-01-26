/**
 * Error classes for the strings module.
 */

/**
 * ValueError - raised when a substring is not found in index/rindex operations.
 * Follows NumPy's convention.
 */
export class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
