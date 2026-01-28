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
