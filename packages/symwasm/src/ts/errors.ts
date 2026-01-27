/**
 * Error thrown when calling a function that has not been implemented yet.
 */
export class NotImplementedError extends Error {
  constructor(name: string) {
    super(name + " is not yet implemented. Coming soon in symwasm.");
    this.name = "NotImplementedError";
  }
}
