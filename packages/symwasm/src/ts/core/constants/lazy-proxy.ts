/**
 * Helper for creating lazy constant proxies.
 * @module core/constants/lazy-proxy
 * @internal
 */

import { Expr } from '../expr.js';

/**
 * Helper to create a lazy proxy for a constant.
 * The proxy delegates all property access to the lazily-initialized constant.
 */
export function lazyConstantProxy<T extends Expr>(getter: () => T): T {
  return new Proxy({} as T, {
    get(_target, prop) {
      const instance = getter();
      const value = (instance as any)[prop];
      if (typeof value === 'function') {
        return value.bind(instance);
      }
      return value;
    },
  });
}
