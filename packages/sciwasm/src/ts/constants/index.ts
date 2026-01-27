/**
 * Physical and mathematical constants.
 * @module constants
 */

import { NotImplementedError } from '../errors.js';

/** Speed of light in vacuum (m/s). */
export const c = 299792458;

/** Planck constant (J·s). */
export const h = 6.62607015e-34;

/** Reduced Planck constant (J·s). */
export const hbar = 1.054571817e-34;

/** Newtonian constant of gravitation (m³/(kg·s²)). */
export const G = 6.67430e-11;

/** Standard acceleration of gravity (m/s²). */
export const g = 9.80665;

/** Elementary charge (C). */
export const e = 1.602176634e-19;

/** Boltzmann constant (J/K). */
export const k = 1.380649e-23;

/** Avogadro constant (1/mol). */
export const N_A = 6.02214076e23;

/** Molar gas constant (J/(mol·K)). */
export const R = 8.314462618;

/** Stefan-Boltzmann constant (W/(m²·K⁴)). */
export const sigma = 5.670374419e-8;

/** Speed of light alias. */
export const speed_of_light = c;

/** Planck constant alias. */
export const Planck = h;

/** Boltzmann constant alias. */
export const Boltzmann = k;

/**
 * Look up a physical constant by name.
 * Mirrors scipy.constants.physical_constants.
 */
export function physical_constants(
  _name: string
): { value: number; unit: string; uncertainty: number } {
  throw new NotImplementedError('sciwasm.constants.physical_constants');
}
