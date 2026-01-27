/**
 * Signal processing functions.
 * @module signal
 */

import { NotImplementedError } from '../errors.js';

/**
 * Convolve two 1-D arrays.
 * Mirrors scipy.signal.convolve.
 */
export function convolve(
  _in1: number[],
  _in2: number[],
  _mode?: 'full' | 'valid' | 'same'
): number[] {
  throw new NotImplementedError('sciwasm.signal.convolve');
}

/**
 * Convolve two N-D arrays using FFT.
 * Mirrors scipy.signal.fftconvolve.
 */
export function fftconvolve(
  _in1: number[],
  _in2: number[],
  _mode?: 'full' | 'valid' | 'same'
): number[] {
  throw new NotImplementedError('sciwasm.signal.fftconvolve');
}

/**
 * Butterworth digital and analog filter design.
 * Mirrors scipy.signal.butter.
 */
export function butter(
  _N: number,
  _Wn: number | [number, number],
  _options?: { btype?: 'low' | 'high' | 'band' | 'bandstop'; output?: 'ba' | 'sos' }
): { sos: number[][] } | { b: number[]; a: number[] } {
  throw new NotImplementedError('sciwasm.signal.butter');
}

/**
 * Filter data along one dimension using cascaded second-order sections.
 * Mirrors scipy.signal.sosfilt.
 */
export function sosfilt(_sos: number[][], _x: number[]): number[] {
  throw new NotImplementedError('sciwasm.signal.sosfilt');
}

/**
 * FIR filter design using the window method.
 * Mirrors scipy.signal.firwin.
 */
export function firwin(
  _numtaps: number,
  _cutoff: number | number[],
  _options?: { window?: string; pass_zero?: boolean; fs?: number }
): number[] {
  throw new NotImplementedError('sciwasm.signal.firwin');
}

/**
 * Estimate power spectral density using Welch's method.
 * Mirrors scipy.signal.welch.
 */
export function welch(
  _x: number[],
  _options?: { fs?: number; window?: string; nperseg?: number }
): { f: number[]; Pxx: number[] } {
  throw new NotImplementedError('sciwasm.signal.welch');
}

/**
 * Compute a spectrogram.
 * Mirrors scipy.signal.spectrogram.
 */
export function spectrogram(
  _x: number[],
  _options?: { fs?: number; window?: string; nperseg?: number }
): { f: number[]; t: number[]; Sxx: number[][] } {
  throw new NotImplementedError('sciwasm.signal.spectrogram');
}
