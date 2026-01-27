import { describe, it, expect } from 'vitest';
import * as sci from '../../src/ts/index';
import { NotImplementedError } from '../../src/ts/errors';

describe('sciwasm module exports', () => {
  it('exports all 12 modules', () => {
    expect(sci.optimize).toBeDefined();
    expect(sci.integrate).toBeDefined();
    expect(sci.interpolate).toBeDefined();
    expect(sci.stats).toBeDefined();
    expect(sci.signal).toBeDefined();
    expect(sci.spatial).toBeDefined();
    expect(sci.special).toBeDefined();
    expect(sci.sparse).toBeDefined();
    expect(sci.cluster).toBeDefined();
    expect(sci.io).toBeDefined();
    expect(sci.ndimage).toBeDefined();
    expect(sci.constants).toBeDefined();
  });

  it('exports NotImplementedError', () => {
    expect(sci.NotImplementedError).toBeDefined();
  });

  it('optimize.minimize throws NotImplementedError', () => {
    expect(() => sci.optimize.minimize(() => 0, [0])).toThrow(NotImplementedError);
  });

  it('integrate.quad throws NotImplementedError', () => {
    expect(() => sci.integrate.quad(Math.sin, 0, 1)).toThrow(NotImplementedError);
  });

  it('stats.describe throws NotImplementedError', () => {
    expect(() => sci.stats.describe([1, 2, 3])).toThrow(NotImplementedError);
  });

  it('special.gamma throws NotImplementedError', () => {
    expect(() => sci.special.gamma(1)).toThrow(NotImplementedError);
  });

  it('constants exports numeric values', () => {
    expect(sci.constants.c).toBe(299792458);
    expect(sci.constants.h).toBeCloseTo(6.626e-34, 36);
    expect(sci.constants.k).toBeCloseTo(1.381e-23, 26);
  });

  it('spatial.KDTree constructor throws NotImplementedError', () => {
    expect(() => new sci.spatial.KDTree([[0, 0]])).toThrow(NotImplementedError);
  });

  it('signal.convolve throws NotImplementedError', () => {
    expect(() => sci.signal.convolve([1], [1])).toThrow(NotImplementedError);
  });

  it('sparse.csr_matrix throws NotImplementedError', () => {
    expect(() => sci.sparse.csr_matrix([[1, 0], [0, 1]])).toThrow(NotImplementedError);
  });

  it('cluster.kmeans throws NotImplementedError', () => {
    expect(() => sci.cluster.kmeans([[1, 2]], 1)).toThrow(NotImplementedError);
  });

  it('io.loadmat throws NotImplementedError', () => {
    expect(() => sci.io.loadmat('test.mat')).toThrow(NotImplementedError);
  });

  it('ndimage.convolve throws NotImplementedError', () => {
    expect(() => sci.ndimage.convolve([[1]], [[1]])).toThrow(NotImplementedError);
  });

  it('interpolate.interp1d throws NotImplementedError', () => {
    expect(() => sci.interpolate.interp1d([1, 2], [3, 4])).toThrow(NotImplementedError);
  });
});
