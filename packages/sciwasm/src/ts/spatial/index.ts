/**
 * Spatial data structures and algorithms.
 * @module spatial
 */

import { NotImplementedError } from '../errors.js';

/**
 * kd-tree for quick nearest-neighbor lookup.
 * Mirrors scipy.spatial.KDTree.
 */
export class KDTree {
  constructor(_data: number[][]) {
    throw new NotImplementedError('sciwasm.spatial.KDTree');
  }

  /** Query the tree for nearest neighbors. */
  query(_x: number[] | number[][], _k?: number): { distances: number[]; indices: number[] } {
    throw new NotImplementedError('sciwasm.spatial.KDTree.query');
  }

  /** Find all points within distance r. */
  query_ball_point(_x: number[], _r: number): number[] {
    throw new NotImplementedError('sciwasm.spatial.KDTree.query_ball_point');
  }
}

/**
 * Delaunay tessellation in N dimensions.
 * Mirrors scipy.spatial.Delaunay.
 */
export class Delaunay {
  /** Indices of simplices. */
  simplices: number[][] = [];

  constructor(_points: number[][]) {
    throw new NotImplementedError('sciwasm.spatial.Delaunay');
  }
}

/**
 * Convex hull in N dimensions.
 * Mirrors scipy.spatial.ConvexHull.
 */
export class ConvexHull {
  /** Indices of vertices. */
  vertices: number[] = [];
  /** Area (2D) or surface area (3D). */
  area: number = 0;
  /** Volume (only for 3D+). */
  volume: number = 0;

  constructor(_points: number[][]) {
    throw new NotImplementedError('sciwasm.spatial.ConvexHull');
  }
}

/**
 * Voronoi diagram in N dimensions.
 * Mirrors scipy.spatial.Voronoi.
 */
export class Voronoi {
  vertices: number[][] = [];
  regions: number[][] = [];

  constructor(_points: number[][]) {
    throw new NotImplementedError('sciwasm.spatial.Voronoi');
  }
}

/**
 * Pairwise distance computation.
 * Mirrors scipy.spatial.distance module.
 */
export const distance = {
  /** Euclidean distance between two 1-D arrays. */
  euclidean(_u: number[], _v: number[]): number {
    throw new NotImplementedError('sciwasm.spatial.distance.euclidean');
  },

  /** Cosine distance. */
  cosine(_u: number[], _v: number[]): number {
    throw new NotImplementedError('sciwasm.spatial.distance.cosine');
  },

  /** Compute distance matrix between each pair of row vectors. */
  cdist(_XA: number[][], _XB: number[][], _metric?: string): number[][] {
    throw new NotImplementedError('sciwasm.spatial.distance.cdist');
  },

  /** Compute condensed distance matrix. */
  pdist(_X: number[][], _metric?: string): number[] {
    throw new NotImplementedError('sciwasm.spatial.distance.pdist');
  },
};
