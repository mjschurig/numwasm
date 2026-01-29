import { describe, it, expect, beforeAll } from 'vitest';
import { spatial, loadWasmModule } from '../../../src/ts/index.js';

beforeAll(async () => {
  await loadWasmModule();
});

describe('KDTree', () => {
  describe('Construction', () => {
    it('constructs tree from 2D data', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);
      expect(tree).toBeDefined();
      expect(tree.size).toBe(3);
      expect(tree.ndim).toBe(2);
      tree.dispose();
    });

    it('constructs tree from 3D data', () => {
      const data = [[0, 0, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3]];
      const tree = new spatial.KDTree(data);
      expect(tree.size).toBe(4);
      expect(tree.ndim).toBe(3);
      tree.dispose();
    });

    it('throws error for empty data', () => {
      expect(() => new spatial.KDTree([])).toThrow('non-empty');
    });

    it('throws error for inconsistent dimensionality', () => {
      const data = [[0, 0], [1, 1, 1], [2, 2]];  // Wrong!
      expect(() => new spatial.KDTree(data)).toThrow('dimensionality');
    });

    it('accepts custom leafsize', () => {
      const data = [[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]];
      const tree = new spatial.KDTree(data, 2);  // leafsize=2
      expect(tree).toBeDefined();
      tree.dispose();
    });
  });

  describe('query() - Single nearest neighbor', () => {
    it('finds single nearest neighbor in 2D', () => {
      const data = [[0, 0], [1, 1], [2, 2], [3, 3]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1.1, 1.1], 1);
      expect(result.indices).toEqual([1]);
      expect(result.distances[0]).toBeCloseTo(0.141, 2);

      tree.dispose();
    });

    it('finds nearest neighbor at exact data point', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1, 1], 1);
      expect(result.indices).toEqual([1]);
      expect(result.distances[0]).toBeCloseTo(0, 10);

      tree.dispose();
    });

    it('finds nearest neighbor in 1D', () => {
      const data = [[0], [1], [2], [3], [4]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([2.3], 1);
      expect(result.indices).toEqual([2]);
      expect(result.distances[0]).toBeCloseTo(0.3, 10);

      tree.dispose();
    });

    it('finds nearest neighbor in 3D', () => {
      const data = [[0, 0, 0], [1, 1, 1], [2, 2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1.1, 0.9, 1.0], 1);
      expect(result.indices).toEqual([1]);
      expect(result.distances[0]).toBeCloseTo(0.141, 2);

      tree.dispose();
    });
  });

  describe('query() - Multiple nearest neighbors', () => {
    it('finds k nearest neighbors (k=3)', () => {
      const data = [[0, 0], [1, 1], [2, 2], [3, 3]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1.1, 1.1], 3);
      expect(result.indices).toHaveLength(3);
      expect(result.indices[0]).toBe(1);  // Closest
      // Other two should be points 0 and 2
      expect(result.indices).toContain(0);
      expect(result.indices).toContain(2);

      // Distances should be sorted (ascending)
      expect(result.distances[0]).toBeLessThan(result.distances[1]);
      expect(result.distances[1]).toBeLessThan(result.distances[2]);

      tree.dispose();
    });

    it('finds all neighbors when k >= n', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1, 1], 10);  // k > n
      expect(result.indices).toHaveLength(10);
      // First 3 should be real points
      expect(result.indices.slice(0, 3).sort()).toEqual([0, 1, 2]);
      // Rest should be placeholder (-1 or tree.size)

      tree.dispose();
    });

    it('handles k=1 correctly', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([0.1, 0.1], 1);
      expect(result.indices).toEqual([0]);

      tree.dispose();
    });
  });

  describe('query() - Batch queries', () => {
    it('handles batch queries (multiple points)', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([[0, 0], [1, 1], [2, 2]], 1);
      expect(result.indices).toEqual([[0], [1], [2]]);
      expect(result.distances).toHaveLength(3);
      expect(result.distances[0][0]).toBeCloseTo(0, 10);
      expect(result.distances[1][0]).toBeCloseTo(0, 10);
      expect(result.distances[2][0]).toBeCloseTo(0, 10);

      tree.dispose();
    });

    it('handles batch queries with k > 1', () => {
      const data = [[0, 0], [1, 1], [2, 2], [3, 3]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([[0.5, 0.5], [2.5, 2.5]], 2);
      expect(result.indices).toHaveLength(2);
      expect(result.indices[0]).toHaveLength(2);
      expect(result.indices[1]).toHaveLength(2);

      tree.dispose();
    });
  });

  describe('query() - Options', () => {
    it('respects distance_upper_bound option', () => {
      const data = [[0, 0], [10, 10]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([0, 0], 2, { distance_upper_bound: 5 });
      expect(result.indices[0]).toBe(0);  // First neighbor within bound
      // Second neighbor should be beyond bound (index -1 or infinity distance)
      expect(result.distances[1]).toBeGreaterThan(5);

      tree.dispose();
    });

    it('uses Euclidean distance by default (p=2)', () => {
      const data = [[0, 0], [1, 0], [0, 1]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([0, 0], 2);
      // Both [1,0] and [0,1] are equidistant under L2 norm (distance = 1.0)
      expect(result.distances[0]).toBeCloseTo(1.0, 10);
      expect(result.distances[1]).toBeCloseTo(1.0, 10);

      tree.dispose();
    });

    it('throws error for query with wrong dimensionality', () => {
      const data = [[0, 0], [1, 1]];
      const tree = new spatial.KDTree(data);

      expect(() => tree.query([0, 0, 0], 1)).toThrow('dimensionality');

      tree.dispose();
    });
  });

  describe('query_ball_point() - Single radius query', () => {
    it('finds neighbors within radius', () => {
      const data = [[0, 0], [1, 0], [2, 0], [3, 0]];
      const tree = new spatial.KDTree(data);

      const neighbors = tree.query_ball_point([1.5, 0], 0.6);
      expect(neighbors).toHaveLength(1);  // Single query returns array with one result
      expect(neighbors[0].sort()).toEqual([1, 2]);

      tree.dispose();
    });

    it('finds all neighbors in large radius', () => {
      const data = [[0, 0], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const neighbors = tree.query_ball_point([1, 1], 10);
      expect(neighbors).toHaveLength(1);
      expect(neighbors[0].sort()).toEqual([0, 1, 2]);

      tree.dispose();
    });

    it('returns empty array when no neighbors in radius', () => {
      const data = [[0, 0], [10, 10]];
      const tree = new spatial.KDTree(data);

      const neighbors = tree.query_ball_point([5, 5], 1.0);
      expect(neighbors).toHaveLength(1);
      expect(neighbors[0]).toEqual([]);

      tree.dispose();
    });

    it('includes exact point at radius boundary', () => {
      const data = [[0, 0], [1, 0], [2, 0]];
      const tree = new spatial.KDTree(data);

      // Query at [0.5, 0] with radius 0.5 should include [0,0] and [1,0]
      const neighbors = tree.query_ball_point([0.5, 0], 0.5);
      expect(neighbors[0].sort()).toEqual([0, 1]);

      tree.dispose();
    });
  });

  describe('query_ball_point() - Batch queries', () => {
    it('handles batch queries with single radius', () => {
      const data = [[0, 0], [1, 0], [2, 0], [3, 0]];
      const tree = new spatial.KDTree(data);

      const neighbors = tree.query_ball_point([[0.5, 0], [2.5, 0]], 0.6);
      expect(neighbors).toHaveLength(2);
      expect(neighbors[0].sort()).toEqual([0, 1]);
      expect(neighbors[1].sort()).toEqual([2, 3]);

      tree.dispose();
    });

    it('handles batch queries with array of radii', () => {
      const data = [[0, 0], [1, 0], [2, 0], [3, 0]];
      const tree = new spatial.KDTree(data);

      const neighbors = tree.query_ball_point([[0, 0], [2, 0]], [1.5, 0.6]);
      expect(neighbors).toHaveLength(2);
      expect(neighbors[0].sort()).toEqual([0, 1]);  // Radius 1.5 from [0,0]
      expect(neighbors[1].sort()).toEqual([2]);     // Radius 0.6 from [2,0]

      tree.dispose();
    });

    it('throws error for mismatched radius array length', () => {
      const data = [[0, 0], [1, 1]];
      const tree = new spatial.KDTree(data);

      expect(() => tree.query_ball_point([[0, 0], [1, 1]], [1.0, 2.0, 3.0])).toThrow('length');

      tree.dispose();
    });
  });

  describe('Memory management', () => {
    it('allows dispose() to be called', () => {
      const data = [[0, 0], [1, 1]];
      const tree = new spatial.KDTree(data);

      tree.dispose();
      expect(() => tree.dispose()).not.toThrow();  // Double dispose is safe
    });

    it('throws error when querying disposed tree', () => {
      const data = [[0, 0], [1, 1]];
      const tree = new spatial.KDTree(data);
      tree.dispose();

      expect(() => tree.query([0, 0], 1)).toThrow('disposed');
    });
  });

  describe('Edge cases', () => {
    it('handles single point', () => {
      const data = [[42, 17]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([40, 15], 1);
      expect(result.indices).toEqual([0]);

      tree.dispose();
    });

    it('handles duplicate points', () => {
      const data = [[1, 1], [1, 1], [1, 1], [2, 2]];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1, 1], 3);
      // Should find the three [1,1] points with distance 0
      expect(result.distances[0]).toBeCloseTo(0, 10);
      expect(result.distances[1]).toBeCloseTo(0, 10);
      expect(result.distances[2]).toBeCloseTo(0, 10);

      tree.dispose();
    });

    it('handles high-dimensional data', () => {
      const data = [
        [1, 2, 3, 4, 5],
        [2, 3, 4, 5, 6],
        [3, 4, 5, 6, 7]
      ];
      const tree = new spatial.KDTree(data);

      const result = tree.query([1.1, 2.1, 3.1, 4.1, 5.1], 1);
      expect(result.indices[0]).toBe(0);

      tree.dispose();
    });

    it('handles large dataset', () => {
      // Generate 1000 random 2D points
      const data: number[][] = [];
      for (let i = 0; i < 1000; i++) {
        data.push([Math.random() * 100, Math.random() * 100]);
      }

      const tree = new spatial.KDTree(data);

      const result = tree.query([50, 50], 10);
      expect(result.indices).toHaveLength(10);
      expect(result.distances).toHaveLength(10);

      // Verify sorted
      for (let i = 1; i < 10; i++) {
        expect(result.distances[i]).toBeGreaterThanOrEqual(result.distances[i - 1]);
      }

      tree.dispose();
    });
  });
});
