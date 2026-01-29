import { getWasmModule } from '../wasm-loader.js';
import type {
  QueryOptions,
  BallQueryOptions,
  QueryResult,
  QueryResultSingle,
  QueryResultBatch
} from './types.js';

/**
 * Helper: Copy JavaScript array to WASM heap
 */
function toWasmF64(wasm: any, arr: Float64Array | number[]): number {
  const typedArr = arr instanceof Float64Array ? arr : new Float64Array(arr);
  const ptr = wasm._malloc(typedArr.byteLength);
  wasm.HEAPF64.set(typedArr, ptr >> 3);
  return ptr;
}

/**
 * Helper: Read Float64Array from WASM heap
 */
function fromWasmF64(wasm: any, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

/**
 * Helper: Read Int64Array from WASM heap (stored as pairs of 32-bit ints)
 */
function fromWasmI64(wasm: any, ptr: number, len: number): number[] {
  const result: number[] = [];
  const heap32 = wasm.HEAP32;
  const offset = ptr >> 2;

  for (let i = 0; i < len; i++) {
    // For now, assume values fit in 32 bits (tree indices typically do)
    // Read as int32 pairs but only use lower 32 bits
    const low = heap32[offset + i * 2];
    // Skip high word (heap32[offset + i * 2 + 1]) - for indices this is always 0
    result.push(low);
  }
  return result;
}

/**
 * KD-tree for fast spatial queries
 *
 * A k-d tree (short for k-dimensional tree) is a space-partitioning data structure
 * for organizing points in k-dimensional space. It's particularly efficient for
 * nearest neighbor and range searches.
 *
 * @example
 * ```typescript
 * import { spatial } from 'sciwasm';
 *
 * // Create tree from point cloud
 * const points = [[0, 0], [1, 1], [2, 2], [3, 3]];
 * const tree = new spatial.KDTree(points);
 *
 * // Find 3 nearest neighbors
 * const result = tree.query([1.5, 1.5], 3);
 * console.log(result.indices);  // [1, 2, 0]
 * console.log(result.distances); // [0.707..., 0.707..., 2.121...]
 *
 * // Find all neighbors within radius 1.0
 * const neighbors = tree.query_ball_point([1, 1], 1.0);
 * console.log(neighbors); // [[0, 1]]
 *
 * // Clean up
 * tree.dispose();
 * ```
 */
export class KDTree {
  private treePtr: number;
  private dataPtr: number;
  private n: number;
  private m: number;
  private disposed: boolean = false;

  /**
   * Construct a KD-tree from data
   *
   * @param data Array of points, shape (n, m) where n = number of points, m = dimensionality
   * @param leafsize Number of points at which to stop subdividing (default: 10)
   * @param balanced Use balanced tree with median-based splits (default: true)
   * @param compact Recompute tight bounding boxes at each node (default: true, slower build but faster queries)
   */
  constructor(
    data: number[][],
    leafsize: number = 10,
    balanced: boolean = true,
    compact: boolean = true
  ) {
    const wasm = getWasmModule();

    // Validate input
    if (!data || data.length === 0) {
      throw new Error('KDTree: data must be a non-empty array');
    }

    this.n = data.length;
    this.m = data[0].length;

    // Check all rows have same length
    for (let i = 1; i < this.n; i++) {
      if (data[i].length !== this.m) {
        throw new Error(`KDTree: all points must have the same dimensionality (row 0 has ${this.m}, row ${i} has ${data[i].length})`);
      }
    }

    // Flatten data to row-major order
    const flatData = new Float64Array(this.n * this.m);
    for (let i = 0; i < this.n; i++) {
      for (let j = 0; j < this.m; j++) {
        flatData[i * this.m + j] = data[i][j];
      }
    }

    // Allocate and copy data to WASM
    this.dataPtr = toWasmF64(wasm, flatData);

    // Build tree
    this.treePtr = wasm._kdtree_build(
      this.dataPtr,
      this.n,
      this.m,
      leafsize,
      balanced ? 1 : 0,
      compact ? 1 : 0
    );

    if (this.treePtr === 0) {
      wasm._free(this.dataPtr);
      throw new Error('KDTree: failed to build tree (allocation error)');
    }
  }

  /**
   * Query the tree for the k nearest neighbors
   *
   * @param x Query point (1D array of length m) or array of query points (2D array of shape n_queries x m)
   * @param k Number of nearest neighbors to return (default: 1)
   * @param options Query options (eps, p, distance_upper_bound)
   * @returns Object with distances and indices arrays
   */
  query(
    x: number[] | number[][],
    k: number = 1,
    options?: QueryOptions
  ): QueryResult {
    if (this.disposed) {
      throw new Error('KDTree: cannot query disposed tree');
    }

    const wasm = getWasmModule();
    const eps = options?.eps ?? 0;
    const p = options?.p ?? 2;
    const distance_upper_bound = options?.distance_upper_bound ?? Infinity;

    // Determine if single or batch query
    const isBatch = Array.isArray(x[0]);
    let n_queries: number;
    let queryData: Float64Array;

    if (isBatch) {
      const queries = x as number[][];
      n_queries = queries.length;

      // Validate query dimensionality
      for (let i = 0; i < n_queries; i++) {
        if (queries[i].length !== this.m) {
          throw new Error(`KDTree.query: query point ${i} has wrong dimensionality (expected ${this.m}, got ${queries[i].length})`);
        }
      }

      // Flatten queries
      queryData = new Float64Array(n_queries * this.m);
      for (let i = 0; i < n_queries; i++) {
        for (let j = 0; j < this.m; j++) {
          queryData[i * this.m + j] = queries[i][j];
        }
      }
    } else {
      const query = x as number[];
      n_queries = 1;

      if (query.length !== this.m) {
        throw new Error(`KDTree.query: query point has wrong dimensionality (expected ${this.m}, got ${query.length})`);
      }

      queryData = new Float64Array(query);
    }

    // Allocate memory for query and results
    const xPtr = toWasmF64(wasm, queryData);
    const distPtr = wasm._malloc(n_queries * k * 8);  // 8 bytes per double
    const idxPtr = wasm._malloc(n_queries * k * 8);   // 8 bytes per int64

    try {
      // Call WASM function
      const result = wasm._kdtree_query_knn(
        this.treePtr,
        xPtr,
        n_queries,
        k,
        p,
        eps,
        distance_upper_bound,
        distPtr,
        idxPtr
      );

      if (result !== 0) {
        throw new Error('KDTree.query: query failed');
      }

      // Read results
      const distances = fromWasmF64(wasm, distPtr, n_queries * k);
      const indices = fromWasmI64(wasm, idxPtr, n_queries * k);

      if (isBatch) {
        // Reshape to 2D arrays
        const distArr: number[][] = [];
        const idxArr: number[][] = [];
        for (let i = 0; i < n_queries; i++) {
          distArr.push(Array.from(distances.slice(i * k, (i + 1) * k)));
          idxArr.push(indices.slice(i * k, (i + 1) * k));
        }
        return { distances: distArr, indices: idxArr } as QueryResultBatch;
      } else {
        // Return 1D arrays
        return {
          distances: Array.from(distances),
          indices: indices
        } as QueryResultSingle;
      }
    } finally {
      // Clean up
      wasm._free(xPtr);
      wasm._free(distPtr);
      wasm._free(idxPtr);
    }
  }

  /**
   * Query the tree for all neighbors within distance r
   *
   * @param x Query point (1D array) or array of query points (2D array)
   * @param r Radius (single number) or array of radii (one per query point)
   * @param options Query options (eps, p, return_sorted, return_length)
   * @returns Array of neighbor index arrays (one array per query point), or array of counts if return_length=true
   */
  query_ball_point(
    x: number[] | number[][],
    r: number | number[],
    options?: BallQueryOptions
  ): number[][] | number[] {
    if (this.disposed) {
      throw new Error('KDTree: cannot query disposed tree');
    }

    const wasm = getWasmModule();
    const eps = options?.eps ?? 0;
    const p = options?.p ?? 2;
    const return_sorted = options?.return_sorted ?? false;
    const return_length = options?.return_length ?? false;

    // Determine if single or batch query
    const isBatch = Array.isArray(x[0]);
    let n_queries: number;
    let queryData: Float64Array;
    let radiusData: Float64Array;

    if (isBatch) {
      const queries = x as number[][];
      n_queries = queries.length;

      // Validate query dimensionality
      for (let i = 0; i < n_queries; i++) {
        if (queries[i].length !== this.m) {
          throw new Error(`KDTree.query_ball_point: query point ${i} has wrong dimensionality (expected ${this.m}, got ${queries[i].length})`);
        }
      }

      // Flatten queries
      queryData = new Float64Array(n_queries * this.m);
      for (let i = 0; i < n_queries; i++) {
        for (let j = 0; j < this.m; j++) {
          queryData[i * this.m + j] = queries[i][j];
        }
      }

      // Handle radius
      if (Array.isArray(r)) {
        if (r.length !== n_queries) {
          throw new Error(`KDTree.query_ball_point: radius array length (${r.length}) must match number of queries (${n_queries})`);
        }
        radiusData = new Float64Array(r);
      } else {
        radiusData = new Float64Array(n_queries).fill(r);
      }
    } else {
      const query = x as number[];
      n_queries = 1;

      if (query.length !== this.m) {
        throw new Error(`KDTree.query_ball_point: query point has wrong dimensionality (expected ${this.m}, got ${query.length})`);
      }

      queryData = new Float64Array(query);
      radiusData = new Float64Array([Array.isArray(r) ? r[0] : r]);
    }

    // Allocate memory
    const xPtr = toWasmF64(wasm, queryData);
    const rPtr = toWasmF64(wasm, radiusData);
    const countsPtr = wasm._malloc(n_queries * 8);  // 8 bytes per int64
    const resultsPtrPtr = wasm._malloc(8);          // Pointer to results array

    try {
      // Call WASM function
      const totalResults = wasm._kdtree_query_ball_point(
        this.treePtr,
        xPtr,
        n_queries,
        rPtr,
        p,
        eps,
        return_length ? 1 : 0,
        return_sorted ? 1 : 0,
        resultsPtrPtr,
        countsPtr
      );

      if (totalResults < 0) {
        throw new Error('KDTree.query_ball_point: query failed');
      }

      // Read counts
      const counts = fromWasmI64(wasm, countsPtr, n_queries);

      // Read pointer to results array
      const resultsPtr = wasm.HEAP32[resultsPtrPtr >> 2];

      // Read results (flat array with format: [count1, idx1_1, ..., count2, idx2_1, ...])
      const flatResults = fromWasmI64(wasm, resultsPtr, totalResults);

      // Parse flat results into arrays
      const results: number[][] = [];
      let pos = 0;

      for (let i = 0; i < n_queries; i++) {
        const count = flatResults[pos];
        pos++;

        const neighbors: number[] = [];
        for (let j = 0; j < count; j++) {
          neighbors.push(flatResults[pos]);
          pos++;
        }

        results.push(neighbors);
      }

      // Free the results array allocated in C++
      wasm._free(resultsPtr);

      if (return_length) {
        return counts;
      }

      return isBatch ? results : results[0] ? [results[0]] : [[]];
    } finally {
      // Clean up
      wasm._free(xPtr);
      wasm._free(rPtr);
      wasm._free(countsPtr);
      wasm._free(resultsPtrPtr);
    }
  }

  /**
   * Free all memory associated with this tree
   *
   * After calling dispose(), this tree cannot be used for queries.
   */
  dispose(): void {
    if (this.disposed) return;

    const wasm = getWasmModule();

    // Free tree (this also frees internal allocations)
    wasm._kdtree_free(this.treePtr);

    // Free data array
    wasm._free(this.dataPtr);

    this.disposed = true;
  }

  /**
   * Get the number of points in the tree
   */
  get size(): number {
    return this.n;
  }

  /**
   * Get the dimensionality of points in the tree
   */
  get ndim(): number {
    return this.m;
  }
}
