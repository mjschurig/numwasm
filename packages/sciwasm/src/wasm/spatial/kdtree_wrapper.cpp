#include "ckdtree_decl.h"
#include <emscripten.h>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <cmath>

extern "C" {

/**
 * Build a KD-tree from data
 *
 * @param data Flat array of points (n * m doubles, row-major)
 * @param n Number of points
 * @param m Dimensionality
 * @param leafsize Leaf size threshold (stop subdividing when <= leafsize points)
 * @param balanced Use balanced tree (median-based splits)
 * @param compact Recompute bounding boxes at each node
 * @return Pointer to allocated ckdtree struct, or NULL on error
 */
EMSCRIPTEN_KEEPALIVE
ckdtree* kdtree_build(
    double* data,
    int64_t n,
    int64_t m,
    int64_t leafsize,
    int balanced,
    int compact
) {
    // Allocate ckdtree struct
    ckdtree* tree = (ckdtree*)malloc(sizeof(ckdtree));
    if (!tree) return NULL;

    // Initialize tree metadata
    tree->raw_data = data;
    tree->n = n;
    tree->m = m;
    tree->leafsize = leafsize;

    // Allocate indices array (permutation of 0..n-1)
    tree->raw_indices = (ckdtree_intp_t*)malloc(n * sizeof(ckdtree_intp_t));
    if (!tree->raw_indices) {
        free(tree);
        return NULL;
    }
    for (int64_t i = 0; i < n; i++) {
        tree->raw_indices[i] = i;
    }

    // Allocate bounding box arrays
    tree->raw_maxes = (double*)malloc(m * sizeof(double));
    tree->raw_mins = (double*)malloc(m * sizeof(double));
    if (!tree->raw_maxes || !tree->raw_mins) {
        free(tree->raw_indices);
        free(tree->raw_maxes);
        free(tree->raw_mins);
        free(tree);
        return NULL;
    }

    // Compute initial bounding box
    for (int64_t i = 0; i < m; i++) {
        tree->raw_maxes[i] = data[i];
        tree->raw_mins[i] = data[i];
    }
    for (int64_t j = 1; j < n; j++) {
        for (int64_t i = 0; i < m; i++) {
            double val = data[j * m + i];
            if (val > tree->raw_maxes[i]) tree->raw_maxes[i] = val;
            if (val < tree->raw_mins[i]) tree->raw_mins[i] = val;
        }
    }

    // No periodic boundaries
    tree->raw_boxsize_data = NULL;

    // Allocate tree buffer
    tree->tree_buffer = new std::vector<ckdtreenode>();
    tree->tree_buffer->reserve(2 * n - 1);  // Max nodes in balanced tree

    // Build the tree
    int result = build_ckdtree(tree, 0, n, tree->raw_maxes, tree->raw_mins,
                               balanced, compact);

    if (result != 0) {
        // Build failed
        delete tree->tree_buffer;
        free(tree->raw_indices);
        free(tree->raw_maxes);
        free(tree->raw_mins);
        free(tree);
        return NULL;
    }

    // Set ctree to root
    tree->ctree = &((*tree->tree_buffer)[0]);
    tree->size = tree->tree_buffer->size();

    return tree;
}

/**
 * Query k nearest neighbors
 *
 * @param tree Pointer to ckdtree
 * @param x Query points (n_queries * m doubles, row-major)
 * @param n_queries Number of query points
 * @param k Number of neighbors to find
 * @param p Minkowski p-norm (1, 2, or infinity)
 * @param eps Approximation factor (0 = exact)
 * @param distance_upper_bound Maximum search distance (infinity = no limit)
 * @param distances_out Output array for distances (n_queries * k doubles)
 * @param indices_out Output array for indices (n_queries * k int64s)
 * @return 0 on success, non-zero on error
 */
EMSCRIPTEN_KEEPALIVE
int kdtree_query_knn(
    const ckdtree* tree,
    const double* x,
    int64_t n_queries,
    int64_t k,
    double p,
    double eps,
    double distance_upper_bound,
    double* distances_out,
    int64_t* indices_out
) {
    if (!tree || !x || !distances_out || !indices_out) return -1;

    // For each query point, we need to call query_knn
    // The scipy API expects k as an array, so we create one
    std::vector<ckdtree_intp_t> k_array(n_queries, k);

    // Call the C++ query function
    int result = query_knn(
        tree,
        distances_out,           // dd: output distances
        indices_out,             // ii: output indices
        x,                       // xx: query points
        n_queries,               // n: number of queries
        k_array.data(),          // k: array of k values
        n_queries,               // nk: length of k array
        k,                       // kmax: maximum k
        eps,                     // eps: approximation tolerance
        p,                       // p: Minkowski p-norm
        distance_upper_bound     // distance_upper_bound
    );

    return result;
}

/**
 * Query all neighbors within radius r
 *
 * @param tree Pointer to ckdtree
 * @param x Query points (n_queries * m doubles, row-major)
 * @param n_queries Number of query points
 * @param r Radius (single value or n_queries values)
 * @param p Minkowski p-norm (1, 2, or infinity)
 * @param eps Approximation factor (0 = exact)
 * @param return_length If true, return counts instead of indices
 * @param sort_output If true, sort results by distance
 * @param results_out Output array pointer (will be allocated)
 * @param counts_out Output array for result counts per query
 * @return Total number of results, or -1 on error
 */
EMSCRIPTEN_KEEPALIVE
int64_t kdtree_query_ball_point(
    const ckdtree* tree,
    const double* x,
    int64_t n_queries,
    const double* r,
    double p,
    double eps,
    int return_length,
    int sort_output,
    int64_t** results_out,
    int64_t* counts_out
) {
    if (!tree || !x || !r || !results_out || !counts_out) return -1;

    // Results vector - will be filled by query_ball_point
    std::vector<ckdtree_intp_t> results;

    // Call the C++ query function
    int result = query_ball_point(
        tree,
        x,                                    // x: query points
        r,                                    // r: radii
        p,                                    // p: Minkowski p-norm
        eps,                                  // eps: approximation tolerance
        n_queries,                            // n_queries
        &results,                             // results: output vector
        return_length != 0,                   // return_length
        sort_output != 0                      // sort_output
    );

    if (result != 0) return -1;

    // Parse results vector
    // Format: [count1, idx1_1, idx1_2, ..., count2, idx2_1, idx2_2, ...]
    int64_t total_results = 0;
    size_t pos = 0;

    for (int64_t i = 0; i < n_queries; i++) {
        if (pos >= results.size()) {
            counts_out[i] = 0;
        } else {
            counts_out[i] = results[pos];
            total_results += results[pos];
            pos += 1 + results[pos];
        }
    }

    // Allocate flat output array
    *results_out = (int64_t*)malloc(results.size() * sizeof(int64_t));
    if (!*results_out) return -1;

    // Copy results
    for (size_t i = 0; i < results.size(); i++) {
        (*results_out)[i] = results[i];
    }

    return results.size();
}

/**
 * Free a KD-tree and all associated memory
 *
 * @param tree Pointer to ckdtree to free
 */
EMSCRIPTEN_KEEPALIVE
void kdtree_free(ckdtree* tree) {
    if (!tree) return;

    // Free tree buffer
    if (tree->tree_buffer) {
        delete tree->tree_buffer;
        tree->tree_buffer = NULL;
    }

    // Free indices (we allocated this)
    if (tree->raw_indices) {
        free(tree->raw_indices);
        tree->raw_indices = NULL;
    }

    // Free bounding boxes (we allocated these)
    if (tree->raw_maxes) {
        free(tree->raw_maxes);
        tree->raw_maxes = NULL;
    }
    if (tree->raw_mins) {
        free(tree->raw_mins);
        tree->raw_mins = NULL;
    }

    // Note: raw_data is owned by JavaScript, don't free it

    // Free the tree struct itself
    free(tree);
}

}  // extern "C"
