/*****************************************************************************
 * netboost - in-memory MC-UPGMA clustering core
 *
 * Portable C++/Rcpp reimplementation of the single-round, memory-resident
 * core of the MC-UPGMA clustering package by Yaniv Loewenstein:
 *
 *   Loewenstein Y, Portugaly E, Fromer M, Linial M. Efficient algorithms for
 *   accurate hierarchical clustering of huge datasets: tackling the entire
 *   protein space. Bioinformatics. 2008 Jul 1;24(13):i41-9.
 *
 * netboost previously shelled out to Loewenstein's bundled C++/Perl/GNU-make
 * suite (src/mcupgma/) for the sparse-distance clustering step. That suite is
 * unix-only (Perl driver + Makefiles + external binaries). This file provides
 * a from-scratch reimplementation of the exact clustering performed by
 * netboost's invocation
 *
 *   cluster.pl -max_distance 1 -max_singleton N -heap_size 1e7 ...
 *
 * which, because the per-round heap (1e7) far exceeds the number of edges and
 * all distances lie in [0, max_distance], reduces to a single in-memory round
 * of average-linkage (UPGMA) agglomeration on the sparse edge graph, treating
 * every missing pair as distance max_distance (psi). The result is a forest
 * (one tree per connected component of the edge graph); isolated features stay
 * singletons.
 *
 * It reproduces, bit-for-bit, the merge sequence, cluster-id numbering and the
 * 6-significant-figure distance values that the original pipeline wrote to its
 * tree file and that netboost read back with read.table().
 *
 * This reimplementation is distributed under the same terms as the rest of
 * netboost (GPL-3), and is a derivative work of the GPL-licensed MC-UPGMA.
 *****************************************************************************/

// [[Rcpp::interfaces(r)]]

#include <Rcpp.h>

#include <vector>
#include <unordered_map>
#include <set>
#include <sstream>
#include <locale>
#include <cstdlib>
#include <cmath>

using namespace Rcpp;

namespace {

typedef unsigned int tClusterId;

/*
 * An active edge between two currently-alive clusters (low < high), holding the
 * full-precision (size-weighted average-linkage) distance. Edges are ordered for
 * selection exactly as the original tournament heap orders them: by minimal
 * distance, ties broken by the LARGEST (low, high) pair lexicographically
 * (this is what HierarchicalClustering's ubEdgeHeap.top() returns).
 */
struct tEdge {
    double      dist;
    tClusterId  low;
    tClusterId  high;
};

struct fEdgeOrder {
    bool operator()(tEdge const & a, tEdge const & b) const {
        if (a.dist < b.dist) return true;
        if (b.dist < a.dist) return false;
        if (a.low != b.low)  return a.low > b.low;
        return a.high > b.high;
    }
};

typedef std::unordered_map<tClusterId, double> tAdj;

/*
 * Round a full-precision distance through the exact lossy path the old pipeline
 * used: default C++ ostream formatting (precision 6 -> 6 significant figures)
 * to text, then strtod back to a double. imbue(classic()) forces a '.' decimal
 * separator regardless of LC_NUMERIC (a benign hardening; the original only
 * ever produced correct output under a C/POSIX numeric locale anyway).
 */
inline double round6(double x) {
    std::ostringstream oss;
    oss.imbue(std::locale::classic());
    oss << x;
    return std::strtod(oss.str().c_str(), nullptr);
}

} // anonymous namespace

//' @title In-memory MC-UPGMA clustering (single round)
//'
//' @description Portable reimplementation of the sparse average-linkage
//'   (UPGMA) clustering performed by netboost. Missing pairs are treated as
//'   distance \code{max_distance}; merging proceeds in order of increasing
//'   distance (ties broken by the largest cluster-id pair) until no edge with
//'   distance <= \code{max_distance} remains, yielding a forest.
//'
//' @param low Integer vector, smaller cluster id of each input edge (1-based).
//' @param high Integer vector, larger cluster id of each input edge (1-based).
//' @param dist Numeric vector of edge distances (in [0, max_distance]).
//' @param max_singleton Numeric. Maximum singleton id; new (merged) cluster ids
//'   start at \code{max_singleton + 1}.
//' @param max_distance Numeric. Upper distance bound (psi); also the distance
//'   assigned to missing pairs.
//' @return Numeric matrix with columns cluster_id1, cluster_id2, distance,
//'   cluster_id3 (one row per merge, in merge order).
// [[Rcpp::export(name = "cpp_mcupgma")]]
NumericMatrix cpp_mcupgma(IntegerVector low,
                          IntegerVector high,
                          NumericVector dist,
                          double max_singleton,
                          double max_distance) {
    const R_xlen_t ne = dist.size();
    if (low.size() != ne || high.size() != ne)
        Rcpp::stop("cpp_mcupgma: low, high and dist must have equal length.");

    const tClusterId n  = static_cast<tClusterId>(max_singleton);
    const double psi    = max_distance;

    // Per-cluster size (sizes start at 1 for singletons). Indexed by cluster id;
    // capacity covers all merge ids (<= 2*n - 1).
    std::vector<tClusterId> csize(static_cast<size_t>(n) * 2 + 2, 0);
    for (tClusterId i = 1; i <= n; ++i) csize[i] = 1;

    // Adjacency: cluster id -> (neighbour id -> full-precision distance).
    std::unordered_map<tClusterId, tAdj> adj;
    // Selection structure: all active edges, ordered as the original heap.
    std::set<tEdge, fEdgeOrder> active;

    // Load input edges (in input row order). Self/duplicate edges are the
    // caller's responsibility (nb_filter guarantees uniqueness), mirroring the
    // original tool's precondition.
    for (R_xlen_t e = 0; e < ne; ++e) {
        tClusterId a = static_cast<tClusterId>(low[e]);
        tClusterId b = static_cast<tClusterId>(high[e]);
        double     d = dist[e];
        if (a == b) continue;
        if (a > b) std::swap(a, b);
        adj[a][b] = d;
        adj[b][a] = d;
        tEdge ed; ed.dist = d; ed.low = a; ed.high = b;
        active.insert(ed);
    }

    // Output merge events (low, high, rounded distance, merged id).
    std::vector<double> out_low, out_high, out_dist, out_merged;
    out_low.reserve(static_cast<size_t>(n));
    out_high.reserve(static_cast<size_t>(n));
    out_dist.reserve(static_cast<size_t>(n));
    out_merged.reserve(static_cast<size_t>(n));

    tClusterId nextId = n + 1;

    while (!active.empty()) {
        const tEdge top = *active.begin();
        if (top.dist > psi) break;          // no mergeable edge remains

        const tClusterId lo = top.low;
        const tClusterId hi = top.high;
        const double     md = top.dist;
        const tClusterId k  = nextId++;
        const double nL = static_cast<double>(csize[lo]);
        const double nH = static_cast<double>(csize[hi]);
        const tClusterId denom = csize[lo] + csize[hi];

        // Snapshot the merged clusters' neighbour maps so we can freely mutate
        // adj (which may rehash) without invalidating references.
        tAdj Lmap; Lmap.swap(adj[lo]);
        tAdj Hmap; Hmap.swap(adj[hi]);
        tAdj Nmap;                          // neighbours of the new cluster k
        Nmap.reserve(Lmap.size() + Hmap.size());

        // Neighbours of the low cluster (combined with high's where shared).
        for (tAdj::const_iterator it = Lmap.begin(); it != Lmap.end(); ++it) {
            const tClusterId j = it->first;
            if (j == hi) continue;
            const double dlj = it->second;          // d(low, j)  (first operand)
            tAdj::const_iterator hj = Hmap.find(j);
            const double dhj = (hj != Hmap.end()) ? hj->second : psi; // d(high,j)/psi
            const double nd  = (dlj * nL + dhj * nH) / static_cast<double>(denom);

            Nmap[j] = nd;
            // Rewire j: drop low (and high if present), point to k.
            tAdj & Jmap = adj[j];
            active.erase(tEdge{dlj, std::min(lo, j), std::max(lo, j)});
            Jmap.erase(lo);
            if (hj != Hmap.end()) {
                active.erase(tEdge{hj->second, std::min(hi, j), std::max(hi, j)});
                Jmap.erase(hi);
            }
            Jmap[k] = nd;
            active.insert(tEdge{nd, std::min(j, k), std::max(j, k)});
        }
        // Neighbours of the high cluster not shared with low (low side missing).
        for (tAdj::const_iterator it = Hmap.begin(); it != Hmap.end(); ++it) {
            const tClusterId j = it->first;
            if (j == lo) continue;
            if (Lmap.find(j) != Lmap.end()) continue;   // already handled above
            const double dhj = it->second;              // d(high, j) (second operand)
            const double nd  = (psi * nL + dhj * nH) / static_cast<double>(denom);

            Nmap[j] = nd;
            tAdj & Jmap = adj[j];
            active.erase(tEdge{dhj, std::min(hi, j), std::max(hi, j)});
            Jmap.erase(hi);
            Jmap[k] = nd;
            active.insert(tEdge{nd, std::min(j, k), std::max(j, k)});
        }

        // Remove the merging edge itself and retire the merged clusters.
        active.erase(tEdge{md, lo, hi});
        adj.erase(lo);
        adj.erase(hi);
        csize[lo] = 0;
        csize[hi] = 0;
        csize[k]  = denom;
        adj[k].swap(Nmap);

        out_low.push_back(static_cast<double>(lo));
        out_high.push_back(static_cast<double>(hi));
        out_dist.push_back(round6(md));
        out_merged.push_back(static_cast<double>(k));
    }

    const R_xlen_t m = static_cast<R_xlen_t>(out_low.size());
    NumericMatrix forest(m, 4);
    for (R_xlen_t r = 0; r < m; ++r) {
        forest(r, 0) = out_low[r];
        forest(r, 1) = out_high[r];
        forest(r, 2) = out_dist[r];
        forest(r, 3) = out_merged[r];
    }
    colnames(forest) = CharacterVector::create(
        "cluster_id1", "cluster_id2", "distance", "cluster_id3");
    return forest;
}
