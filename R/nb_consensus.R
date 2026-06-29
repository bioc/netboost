###############################################################################
# Multi-dataset consensus network analysis for the netboost ecosystem
#
# 10-step workflow:
#   1+2. Filtered adjacency per dataset (netboost::nb_filter)
#   3+4. Direction filter + consensus at the adjacency level
#   5.   Soft thresholding
#   6.   Consensus TOM (nb_dist per dataset, then pmax/pmin)
#   7.   Clustering (hclust + cutreeDynamic)
#   8.   Module eigengenes (aligned to mean expression per dataset)
#   9.   Iterative hclust-based cross-dataset module merging
#   10.  Module retention
#
# Output: module assignments, per-dataset eigengenes, PC weights,
#         per-dataset filters, rotation, var_explained, PDF report.
#
# Compatible with: nb_transfer, nb_plot_dendro (via nb_summary fields)
###############################################################################

# NOTE (packaging): no top-level library()/allowWGCNAThreads()/options() calls.
# Package code must not mutate global state at load time. WGCNA threading is
# enabled inside nb_consensus() only when cores > 1. Dependencies are declared
# via the @importFrom tags below and the package DESCRIPTION Imports field.

# Sign helper: TRUE iff cor(x, y) is a real, negative number. Guards the
# eigengene/PC sign-alignment steps against zero-variance inputs, where
# stats::cor() returns NA (and `if (NA < 0)` would error).
#' @noRd
.neg_cor <- function(x, y, use = "complete.obs") {
    r <- suppressWarnings(stats::cor(x, y, use = use))
    isTRUE(r < 0)
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Multi-dataset consensus network analysis
#'
#' Builds a consensus gene/feature network across two or more datasets that
#' share the same features. Combines a per-dataset boosting/correlation filter
#' (\code{\link{nb_filter}}), a consensus topological overlap matrix computed
#' from per-dataset \code{\link{nb_dist}} distances, dynamic-tree clustering,
#' iterative cross-dataset module merging, and module retention based on
#' cross-dataset PC1 concordance.
#'
#' @param datan_list List of at least two data frames or matrices
#'   (rows = samples, columns = features). All datasets must have identical
#'   column names in the same order.
#' @param filter_method Edge pre-filter passed to \code{\link{nb_filter}};
#'   one of \code{"spearman"}, \code{"pearson"}, \code{"kendall"},
#'   \code{"boosting"}, or \code{"skip"} (no filtering).
#' @param soft_power Optional non-negative soft-thresholding power applied to
#'   the consensus adjacency. \code{NULL} disables soft thresholding.
#' @param consensus_method How per-dataset edge weights are aggregated:
#'   \code{"min"} (edge must exist in all datasets), \code{"max"} (in at least
#'   two), or \code{"median"} (in at least half).
#' @param network_type \code{"unsigned"} (default) or \code{"signed"}. When
#'   \code{"signed"}, only positively co-directional gene pairs are kept;
#'   \code{"unsigned"} keeps both directions (using \code{|cor|}). Independent
#'   of \code{filter_method} -- both network types work with every filter
#'   (boosting and the correlation-based filters alike).
#' @param filter_dir Logical; for \code{network_type = "unsigned"}, require
#'   the correlation sign to agree across datasets.
#' @param min_cluster_size Minimum module size for \code{cutreeDynamic} and for
#'   final module retention.
#' @param ME_diss_thres Module-eigengene dissimilarity threshold (1 - cor) for
#'   merging modules; must be in (0, 1).
#' @param merge_criterion \code{"all"} (use the minimum cross-dataset ME
#'   correlation) or \code{"any"} (use the maximum) when merging.
#' @param module_retention \code{"all"}, \code{"any"}, or \code{"median"};
#'   how cross-dataset PC1 correlations are summarised for retention.
#' @param min_pc1_cor Minimum cross-dataset PC1 loading correlation for a
#'   module to be retained.
#' @param deepSplit Sensitivity (0-4) passed to \code{cutreeDynamic}.
#' @param method Correlation method for adjacency and \code{\link{nb_dist}};
#'   one of \code{"pearson"}, \code{"kendall"}, \code{"spearman"}.
#' @param scale Logical; scale features before PCA/eigengene computation.
#' @param stepno Boosting steps passed to \code{\link{nb_filter}}.
#' @param n_pc Number of principal components for module eigengenes.
#' @param robust_PCs Logical; use rank-based (robust) PCA. When \code{TRUE} the
#'   stored \code{rotation} is computed on rank-transformed data so it is
#'   consistent with the rank-based projection in \code{nb_consensus_transfer}.
#' @param nb_min_varExpl Minimum variance explained for an eigengene to be
#'   considered valid (passed to \code{\link{nb_moduleEigengenes}}).
#' @param cores Number of CPU cores.
#' @param pdf_report Optional path; if set, a multi-page PDF report is written.
#' @param mask_cache_dir Optional directory for caching filter masks across
#'   runs that differ only in soft power.
#' @param verbose Logical; emit progress messages.
#'
#' @return Invisibly, a list compatible with \code{\link{nb_summary}} fields
#'   (\code{names}, \code{colors}, \code{MEs}, \code{var_explained},
#'   \code{rotation}, \code{filter}) plus consensus-specific fields
#'   (\code{colors_initial}, \code{colors_merged}, \code{modules_retained},
#'   \code{ME_all_datasets}, \code{pc1_weights}, \code{filter_stats},
#'   \code{params}, and module counts).
#'
#' @seealso \code{\link{nb_consensus_transfer}}, \code{\link{nb_filter}},
#'   \code{\link{nb_dist}}, \code{\link{nb_moduleEigengenes}}
#'
#' @examples
#' set.seed(1)
#' mk <- function(n) {
#'     m <- matrix(rnorm(n * 12), n, 12)
#'     m[, 1:4] <- m[, 1:4] + rnorm(n)      # correlated block 1
#'     m[, 5:8] <- m[, 5:8] + rnorm(n)      # correlated block 2
#'     colnames(m) <- paste0("g", 1:12)
#'     m
#' }
#' ds <- list(mk(40), mk(45))
#' res <- nb_consensus(ds, filter_method = "pearson", min_cluster_size = 3,
#'                     stepno = 5L, verbose = FALSE)
#' table(res$colors)
#'
#' @importFrom stats hclust cutree as.dist prcomp cor lm setNames
#' @importFrom grDevices pdf dev.off colorRampPalette adjustcolor
#' @importFrom graphics barplot text pie plot.new title par image axis abline mtext plot
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
nb_consensus <- function(
    datan_list,
    filter_method    = c("spearman", "pearson", "kendall", "boosting", "skip"),
    soft_power       = NULL,
    consensus_method = c("min", "max", "median"),
    network_type     = "unsigned",
    filter_dir       = TRUE,
    min_cluster_size = 15L,
    ME_diss_thres    = 0.25,
    merge_criterion  = c("all", "any"),
    module_retention = c("all", "any", "median"),
    min_pc1_cor      = 0.8,
    deepSplit        = 2,
    method           = c("pearson", "kendall", "spearman"),
    scale            = TRUE,
    stepno           = 20L,
    n_pc             = 1,
    robust_PCs       = FALSE,
    nb_min_varExpl   = 0.5,
    cores            = 1L,
    pdf_report       = NULL,
    mask_cache_dir   = NULL,
    verbose          = TRUE
) {
    filter_method    <- match.arg(filter_method)
    consensus_method <- match.arg(consensus_method)
    # Honor caller-provided network_type ("unsigned" or "signed").
    network_type <- match.arg(network_type, c("unsigned", "signed"))
    method       <- match.arg(method)
    merge_criterion  <- match.arg(merge_criterion)
    module_retention <- match.arg(module_retention)

    # Note: direction filter checks correlation sign agreement across datasets


    # ==================================================================
    # Input validation
    # ==================================================================
    if (!is.list(datan_list))
        stop("datan_list must be a list of data frames or matrices.")
    if (length(datan_list) < 2)
        stop("datan_list must contain at least 2 datasets.")

    n_ds <- length(datan_list)

    for (i in seq_len(n_ds)) {
        d <- datan_list[[i]]
        if (!is.data.frame(d) && !is.matrix(d))
            stop(sprintf("Dataset %d is not a data.frame or matrix.", i))
        if (nrow(d) < 3)
            stop(sprintf("Dataset %d has only %d samples (need >= 3).", i, nrow(d)))
        if (ncol(d) < 2)
            stop(sprintf("Dataset %d has only %d features (need >= 2).", i, ncol(d)))
        if (is.null(colnames(d)))
            stop(sprintf("Dataset %d has no column names.", i))
        # Numeric check (data frames AND matrices -- a character matrix would
        # otherwise slip through and fail later inside cor() with a cryptic error)
        if (is.data.frame(d)) {
            non_num <- !vapply(d, is.numeric, logical(1))
            if (any(non_num))
                stop(sprintf("Dataset %d has non-numeric columns: %s",
                             i, paste(names(non_num)[non_num], collapse = ", ")))
        } else if (!is.numeric(d)) {
            stop(sprintf("Dataset %d is a non-numeric matrix (storage mode '%s').",
                         i, storage.mode(d)))
        }
        # Duplicate feature names corrupt the gene -> module mapping
        dup <- unique(colnames(d)[duplicated(colnames(d))])
        if (length(dup))
            stop(sprintf("Dataset %d has duplicate feature names: %s",
                         i, paste(dup[seq_len(min(5L, length(dup)))], collapse = ", ")))
        # NA/NaN/Inf are not supported by nb_filter/nb_dist; fail fast and clearly
        # rather than deep inside the C++ filter with an opaque message.
        dm <- as.matrix(d)
        if (anyNA(dm))
            stop(sprintf("Dataset %d contains NA/NaN values; nb_consensus does not support missing data.", i))
        if (!all(is.finite(dm)))
            stop(sprintf("Dataset %d contains non-finite values (Inf); remove or impute them first.", i))
        # Zero-variance features break correlation and eigengene computation
        sds <- apply(dm, 2L, stats::sd)
        zv  <- which(sds == 0)
        if (length(zv))
            stop(sprintf("Dataset %d has %d zero-variance feature(s): %s",
                         i, length(zv),
                         paste(colnames(dm)[zv][seq_len(min(5L, length(zv)))],
                               collapse = ", ")))
    }

    p      <- ncol(datan_list[[1]])
    fnames <- colnames(datan_list[[1]])
    for (i in 2:n_ds) {
        if (ncol(datan_list[[i]]) != p)
            stop(sprintf("Dataset %d has %d features but dataset 1 has %d.",
                         i, ncol(datan_list[[i]]), p))
        if (!all(colnames(datan_list[[i]]) == fnames))
            stop(sprintf("Column names in dataset %d do not match dataset 1.", i))
    }

    if (!is.null(soft_power)) {
        if (!is.numeric(soft_power) || soft_power < 0)
            stop("soft_power must be a non-negative number or NULL.")
    }
    if (!is.numeric(ME_diss_thres) || ME_diss_thres <= 0 || ME_diss_thres >= 1)
        stop("ME_diss_thres must be between 0 and 1 (exclusive).")
    if (!is.numeric(min_cluster_size) || min_cluster_size < 2)
        stop("min_cluster_size must be >= 2.")
    # A correlation can never exceed [-1, 1]; an out-of-range threshold (e.g. a
    # typo of 8 for 0.8) would silently retain nothing / everything.
    if (!is.numeric(min_pc1_cor) || length(min_pc1_cor) != 1 ||
        min_pc1_cor < -1 || min_pc1_cor > 1)
        stop("min_pc1_cor must be a single number in [-1, 1].")
    if (!is.numeric(n_pc) || length(n_pc) != 1 || n_pc < 1 ||
        n_pc != round(n_pc))
        stop("n_pc must be a positive integer.")
    if (!is.numeric(nb_min_varExpl) || length(nb_min_varExpl) != 1 ||
        nb_min_varExpl < 0 || nb_min_varExpl > 1)
        stop("nb_min_varExpl must be between 0 and 1.")

    cores <- as.integer(cores)
    if (is.na(cores) || cores < 1)
        stop("cores must be a positive integer.")
    # nb_filter requires an integer stepno; coerce here so callers may pass a
    # plain numeric (e.g. stepno = 20) without tripping nb_filter's type check.
    stepno <- as.integer(stepno)
    if (is.na(stepno) || stepno < 1L)
        stop("stepno must be a positive integer.")
    if (cores > 1) WGCNA::allowWGCNAThreads(nThreads = cores)
    msg <- if (verbose) message else function(...) invisible(NULL)

    msg("=== Netboost Consensus ===")
    msg(sprintf("  %d datasets, %d features, %d cores", n_ds, p, cores))
    for (i in seq_len(n_ds))
        msg(sprintf("  Dataset %d: %d samples", i, nrow(datan_list[[i]])))

    # Tracking for report
    fstats <- list(n_total = p * (p - 1) / 2,
                   n_after_filter = numeric(n_ds),
                   n_dir_filtered = 0,
                   n_consensus = 0)
    BS <- min(40000L, p)   # block size for memory-safe operations

    # ==================================================================
    # STEP 1+2: Filtered adjacency per dataset
    # - Fast BLAS path when no NAs
    # - Progressive active-gene filtering for min consensus
    # - Optional filter mask caching (skip recomputation for new power)
    # ==================================================================
    # Input validation already rejected NA/Inf, so the complete-data BLAS path
    # is always valid here (nb_filter/nb_dist do not support missing data).
    cor_use <- "everything"

    # --- Filter mask cache ---
    # Build a unique key from dataset dims, gene count, and filter settings
    # Fingerprint each dataset by shape AND content (total + first-row sums) so a
    # cached mask cannot be reused across genuinely different datasets that happen
    # to share dimensions, filter_method and stepno -- which would otherwise yield
    # silently wrong filter masks. Sums are finite (NA/Inf rejected in validation).
    # The key is sanitised to remain a valid filename on every platform (Windows
    # disallows ':' and friends that the '%a' hex-float format can introduce).
    ds_sig <- paste(vapply(datan_list, function(d) {
        m <- as.matrix(d)
        sprintf("%dx%d_%a_%a", nrow(m), ncol(m), sum(m), sum(m[1L, ]))
    }, character(1)), collapse = "_")
    cache_key <- gsub("[^A-Za-z0-9_]", "_",
                      sprintf("mask_%s_p%d_%s_step%d",
                              ds_sig, p, filter_method, stepno))
    cache_dir <- if (!is.null(mask_cache_dir)) mask_cache_dir else tempdir()
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    cached <- NULL
    if (file.exists(cache_file)) {
        msg(sprintf("  Loading cached filter masks: %s", cache_file))
        cached <- readRDS(cache_file)
    }

    msg("\nStep 1-2: Computing filtered adjacencies...")
    adj_raw <- list()
    filter_masks <- list()       # store for caching
    filter_per_ds <- list()      # store sparse filters for return value
    active_genes <- seq_len(p)   # progressive filter: active gene indices

    for (i in seq_len(n_ds)) {
        n_i <- nrow(datan_list[[i]])

        # Progressive filtering: for min consensus, skip inactive genes
        use_subset <- (consensus_method == "min" && i > 1 &&
                       length(active_genes) < p)
        if (use_subset) {
            p_i <- length(active_genes)
            dat_i <- datan_list[[i]][, active_genes, drop = FALSE]
            msg(sprintf("  DS%d: progressive -- %d/%d active genes (%.0f%% saved)",
                        i, p_i, p, 100 * (1 - (p_i / p)^2)))
        } else {
            p_i <- p
            dat_i <- datan_list[[i]]
        }

        if (filter_method != "skip") {
            # --- nb_filter for all non-skip methods ---
            if (!is.null(cached) && length(cached) >= i) {
                msg(sprintf("  DS%d: using cached %s mask", i, filter_method))
                mask_full <- cached[[i]]
                ut_idx <- which(mask_full & upper.tri(mask_full),
                                arr.ind = TRUE)
                filter_per_ds[[i]] <- unname(ut_idx)
                if (use_subset) {
                    mask <- mask_full[active_genes, active_genes]
                } else {
                    mask <- mask_full
                }
                rm(mask_full)
            } else {
                msg(sprintf("  DS%d: nb_filter (method=%s, stepno=%d)...",
                            i, filter_method, stepno))
                nb_edges <- netboost::nb_filter(
                    as.data.frame(dat_i),
                    filter_method = filter_method,
                    stepno = stepno,
                    until  = 0L,
                    mode   = 2L,
                    cores  = cores,
                    verbose = verbose
                )
                # Store sparse filter. Under progressive ("min") subsetting,
                # nb_edges indexes the active-gene subset; map back to full
                # feature indices so the returned filter is consistent with the
                # non-progressive and cached paths (which use full indices).
                if (use_subset && nrow(nb_edges) > 0) {
                    nb_edges_full <- nb_edges
                    nb_edges_full[, 1] <- active_genes[nb_edges[, 1]]
                    nb_edges_full[, 2] <- active_genes[nb_edges[, 2]]
                    filter_per_ds[[i]] <- nb_edges_full
                } else {
                    filter_per_ds[[i]] <- nb_edges
                }
                mask <- matrix(FALSE, p_i, p_i)
                if (nrow(nb_edges) > 0) {
                    mask[cbind(nb_edges[, 1], nb_edges[, 2])] <- TRUE
                    mask[cbind(nb_edges[, 2], nb_edges[, 1])] <- TRUE
                }
                diag(mask) <- FALSE
                rm(nb_edges); gc()
                if (use_subset) {
                    mask_full <- matrix(FALSE, p, p)
                    mask_full[active_genes, active_genes] <- mask
                    filter_masks[[i]] <- mask_full
                    rm(mask_full)
                } else {
                    filter_masks[[i]] <- mask
                }
            }

            msg(sprintf("  DS%d: %s adjacency...", i, method))
            r_a <- WGCNA::cor(dat_i, method = method, use = cor_use)
            adj_sub <- r_a * mask
            diag(adj_sub) <- 0
            rm(r_a, mask); gc()

        } else {
            # --- skip: no filter, just correlation adjacency ---
            msg(sprintf("  DS%d: %s (no filter)...", i, method))
            adj_sub <- WGCNA::cor(dat_i, method = method, use = cor_use)
            diag(adj_sub) <- 0
            filter_per_ds[[i]] <- matrix(integer(0), ncol = 2)
        }

        # Keep signed adjacencies here -- abs() is applied after direction
        # filtering in Step 3-4 (consensus already takes abs() internally).

        # Expand subset back to full p x p if needed
        if (use_subset) {
            adj_raw[[i]] <- matrix(0, p, p)
            adj_raw[[i]][active_genes, active_genes] <- adj_sub
            rm(adj_sub); gc()
        } else {
            adj_raw[[i]] <- adj_sub
        }

        rownames(adj_raw[[i]]) <- fnames
        colnames(adj_raw[[i]]) <- fnames
        fstats$n_after_filter[i] <- sum(adj_raw[[i]][upper.tri(adj_raw[[i]])] != 0)
        msg(sprintf("  DS%d: %d edges after filter", i, fstats$n_after_filter[i]))

        # Update active genes for progressive filtering
        if (consensus_method == "min") {
            has_edge <- logical(p)
            for (s in seq(1, p, BS)) {
                e <- min(s + BS - 1, p)
                has_edge[s:e] <- rowSums(adj_raw[[i]][s:e, , drop = FALSE] != 0) > 0
            }
            active_genes <- intersect(active_genes, which(has_edge))
            msg(sprintf("  Active genes after DS%d: %d/%d", i, length(active_genes), p))
        }
    }

    # Save filter masks to cache
    if (is.null(cached) && length(filter_masks) > 0) {
        dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
        saveRDS(filter_masks, cache_file)
        msg(sprintf("  Saved filter masks: %s", cache_file))
    }
    rm(filter_masks); gc()

    # ==================================================================
    # STEP 3+4: Direction filter + consensus (block-based)
    # ==================================================================
    msg(sprintf("\nStep 3-4: Direction filter (%s) + consensus (%s)...",
                if (filter_dir) "on" else "off", consensus_method))

    adj_consensus <- matrix(0, p, p)
    n_dir_filt <- 0

    for (s in seq(1, p, BS)) {
        e <- min(s + BS - 1, p)
        blocks <- lapply(adj_raw, function(a) a[s:e, , drop = FALSE])
        n_nz   <- Reduce("+", lapply(blocks, function(b) (b != 0) * 1L))
        # For min consensus:    edge must be present in ALL datasets
        # For max consensus:    edge must be present in at least 2
        # For median consensus: edge must be present in at least half the
        #   datasets (otherwise the across-dataset median of |cor| is 0)
        in_2   <- n_nz >= switch(consensus_method,
                                 min    = n_ds,
                                 max    = 2L,
                                 median = as.integer(ceiling(n_ds / 2)))

        if (network_type == "signed") {
            # Signed: require ALL signs positive across datasets.
            # Anti-correlated gene pairs are dropped from the consensus.
            ss <- Reduce("+", lapply(blocks, sign))
            pos_dir <- ss == n_nz
            n_dir_filt <- n_dir_filt + sum(in_2 & !pos_dir)
            valid <- in_2 & pos_dir
        } else if (filter_dir) {
            ss <- Reduce("+", lapply(blocks, sign))
            same_dir <- abs(ss) == n_nz
            n_dir_filt <- n_dir_filt + sum(in_2 & !same_dir)
            valid <- in_2 & same_dir
        } else {
            valid <- in_2
        }

        if (consensus_method == "min") {
            ab <- lapply(blocks, function(b) { a <- abs(b); a[a == 0] <- Inf; a })
            ca <- Reduce(pmin, ab)
            ca[is.infinite(ca)] <- 0
        } else if (consensus_method == "max") {
            ca <- Reduce(pmax, lapply(blocks, abs))
        } else {
            # median consensus: per-edge median of |cor| across ALL datasets,
            # treating an absent edge as 0. With the in_2 (>= half present)
            # gate above, ca > 0 exactly when a majority of datasets share the
            # edge -- a robust middle ground between min (all) and max (any).
            # Memory-light: chunk columns and stack as a 2-D matrix for
            # matrixStats::rowMedians (avoids a p x p x n_ds array).
            nb_r   <- nrow(blocks[[1]]); nb_c <- ncol(blocks[[1]])
            ca     <- matrix(0, nb_r, nb_c)
            CH     <- 256L
            for (cs in seq(1, nb_c, CH)) {
                ce  <- min(cs + CH - 1L, nb_c)
                len <- nb_r * (ce - cs + 1L)
                stak <- vapply(blocks,
                               function(b) as.vector(abs(b[, cs:ce, drop = FALSE])),
                               numeric(len))
                # base-R row medians (avoids a matrixStats dependency; identical
                # result to matrixStats::rowMedians for finite input)
                ca[, cs:ce] <- matrix(apply(stak, 1L, stats::median),
                                      nb_r, ce - cs + 1L)
            }
        }

        adj_consensus[s:e, ] <- ca * valid
    }

    diag(adj_consensus) <- 0
    rownames(adj_consensus) <- fnames
    colnames(adj_consensus) <- fnames
    fstats$n_dir_filtered <- as.integer(n_dir_filt / 2)
    fstats$n_consensus    <- sum(adj_consensus[upper.tri(adj_consensus)] != 0)

    # Collect example edges for scatter plots (before freeing adj_raw)
    example_edges <- .sample_example_edges(adj_raw, adj_consensus, fnames, n_examples = 6)

    rm(adj_raw); gc()

    msg(sprintf("  Direction-filtered: %d", fstats$n_dir_filtered))
    msg(sprintf("  Consensus edges:    %d", fstats$n_consensus))

    if (fstats$n_consensus == 0)
        stop("No consensus edges. Try a different filter_method or disable direction filter.")

    # ==================================================================
    # STEP 5: Soft power
    # ==================================================================
    if (!is.null(soft_power)) {
        msg(sprintf("\nStep 5: Soft power = %g (unsigned)", soft_power))
        for (s in seq(1, p, BS)) {
            e <- min(s + BS - 1, p)
            adj_consensus[s:e, ] <- adj_consensus[s:e, ] ^ soft_power
        }
    } else {
        msg("\nStep 5: No soft power.")
    }

    # ==================================================================
    # STEP 6: Consensus TOM (nb_dist per dataset, then pmax/pmin)
    # ==================================================================
    msg("\nStep 6: Consensus TOM (nb_dist per dataset)...")

    ut_idx <- which(adj_consensus != 0 & upper.tri(adj_consensus), arr.ind = TRUE)
    consensus_filter <- unname(ut_idx)
    colnames(consensus_filter) <- c("cluster_id1", "cluster_id2")
    msg(sprintf("  Filter edges for nb_dist: %s",
                format(nrow(consensus_filter), big.mark = ",")))

    sp <- if (!is.null(soft_power)) soft_power else 1L

    dist_per_ds <- list()
    for (ds in seq_len(n_ds)) {
        msg(sprintf("  DS%d: nb_dist (soft_power=%g, method=%s)...", ds, sp, method))
        dist_per_ds[[ds]] <- netboost::nb_dist(
            filter     = consensus_filter,
            datan      = as.data.frame(datan_list[[ds]]),
            soft_power = sp,
            cores      = cores,
            method     = method
        )
    }

    # Consensus TOM distance:
    #   "min" consensus    -> pmax of distances (weakest link = most conservative)
    #   "max" consensus    -> pmin of distances (strongest link in any dataset)
    #   "median" consensus -> per-pair median distance across datasets (robust)
    msg("  Aggregating consensus TOM distances...")
    if (consensus_method == "min") {
        consensus_dist <- do.call(pmax, dist_per_ds)
    } else if (consensus_method == "max") {
        consensus_dist <- do.call(pmin, dist_per_ds)
    } else {
        consensus_dist <- apply(do.call(cbind, dist_per_ds), 1L, stats::median)
    }
    rm(dist_per_ds); gc()

    msg(sprintf("  Consensus TOM distance range: [%.4f, %.4f]",
                min(consensus_dist), max(consensus_dist)))

    # ==================================================================
    # STEP 7: Clustering (hclust + cutreeDynamic)
    # ==================================================================
    msg("\nStep 7: Clustering (hclust + cutreeDynamic)...")

    # Build full distance matrix from sparse nb_dist (non-filter pairs = 1.0)
    dissTOM <- matrix(1.0, p, p)
    rownames(dissTOM) <- fnames; colnames(dissTOM) <- fnames
    diag(dissTOM) <- 0
    # Vectorised symmetric scatter (matrix indexing): identical result to the
    # per-edge loop but avoids an R-level loop over potentially millions of edges.
    dissTOM[consensus_filter]            <- consensus_dist
    dissTOM[consensus_filter[, c(2, 1)]] <- consensus_dist

    geneTree    <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                  deepSplit = deepSplit,
                                  pamRespectsDendro = FALSE,
                                  minClusterSize = min_cluster_size)
    names(dynamicMods) <- fnames
    rm(dissTOM); gc()

    rm(adj_consensus, consensus_dist, consensus_filter); gc()

    n_initial <- length(unique(dynamicMods[dynamicMods > 0]))
    msg(sprintf("  Initial modules: %d", n_initial))
    if (n_initial > 0) {
        tab <- sort(table(dynamicMods[dynamicMods > 0]), decreasing = TRUE)
        msg(sprintf("  Sizes: %s", paste(tab, collapse = ", ")))
    }
    if (n_initial == 0) {
        message("No modules found.")
        return(invisible(list(
            colors = setNames(dynamicMods, fnames),
            n_modules_final = 0)))
    }

    # ==================================================================
    # STEP 8: Module eigengenes (aligned to mean expression per dataset)
    # ==================================================================
    msg("\nStep 8: Module eigengenes...")
    ME_per_ds <- list()
    for (ds in seq_len(n_ds)) {
        dat    <- as.matrix(datan_list[[ds]])
        ME_res <- nb_moduleEigengenes(dat, colors = dynamicMods,
                                       n_pc = n_pc, robust = robust_PCs,
                                       nb_min_varExpl = nb_min_varExpl,
                                       scale = scale)
        ME_per_ds[[ds]] <- ME_res$eigengenes
    }

    # Align MEs to mean expression of module genes (ensures consistent
    # sign orientation across datasets before cross-dataset PC comparison)
    msg("  Aligning MEs to mean expression...")
    mods_align <- sort(unique(dynamicMods[dynamicMods > 0]))
    for (ds in seq_len(n_ds)) {
        dat <- as.matrix(datan_list[[ds]])
        for (mod in mods_align) {
            me_col <- paste0("ME", mod)
            if (!(me_col %in% colnames(ME_per_ds[[ds]]))) next
            genes <- names(dynamicMods[dynamicMods == mod])
            mean_expr <- rowMeans(dat[, genes, drop = FALSE], na.rm = TRUE)
            if (.neg_cor(ME_per_ds[[ds]][, me_col], mean_expr))
                ME_per_ds[[ds]][, me_col] <- -ME_per_ds[[ds]][, me_col]
        }
    }

    # ==================================================================
    # STEP 9: Iterative cross-dataset module merging
    #
    # Iterative hclust: build ME dissimilarity dendrogram, cut at threshold,
    # merge clusters, recompute eigengenes, repeat until stable.
    # ==================================================================
    merge_cor_thres <- 1 - ME_diss_thres
    msg(sprintf("\nStep 9: Iterative hclust merging (diss thres = %.2f, cor > %.2f, %s datasets)...",
                ME_diss_thres, merge_cor_thres, merge_criterion))

    merged_colors <- dynamicMods
    ME_cor_ds <- NULL; ME_tree <- NULL
    n_merges <- 0
    n_iterations <- 0

    repeat {
        n_iterations <- n_iterations + 1
        mods <- sort(unique(merged_colors[merged_colors > 0]))
        if (length(mods) < 2) break

        # Compute MEs per dataset
        ME_ds <- list()
        for (ds in seq_len(n_ds)) {
            dat <- as.matrix(datan_list[[ds]])
            ME_res <- nb_moduleEigengenes(dat, colors = merged_colors,
                                           n_pc = n_pc, robust = robust_PCs,
                                           nb_min_varExpl = nb_min_varExpl,
                                           scale = scale)
            ME_ds[[ds]] <- ME_res$eigengenes
        }

        # Align MEs to mean expression
        for (ds in seq_len(n_ds)) {
            dat <- as.matrix(datan_list[[ds]])
            for (mod in mods) {
                me_col <- paste0("ME", mod)
                if (!(me_col %in% colnames(ME_ds[[ds]]))) next
                genes <- names(merged_colors[merged_colors == mod])
                mean_expr <- rowMeans(dat[, genes, drop = FALSE], na.rm = TRUE)
                if (.neg_cor(ME_ds[[ds]][, me_col], mean_expr))
                    ME_ds[[ds]][, me_col] <- -ME_ds[[ds]][, me_col]
            }
        }

        # Consensus ME dissimilarity
        me_cols <- paste0("ME", mods)
        me_cols <- me_cols[me_cols %in% colnames(ME_ds[[1]])]
        if (length(me_cols) < 2) break

        ME_cor_ds <- lapply(ME_ds, function(MEs)
            abs(WGCNA::cor(MEs[, me_cols, drop = FALSE],
                           use = "pairwise.complete.obs")))

        if (merge_criterion == "all") {
            ME_cor_cross <- Reduce(pmin, ME_cor_ds)
        } else {
            ME_cor_cross <- Reduce(pmax, ME_cor_ds)
        }
        ME_diss <- 1 - ME_cor_cross
        diag(ME_diss) <- 0
        rownames(ME_diss) <- me_cols
        colnames(ME_diss) <- me_cols

        # Hierarchical clustering of MEs and cut at threshold
        ME_tree <- hclust(as.dist(ME_diss), method = "average")
        cut <- cutree(ME_tree, h = ME_diss_thres)

        # Merge clusters from the cut
        merged_any <- FALSE
        for (cl in unique(cut)) {
            members <- me_cols[cut == cl]
            if (length(members) < 2) next
            mod_ids <- as.integer(sub("^ME", "", members))
            keep <- min(mod_ids)
            for (drop in sort(mod_ids[mod_ids != keep])) {
                merged_colors[merged_colors == drop] <- keep
                n_merges <- n_merges + 1
                merged_any <- TRUE
                msg(sprintf("  Iter %d, merge %d: M%d + M%d -> M%d",
                            n_iterations, n_merges, keep, drop, keep))
            }
        }
        if (!merged_any) break
    }

    n_merged <- length(unique(merged_colors[merged_colors > 0]))
    msg(sprintf("  After merging: %d modules (%d merges, %d iterations)",
                n_merged, n_merges, n_iterations))

    # ==================================================================
    # STEP 10: Module retention (size + cross-dataset PC1 cor)
    # ==================================================================
    msg(sprintf("\nStep 10: Retention (%s, PC1 cor >= %.2f, minSz = %d)...",
                module_retention, min_pc1_cor, min_cluster_size))
    all_mods <- sort(unique(merged_colors[merged_colors > 0]))
    retain   <- integer(0)
    # Track which datasets need ME sign flip per module (due to anti-correlated PC1)
    flip_ds <- list()  # flip_ds[[mod]] = vector of dataset indices to flip
    for (mod in all_mods) {
        sz <- sum(merged_colors == mod)
        if (sz < min_cluster_size) {
            msg(sprintf("  Dropping module %d (size %d < %d)", mod, sz, min_cluster_size))
            next
        }

        # Cross-dataset PC1 loading correlation (aligned to mean expression)
        genes <- names(merged_colors[merged_colors == mod])
        pc1_list <- list()
        for (ds in seq_len(n_ds)) {
            ex <- datan_list[[ds]][, genes, drop = FALSE]
            if (ncol(ex) >= 2) {
                # robust_PCs: rank-transform so PC1 loadings are rank-based,
                # consistent with nb_moduleEigengenes(robust=) and the rotation
                # stored for nb_consensus_transfer().
                exm <- if (robust_PCs) apply(ex, 2L, rank) else as.matrix(ex)
                pc <- prcomp(exm, center = TRUE, scale. = scale)
                w  <- pc$rotation[, 1]
                if (.neg_cor(pc$x[, 1], rowMeans(exm, na.rm = TRUE)))
                    w <- -w
                pc1_list[[ds]] <- w
            }
        }

        cors <- c()
        for (a in seq_len(n_ds - 1)) for (b in (a + 1):n_ds) {
            if (!is.null(pc1_list[[a]]) && !is.null(pc1_list[[b]])) {
                v <- !is.na(pc1_list[[a]]) & !is.na(pc1_list[[b]])
                if (sum(v) >= 2) {
                    r <- suppressWarnings(cor(pc1_list[[a]][v], pc1_list[[b]][v]))
                    if (!is.na(r)) cors <- c(cors, r)   # skip zero-variance pairs
                }
            }
        }

        # Apply retention criterion: "all" = min, "any" = max
        if (length(cors) > 0) {
            me_cor <- switch(module_retention, all = min(cors), any = max(cors), median = stats::median(cors))
        } else {
            me_cor <- NA_real_
        }

        # Handle flipped eigengenes: if PC1 weights are strongly negatively
        # correlated (|cor| >= threshold but sign is wrong), the ME orientation
        # flipped between cohorts. Fix by flipping PC weights of the smaller
        # dataset(s) to align with the largest, then recompute correlation.
        if (!is.na(me_cor) && me_cor <= -min_pc1_cor) {
            msg(sprintf("  Module %d: strong negative PC1 cor (%.3f), flipping smaller datasets...",
                        mod, me_cor))

            # Find the largest dataset index
            ds_sizes <- vapply(seq_len(n_ds), function(ds) nrow(datan_list[[ds]]), integer(1))
            ref_ds <- which.max(ds_sizes)

            # Flip PC weights of non-reference datasets that are negatively
            # correlated with the reference
            flipped <- c()
            for (ds in seq_len(n_ds)) {
                if (ds == ref_ds || is.null(pc1_list[[ds]]) || is.null(pc1_list[[ref_ds]])) next
                v <- !is.na(pc1_list[[ref_ds]]) & !is.na(pc1_list[[ds]])
                if (sum(v) >= 2) {
                    r <- suppressWarnings(cor(pc1_list[[ref_ds]][v], pc1_list[[ds]][v]))
                    if (isTRUE(r < 0)) {
                        pc1_list[[ds]] <- -pc1_list[[ds]]
                        flipped <- c(flipped, ds)
                        msg(sprintf("    Flipped dataset %d (was r=%.3f with ref dataset %d)",
                                    ds, r, ref_ds))
                    }
                }
            }
            if (length(flipped) > 0) flip_ds[[as.character(mod)]] <- flipped

            # Recompute pairwise correlations after flipping
            cors <- c()
            for (a in seq_len(n_ds - 1)) for (b in (a + 1):n_ds) {
                if (!is.null(pc1_list[[a]]) && !is.null(pc1_list[[b]])) {
                    v <- !is.na(pc1_list[[a]]) & !is.na(pc1_list[[b]])
                    if (sum(v) >= 2) {
                        r <- suppressWarnings(cor(pc1_list[[a]][v], pc1_list[[b]][v]))
                        if (!is.na(r)) cors <- c(cors, r)
                    }
                }
            }
            if (length(cors) > 0) {
                me_cor <- switch(module_retention, all = min(cors), any = max(cors), median = stats::median(cors))
            }
            msg(sprintf("  Module %d: PC1 cor after flip = %.3f", mod, me_cor))
        }

        if (!is.na(me_cor) && me_cor >= min_pc1_cor) {
            retain <- c(retain, mod)
            msg(sprintf("  Retaining module %d (size %d, PC1 cor = %.3f)",
                        mod, sz, me_cor))
        } else {
            msg(sprintf("  Greying module %d (size %d, PC1 cor = %.3f < %.2f)",
                        mod, sz, if (is.na(me_cor)) NA else me_cor, min_pc1_cor))
        }
    }

    final_colors <- merged_colors
    final_colors[!(final_colors %in% retain)] <- 0
    names(final_colors) <- fnames
    n_final <- length(retain)
    msg(sprintf("  Final: %d modules", n_final))

    # ==================================================================
    # Recompute MEs with final colors (aligned to mean expression)
    # ==================================================================
    ME_final <- list()
    var_explained_ds <- list()
    for (ds in seq_len(n_ds)) {
        dat    <- as.matrix(datan_list[[ds]])
        ME_res <- nb_moduleEigengenes(dat, colors = final_colors,
                                       n_pc = n_pc, robust = robust_PCs,
                                       nb_min_varExpl = nb_min_varExpl,
                                       scale = scale)
        ME_final[[ds]] <- ME_res$eigengenes
        var_explained_ds[[ds]] <- ME_res$var_explained
    }

    # Align final MEs to mean expression of module genes
    # For modules where PC1 was flipped between datasets (detected in Step 10),
    # align all datasets to the reference dataset (largest) instead of aligning
    # each independently to mean expression (which can produce inconsistent signs).
    ds_sizes <- vapply(seq_len(n_ds), function(ds) nrow(datan_list[[ds]]), integer(1))
    ref_ds <- which.max(ds_sizes)

    for (mod in sort(unique(final_colors[final_colors > 0]))) {
        me_col <- paste0("ME", mod)
        genes <- names(final_colors[final_colors == mod])
        mod_key <- as.character(mod)

        if (mod_key %in% names(flip_ds)) {
            # Flipped module: align reference dataset to mean expression,
            # then align other datasets to the reference ME via PC weight correlation
            dat_ref <- as.matrix(datan_list[[ref_ds]])
            if (me_col %in% colnames(ME_final[[ref_ds]])) {
                mean_expr_ref <- rowMeans(dat_ref[, genes, drop = FALSE], na.rm = TRUE)
                if (.neg_cor(ME_final[[ref_ds]][, me_col], mean_expr_ref))
                    ME_final[[ref_ds]][, me_col] <- -ME_final[[ref_ds]][, me_col]
            }

            # For flipped datasets: flip the ME sign
            for (ds in flip_ds[[mod_key]]) {
                if (!(me_col %in% colnames(ME_final[[ds]]))) next
                dat_ds <- as.matrix(datan_list[[ds]])
                mean_expr_ds <- rowMeans(dat_ds[, genes, drop = FALSE], na.rm = TRUE)
                # First align to mean expression (standard)
                if (.neg_cor(ME_final[[ds]][, me_col], mean_expr_ds))
                    ME_final[[ds]][, me_col] <- -ME_final[[ds]][, me_col]
                # Then flip (since we know PC weights are anti-correlated with ref)
                ME_final[[ds]][, me_col] <- -ME_final[[ds]][, me_col]
                msg(sprintf("  Flipped ME%d sign for dataset %d (anti-correlated PC weights with ref)",
                            mod, ds))
            }

            # Non-flipped, non-reference datasets: standard alignment
            for (ds in setdiff(seq_len(n_ds), c(ref_ds, flip_ds[[mod_key]]))) {
                if (!(me_col %in% colnames(ME_final[[ds]]))) next
                dat_ds <- as.matrix(datan_list[[ds]])
                mean_expr_ds <- rowMeans(dat_ds[, genes, drop = FALSE], na.rm = TRUE)
                if (.neg_cor(ME_final[[ds]][, me_col], mean_expr_ds))
                    ME_final[[ds]][, me_col] <- -ME_final[[ds]][, me_col]
            }
        } else {
            # Standard alignment: each dataset aligned to its own mean expression
            for (ds in seq_len(n_ds)) {
                if (!(me_col %in% colnames(ME_final[[ds]]))) next
                dat_ds <- as.matrix(datan_list[[ds]])
                mean_expr_ds <- rowMeans(dat_ds[, genes, drop = FALSE], na.rm = TRUE)
                if (.neg_cor(ME_final[[ds]][, me_col], mean_expr_ds))
                    ME_final[[ds]][, me_col] <- -ME_final[[ds]][, me_col]
            }
        }
    }

    # ==================================================================
    # PC weights + rotation (sign-aligned to mean expression per dataset)
    # ==================================================================
    msg("\nPC weights...")
    mods_final <- sort(unique(final_colors[final_colors > 0]))
    pc1_all    <- list()
    rotation_list <- list()

    for (mod in mods_final) {
        genes  <- names(final_colors[final_colors == mod])
        mod_key <- as.character(mod)
        w_list <- list()
        rot_list_ds <- list()
        for (i in seq_len(n_ds)) {
            ex <- datan_list[[i]][, genes, drop = FALSE]
            if (ncol(ex) >= 2) {
                # robust_PCs: rank-transform the discovery data so the stored
                # rotation matches the rank-based projection nb_consensus_transfer()
                # applies (otherwise robust_PCs=TRUE yields a mismatched projection).
                exm <- if (robust_PCs) apply(ex, 2L, rank) else as.matrix(ex)
                pc <- prcomp(exm, center = TRUE, scale. = scale)
                n_avail <- min(n_pc, ncol(pc$rotation))
                rot_ds <- pc$rotation[, seq_len(n_avail), drop = FALSE]
                # PC1 sign: align to mean expression (one cor() call, NA-guarded)
                if (.neg_cor(pc$x[, 1], rowMeans(exm, na.rm = TRUE)))
                    rot_ds[, 1] <- -rot_ds[, 1]
                # PC2+ have no mean-expression anchor; orient deterministically so
                # the largest-magnitude loading is positive (reproducible across
                # datasets instead of prcomp's arbitrary per-dataset sign).
                if (n_avail >= 2) for (k in 2:n_avail) {
                    if (rot_ds[which.max(abs(rot_ds[, k])), k] < 0)
                        rot_ds[, k] <- -rot_ds[, k]
                }
                w <- rot_ds[, 1]

                # For flipped modules: flip PC weights of affected datasets
                # so they align with the reference before averaging
                if (mod_key %in% names(flip_ds) && i %in% flip_ds[[mod_key]]) {
                    w <- -w
                    rot_ds[, 1] <- -rot_ds[, 1]
                }
            } else {
                w <- setNames(rep(NA_real_, length(genes)), genes)
                rot_ds <- matrix(NA_real_, length(genes), n_pc,
                                  dimnames = list(genes, paste0("PC", seq_len(n_pc))))
            }
            w_list[[i]] <- w
            rot_list_ds[[i]] <- rot_ds
        }

        df <- data.frame(gene_id = genes, module = mod, stringsAsFactors = FALSE)
        for (i in seq_len(n_ds))
            df[[paste0("PC1_dataset", i)]] <- w_list[[i]][genes]
        pc1c <- grep("^PC1_dataset", names(df))
        df$PC1_mean <- rowMeans(df[, pc1c, drop = FALSE], na.rm = TRUE)

        cors <- c()
        for (a in seq_len(n_ds - 1)) for (b in (a + 1):n_ds) {
            v <- !is.na(w_list[[a]]) & !is.na(w_list[[b]])
            if (sum(v) >= 2) {
                r <- suppressWarnings(cor(w_list[[a]][v], w_list[[b]][v]))
                if (!is.na(r)) cors <- c(cors, r)
            }
        }
        df$PC1_cor <- if (length(cors) > 0) mean(cors) else NA_real_
        df$flipped_datasets <- if (mod_key %in% names(flip_ds))
            paste(flip_ds[[mod_key]], collapse = ",") else NA_character_
        pc1_all[[mod_key]] <- df

        # Mean rotation across datasets for this module
        mean_rot <- Reduce("+", rot_list_ds) / length(rot_list_ds)
        rownames(mean_rot) <- genes
        rotation_list[[mod_key]] <- mean_rot
    }
    if (length(pc1_all) > 0) {
        all_weights <- do.call(rbind, pc1_all)
        rownames(all_weights) <- NULL
    } else {
        all_weights <- data.frame(gene_id = character(0), module = integer(0),
                                   PC1_mean = numeric(0), PC1_cor = numeric(0))
    }

    # Build combined rotation matrix (all genes, n_pc columns)
    rotation <- matrix(0, nrow = p, ncol = n_pc,
                        dimnames = list(fnames, paste0("PC", seq_len(n_pc))))
    for (mod_chr in names(rotation_list)) {
        rot <- rotation_list[[mod_chr]]
        n_avail <- min(ncol(rot), n_pc)
        rotation[rownames(rot), seq_len(n_avail)] <- rot[, seq_len(n_avail)]
    }

    # ==================================================================
    # Build results
    # ==================================================================
    results <- list(
        # --- nb_summary-compatible fields ---
        names             = fnames,
        colors            = final_colors,
        MEs               = ME_final[[1]],
        var_explained     = var_explained_ds,
        rotation          = rotation,
        filter            = filter_per_ds,

        # --- Consensus-specific fields ---
        colors_initial    = dynamicMods,
        colors_merged     = merged_colors,
        modules_retained  = retain,
        geneTree          = geneTree,
        ME_all_datasets   = ME_final,
        ME_cor_per_dataset = ME_cor_ds,
        ME_tree           = ME_tree,
        pc1_weights       = all_weights,
        filter_stats      = fstats,
        example_edges     = example_edges,
        params = list(
            filter_method = filter_method,
            network_type = network_type, stepno = stepno, n_pc = n_pc,
            robust_PCs = robust_PCs, nb_min_varExpl = nb_min_varExpl,
            soft_power = soft_power, consensus_method = consensus_method,
            filter_dir = filter_dir, min_cluster_size = min_cluster_size,
            ME_diss_thres = ME_diss_thres, merge_criterion = merge_criterion,
            module_retention = module_retention, min_pc1_cor = min_pc1_cor,
            deepSplit = deepSplit,
            method = method, scale = scale, cores = cores
        ),
        n_datasets        = n_ds,
        n_features        = p,
        n_modules_initial = n_initial,
        n_modules_merged  = n_merged,
        n_modules_final   = n_final
    )

    # ==================================================================
    # PDF report
    # ==================================================================
    if (!is.null(pdf_report)) {
        msg(sprintf("\nGenerating report: %s", pdf_report))
        .consensus_report(results, datan_list, pdf_report)
    }

    # ==================================================================
    # Summary
    # ==================================================================
    msg(paste0("\n", paste(rep("=", 60), collapse = "")))
    msg("SUMMARY")
    msg(paste(rep("=", 60), collapse = ""))
    msg(sprintf("  Features:    %d", p))
    msg(sprintf("  In modules:  %d (%.1f%%)",
                sum(final_colors != 0), 100 * sum(final_colors != 0) / p))
    msg(sprintf("  Unassigned:  %d", sum(final_colors == 0)))
    n_greyed <- n_merged - n_final
    msg(sprintf("  Modules: %d (clustering) -> %d (merging) -> %d (retention)",
                n_initial, n_merged, n_final))
    if (n_greyed > 0)
        msg(sprintf("  Greyed:  %d modules (PC1 cor < %.2f)",
                    n_greyed, min_pc1_cor))
    if (n_merges > 0)
        msg(sprintf("  Merged:  %d times (|ME cor| >= %.2f in %s datasets)",
                    n_merges, merge_cor_thres, merge_criterion))
    if (n_final > 0) {
        sz <- table(final_colors[final_colors != 0])
        msg(sprintf("  Sizes: %d - %d (mean %.1f)", min(sz), max(sz), mean(sz)))
    }
    if (nrow(all_weights) > 0 && any(!is.na(all_weights$PC1_cor))) {
        pc1_per_mod <- unique(all_weights[, c("module", "PC1_cor")])
        n_above <- sum(pc1_per_mod$PC1_cor >= min_pc1_cor, na.rm = TRUE)
        n_total_mod <- nrow(pc1_per_mod)
        msg(sprintf("  PC1 cor: mean=%.3f [%.3f, %.3f]",
                    mean(pc1_per_mod$PC1_cor, na.rm = TRUE),
                    min(pc1_per_mod$PC1_cor, na.rm = TRUE),
                    max(pc1_per_mod$PC1_cor, na.rm = TRUE)))
        msg(sprintf("  PC1 concordance >= %.2f: %d / %d modules",
                    min_pc1_cor, n_above, n_total_mod))
    }
    msg(paste(rep("=", 60), collapse = ""))
    invisible(results)
}

# =============================================================================
# TRANSFER FUNCTION
# =============================================================================

#' Transfer consensus module assignments and eigengenes to new data
#'
#' Projects new samples onto the consensus modules discovered by
#' \code{\link{nb_consensus}} using the sign-aligned PC1 rotation weights from
#' the discovery cohorts (no per-cohort sign flipping), giving consistently
#' oriented module eigengenes across validation datasets.
#'
#' @param nb_consensus_result Result list returned by \code{\link{nb_consensus}}.
#' @param new_data Data frame or matrix (rows = samples, cols = features).
#'   Feature names are matched against the consensus result.
#' @param scale Logical; scale and center \code{new_data} before projection.
#' @param robust_PCs Logical; rank-transform \code{new_data} before projection.
#' @return A data frame of projected module eigengenes (samples x modules).
#' @seealso \code{\link{nb_consensus}}
#' @examples
#' set.seed(1)
#' mk <- function(n) {
#'     m <- matrix(rnorm(n * 12), n, 12)
#'     m[, 1:4] <- m[, 1:4] + rnorm(n)      # correlated block 1
#'     m[, 5:8] <- m[, 5:8] + rnorm(n)      # correlated block 2
#'     colnames(m) <- paste0("g", 1:12)
#'     m
#' }
#' res <- nb_consensus(list(mk(40), mk(45)), filter_method = "pearson",
#'                     min_cluster_size = 3, stepno = 5L, verbose = FALSE)
#' projected <- nb_consensus_transfer(res, mk(20))
#' head(projected)
#' @export
nb_consensus_transfer <- function(nb_consensus_result,
                                   new_data,
                                   scale = FALSE,
                                   robust_PCs = FALSE) {
    res   <- nb_consensus_result
    fnames <- res$names
    colors <- res$colors
    rot    <- res$rotation

    stopifnot(!is.null(fnames), !is.null(colors), !is.null(rot))

    # Match features
    new_data <- as.matrix(new_data)
    common   <- intersect(colnames(new_data), fnames)
    if (length(common) == 0)
        stop("No overlapping feature names between new_data and consensus result.")
    if (length(common) < length(fnames))
        message(sprintf("  nb_consensus_transfer: %d / %d features matched",
                        length(common), length(fnames)))

    dat <- new_data[, common, drop = FALSE]
    if (scale) dat <- scale(dat)
    if (robust_PCs) dat <- apply(dat, 2, rank)

    mods <- sort(unique(colors[colors != 0]))
    # Guard the no-module case: paste0("ME", integer(0)) returns "ME" (length 1,
    # the zero-length-argument recycling rule), which would mismatch ncol = 0.
    me_names <- if (length(mods) > 0) paste0("ME", mods) else character(0)
    MEs  <- matrix(NA_real_, nrow = nrow(dat), ncol = length(mods),
                    dimnames = list(rownames(dat), me_names))

    # The consensus rotation weights are already sign-aligned to mean
    # expression in the discovery cohorts (Step 8 of nb_consensus).
    # We use them directly WITHOUT per-cohort sign flipping to ensure
    # consistent ME orientation across all validation cohorts.

    for (j in seq_along(mods)) {
        mod   <- mods[j]
        genes <- names(colors[colors == mod])
        genes <- intersect(genes, common)
        if (length(genes) < 2) next

        # Project using consensus rotation (PC1) -- sign is fixed by discovery
        w <- rot[genes, 1]
        w <- w[!is.na(w)]
        genes_use <- intersect(names(w), colnames(dat))
        if (length(genes_use) < 2) next

        scores <- dat[, genes_use, drop = FALSE] %*% w[genes_use]
        MEs[, j] <- scores
    }

    as.data.frame(MEs)
}

# =============================================================================
# Helper: sample example edges for report scatter plots
# =============================================================================

#' Sample top consensus edges for report scatter plots
#' @noRd
.sample_example_edges <- function(adj_raw, adj_consensus, fnames, n_examples = 6) {
    p <- length(fnames)
    # Get upper-triangle indices of strong consensus edges
    ut <- which(upper.tri(adj_consensus) & adj_consensus != 0, arr.ind = TRUE)
    if (nrow(ut) == 0) return(NULL)

    vals <- abs(adj_consensus[ut])
    # Top edges by consensus strength
    top_idx <- utils::head(order(vals, decreasing = TRUE), n_examples)
    top_edges <- ut[top_idx, , drop = FALSE]

    lapply(seq_len(nrow(top_edges)), function(k) {
        r <- top_edges[k, 1]; cc <- top_edges[k, 2]
        per_ds <- vapply(adj_raw, function(a) a[r, cc], numeric(1))
        list(gene1 = fnames[r], gene2 = fnames[cc],
             consensus_val = adj_consensus[r, cc],
             per_dataset_val = per_ds)
    })
}

# =============================================================================
# PDF report
# =============================================================================

#' Render the multi-page consensus PDF report
#' @noRd
.consensus_report <- function(results, datan_list, pdf_path) {

    dir.create(dirname(pdf_path), showWarnings = FALSE, recursive = TRUE)
    pdf(pdf_path, width = 11, height = 8.5)
    on.exit(dev.off(), add = TRUE)

    fs     <- results$filter_stats
    pa     <- results$params
    n_ds   <- results$n_datasets
    fc     <- results$colors

    # ---- Page 1: Filtering overview ----
    par(mfrow = c(2, 2), mar = c(7, 5, 4, 2), oma = c(0, 0, 0, 0))

    # 1a: Edge counts per dataset + consensus
    edge_vals <- c(fs$n_after_filter, fs$n_consensus)
    edge_labs <- c(paste0("DS", seq_len(n_ds), " filtered"), "Consensus")
    ymax <- max(edge_vals) * 1.15
    bp <- barplot(edge_vals, names.arg = edge_labs, col = c(rep("#4575b4", n_ds), "#91cf60"),
                  main = "Edge Counts", ylab = "Edges", las = 2,
                  cex.names = 0.8, ylim = c(0, ymax))
    text(bp, edge_vals, labels = format(edge_vals, big.mark = ","),
         pos = 3, cex = 0.7)

    # 1b: Edge density per dataset (% of total possible edges)
    pct <- fs$n_after_filter / fs$n_total * 100
    ymax_pct <- max(pct) * 1.25
    bp2 <- barplot(pct, names.arg = paste("DS", seq_len(n_ds)),
            col = c("#2166ac", "#b2182b", "#1b7837", "#762a83")[seq_len(n_ds)],
            main = "Edge Density per Dataset",
            ylab = sprintf("%%  of total possible (%s)", format(fs$n_total, big.mark = ",")),
            ylim = c(0, ymax_pct), cex.names = 0.8)
    text(bp2, pct, labels = sprintf("%.1f%%", pct), pos = 3, cex = 0.8)

    # 1c: Edge fate pie
    cons   <- fs$n_consensus
    dir_f  <- fs$n_dir_filtered
    no_sig <- fs$n_total - max(fs$n_after_filter)
    other  <- fs$n_total - cons - dir_f - no_sig
    vals   <- c(cons, dir_f, max(other, 0), no_sig)
    labs   <- c(sprintf("Consensus\n(%s)", format(cons, big.mark = ",")),
                sprintf("Dir-filtered\n(%s)", format(dir_f, big.mark = ",")),
                sprintf("In 1 DS only\n(%s)", format(max(other, 0), big.mark = ",")),
                sprintf("Not signif.\n(%s)", format(no_sig, big.mark = ",")))
    par(mar = c(4, 2, 4, 2))
    pie(vals, labels = labs, col = c("#91cf60", "#fc8d59", "#fee08b", "#d9d9d9"),
        main = "Edge Fate", cex = 0.8)

    # 1d: Parameters text
    par(mar = c(2, 2, 3, 2))
    plot.new()
    title("Parameters", line = 0.5)
    txt <- c(
        sprintf("Filter: %s (nb_filter, stepno=%d)", pa$filter_method, pa$stepno),
        sprintf("TOM: nb_dist consensus (%s)", pa$consensus_method),
        sprintf("Clustering: hclust + cutreeDynamic"),
        sprintf("Adjacency: %s", pa$method),
        sprintf("Direction filter: %s", pa$filter_dir),
        sprintf("Consensus: %s", pa$consensus_method),
        sprintf("Soft power: %s", ifelse(is.null(pa$soft_power), "none", pa$soft_power)),
        sprintf("ME merge thres: %.2f (cor >= %.2f, %s)",
                pa$ME_diss_thres, 1 - pa$ME_diss_thres, pa$merge_criterion),
        sprintf("Retention: %s, PC1 cor >= %.2f", pa$module_retention, pa$min_pc1_cor),
        sprintf("deepSplit: %d | minSz: %d", pa$deepSplit, pa$min_cluster_size),
        "",
        sprintf("Features: %d", results$n_features),
        sprintf("In modules: %d (%.1f%%)",
                sum(fc != 0), 100 * sum(fc != 0) / results$n_features),
        sprintf("Modules: %d (clustering) -> %d (merging) -> %d (retention)",
                results$n_modules_initial, results$n_modules_merged,
                results$n_modules_final),
        {
            pw <- results$pc1_weights
            if (!is.null(pw) && nrow(pw) > 0 && any(!is.na(pw$PC1_cor))) {
                pc1_mod <- unique(pw[, c("module", "PC1_cor")])
                sprintf("PC1 cor: mean=%.3f [%.3f, %.3f]",
                        mean(pc1_mod$PC1_cor, na.rm = TRUE),
                        min(pc1_mod$PC1_cor, na.rm = TRUE),
                        max(pc1_mod$PC1_cor, na.rm = TRUE))
            } else ""
        }
    )
    for (k in seq_along(txt))
        text(0.05, 0.95 - k * 0.08, txt[k], adj = 0, cex = 0.85, family = "mono")

    # ---- Page 2: ME dendrogram ----
    # Only when there are >= 3 modules; a 1-2 leaf tree carries no information
    # and would otherwise emit a blank / "not available" placeholder page.
    if (!is.null(results$ME_tree) && length(results$ME_tree$order) >= 3) {
        par(mfrow = c(1, 1), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
        tryCatch({
            plot(results$ME_tree, main = "Module Eigengene Dendrogram (pre-merge)",
                 xlab = "", ylab = "Dissimilarity (1 - min|cor|)", sub = "",
                 cex = 0.8)
            abline(h = pa$ME_diss_thres, col = "red", lty = 2, lwd = 2)
            mtext(sprintf("merge threshold = %.2f  (cor = %.2f)",
                           pa$ME_diss_thres, 1 - pa$ME_diss_thres),
                  side = 3, line = -1.5, col = "red", cex = 0.9)
        }, error = function(e) {
            plot.new()
            text(0.5, 0.5, sprintf("Dendrogram not available\n(%d modules)",
                                    results$n_modules_initial), cex = 1.2)
        })
    }

    # ---- Page 3: ME correlation heatmaps per dataset ----
    if (!is.null(results$ME_cor_per_dataset) && length(results$ME_cor_per_dataset) > 0) {
        par(mfrow = c(1, n_ds), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
        for (ds in seq_len(n_ds)) {
            mc <- results$ME_cor_per_dataset[[ds]]
            lbl <- sub("^ME", "", rownames(mc))
            image(seq_len(nrow(mc)), seq_len(ncol(mc)), mc,
                  col = colorRampPalette(c("white", "orange", "red"))(50),
                  zlim = c(0, 1), axes = FALSE,
                  main = sprintf("DS%d: |ME cor|", ds))
            axis(1, seq_len(nrow(mc)), lbl, las = 2, cex.axis = 0.5)
            axis(2, seq_len(ncol(mc)), lbl, las = 2, cex.axis = 0.5)
        }
    }

    # ---- Page 4: Module sizes ----
    par(mfrow = c(1, 1), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
    tab <- sort(table(fc[fc != 0]), decreasing = TRUE)
    if (length(tab) > 0) {
        ymax_tab <- max(tab) * 1.1
        bp <- barplot(tab, main = "Final Module Sizes", xlab = "Module",
                      ylab = "Genes", col = "steelblue", las = 2,
                      ylim = c(0, ymax_tab), cex.names = 0.8)
        text(bp, tab, labels = tab, pos = 3, cex = 0.65)
    } else {
        plot.new()
        text(0.5, 0.5, "No modules retained", cex = 1.5)
    }

    # ---- Page 5: Gene-gene scatters (top consensus edges) ----
    ex <- results$example_edges
    if (!is.null(ex) && length(ex) > 0) {
        n_ex <- min(length(ex), 6)
        par(mfrow = c(n_ex, n_ds), mar = c(4, 4, 2.5, 1.5),
            oma = c(1, 1, 2, 1), mgp = c(2, 0.6, 0))
        for (k in seq_len(n_ex)) {
            g1 <- ex[[k]]$gene1; g2 <- ex[[k]]$gene2
            for (ds in seq_len(n_ds)) {
                x <- datan_list[[ds]][, g1]
                y <- datan_list[[ds]][, g2]
                r <- round(cor(x, y, use = "complete.obs"), 3)
                plot(x, y, pch = 19, cex = 0.7,
                     col = adjustcolor("steelblue", 0.7),
                     main = sprintf("DS%d: r=%.3f", ds, r),
                     xlab = g1, ylab = g2, cex.main = 0.8, cex.lab = 0.65,
                     cex.axis = 0.7)
                abline(lm(y ~ x), col = "red", lty = 2)
            }
        }
        mtext("Top Consensus Edges: Gene-Gene Scatters", outer = TRUE,
              side = 3, line = 0.3, cex = 1.0)
    }

    # ---- Page 6: PC1 weight scatter across datasets ----
    wt <- results$pc1_weights
    if (n_ds == 2 && !is.null(wt) && NROW(wt) > 0) {
        mods_plot <- sort(unique(wt$module))
        n_plot <- min(9, length(mods_plot))
        nr <- ceiling(n_plot / 3)
        par(mfrow = c(nr, 3), mar = c(4, 4, 3, 2), oma = c(1, 1, 2, 1))
        for (mod in mods_plot[seq_len(n_plot)]) {
            mw <- wt[wt$module == mod, ]
            w1 <- mw$PC1_dataset1; w2 <- mw$PC1_dataset2
            r  <- round(cor(w1, w2, use = "complete.obs"), 3)
            plot(w1, w2, pch = 19, cex = 0.5,
                 col = adjustcolor("darkgreen", 0.6),
                 main = sprintf("Mod %s (r=%.3f)", mod, r),
                 xlab = "DS1 weight", ylab = "DS2 weight",
                 cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7)
            abline(0, 1, col = "grey50", lty = 2)
            abline(lm(w2 ~ w1), col = "red")
        }
        mtext("PC1 Weights Across Datasets", outer = TRUE,
              side = 3, line = 0.3, cex = 1.0)
    }

    message(sprintf("  Report: %s", pdf_path))
}
