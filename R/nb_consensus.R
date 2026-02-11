###############################################################################
# nb_consensus.R
# Multi-dataset consensus network analysis for the netboost ecosystem
#
# 10-step workflow:
#   1+2. Filtered adjacency per dataset (nb_filter for
#        boosting/pearson/spearman/kendall, or skip)
#   3+4. Direction filter + consensus at the adjacency level
#   5.   Soft thresholding
#   6.   TOM (WGCNA::TOMsimilarity on consensus adjacency)
#   7.   Hierarchical clustering (hclust + cutreeDynamic)
#   8.   Module eigengenes (aligned to mean expression per dataset)
#   9.   Cross-dataset module merging (pmin/pmax of |ME cor|)
#   10.  Module retention
#
# Output: module assignments, per-dataset eigengenes, PC weights,
#         per-dataset filters, rotation, var_explained, PDF report.
#
# Compatible with: nb_transfer, nb_plot_dendro (via nb_summary fields)
###############################################################################

#' Multi-dataset consensus network analysis
#'
#' Performs a 10-step consensus network analysis across multiple datasets:
#' (1-2) Filtered adjacency per dataset, (3-4) Direction filter + consensus,
#' (5) Soft thresholding, (6) TOM, (7) Hierarchical clustering,
#' (8) Module eigengenes, (9) Cross-dataset module merging,
#' (10) Module retention.
#'
#' @param datan_list List of data frames or matrices (rows = samples,
#'   cols = features). All datasets must share the same column names.
#' @param filter_method Filtering method: \code{"spearman"}, \code{"pearson"},
#'   \code{"kendall"}, \code{"boosting"}, or \code{"skip"}.
#' @param soft_power Numeric soft-thresholding power, or \code{NULL} for none.
#' @param consensus_method How to combine adjacencies across datasets:
#'   \code{"min"} or \code{"max"}.
#' @param network_type Network type (currently only \code{"unsigned"}).
#' @param filter_dir Logical. Apply direction filter (require consistent
#'   correlation sign across datasets)?
#' @param min_cluster_size Integer. Minimum module size.
#' @param ME_diss_thres Numeric (0, 1). Module eigengene dissimilarity
#'   threshold for merging.
#' @param merge_criterion How to evaluate cross-dataset ME correlation for
#'   merging: \code{"all"} (pmin) or \code{"any"} (pmax).
#' @param module_retention Criterion for retaining modules: \code{"any"} or
#'   \code{"all"}.
#' @param deepSplit Integer passed to \code{\link[dynamicTreeCut]{cutreeDynamic}}.
#' @param method Correlation method for adjacency calculation.
#' @param scale Logical. Scale and center data before PCA?
#' @param stepno Integer. Number of boosting steps for \code{\link{nb_filter}}.
#' @param n_pc Integer. Number of principal components.
#' @param robust_PCs Logical. Use rank-based (Spearman) PCA?
#' @param nb_min_varExpl Numeric. Minimum variance explained for module
#'   eigengenes.
#' @param cores Integer. Number of CPU cores.
#' @param pdf_report Path for PDF report, or \code{NULL} to skip.
#' @param mask_cache_dir Directory for caching filter masks, or \code{NULL}
#'   for \code{tempdir()}.
#' @param verbose Logical. Print progress messages?
#'
#' @return A list with components:
#'   \describe{
#'     \item{names}{Feature names.}
#'     \item{colors}{Final module color assignments (WGCNA color labels).}
#'     \item{MEs}{Module eigengenes for the first dataset.}
#'     \item{var_explained}{Variance explained per dataset.}
#'     \item{rotation}{Mean PC rotation matrix across datasets.}
#'     \item{filter}{Per-dataset filter edges.}
#'     \item{colors_initial}{Pre-merge module assignments.}
#'     \item{colors_merged}{Post-merge, pre-retention assignments.}
#'     \item{modules_retained}{Retained module color labels.}
#'     \item{geneTree}{Hierarchical clustering dendrogram.}
#'     \item{ME_all_datasets}{Module eigengenes per dataset.}
#'     \item{ME_cor_per_dataset}{ME correlation matrices per dataset.}
#'     \item{ME_tree}{ME dendrogram.}
#'     \item{pc1_weights}{PC1 weight comparison across datasets.}
#'     \item{filter_stats}{Edge filtering statistics.}
#'     \item{params}{Analysis parameters.}
#'     \item{n_datasets}{Number of datasets.}
#'     \item{n_features}{Number of features.}
#'     \item{n_modules_initial}{Modules before merging.}
#'     \item{n_modules_merged}{Modules after merging.}
#'     \item{n_modules_final}{Modules after retention.}
#'   }
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18', package='netboost')
#' # Create two mock datasets from the example data
#' d1 <- tcga_aml_meth_rna_chr18[1:40, 1:50]
#' d2 <- tcga_aml_meth_rna_chr18[41:80, 1:50]
#' res <- nb_consensus(datan_list = list(d1, d2),
#'     filter_method = "spearman", stepno = 20L,
#'     min_cluster_size = 5L, ME_diss_thres = 0.25,
#'     verbose = FALSE)
#'
#' @export
nb_consensus <- function(
    datan_list,
    filter_method    = c("spearman", "pearson", "kendall", "boosting", "skip"),
    soft_power       = NULL,
    consensus_method = c("min", "max"),
    network_type     = "unsigned",
    filter_dir       = TRUE,
    min_cluster_size = 15L,
    ME_diss_thres    = 0.25,
    merge_criterion  = c("all", "any"),
    module_retention = c("any", "all"),
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
    network_type     <- "unsigned"  # always unsigned
    method       <- match.arg(method)
    module_retention <- match.arg(module_retention)
    merge_criterion  <- match.arg(merge_criterion)

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
            stop(sprintf("Dataset %d has only %d samples (need >= 3).",
                         i, nrow(d)))
        if (ncol(d) < 2)
            stop(sprintf("Dataset %d has only %d features (need >= 2).",
                         i, ncol(d)))
        if (is.null(colnames(d)))
            stop(sprintf("Dataset %d has no column names.", i))
        if (is.data.frame(d)) {
            non_num <- !vapply(d, is.numeric, logical(1))
            if (any(non_num))
                stop(sprintf("Dataset %d has non-numeric columns: %s",
                             i, paste(names(non_num)[non_num], collapse = ", ")))
        }
    }

    p      <- ncol(datan_list[[1]])
    fnames <- colnames(datan_list[[1]])
    for (i in 2:n_ds) {
        if (ncol(datan_list[[i]]) != p)
            stop(sprintf("Dataset %d has %d features but dataset 1 has %d.",
                         i, ncol(datan_list[[i]]), p))
        if (!all(colnames(datan_list[[i]]) == fnames))
            stop(sprintf(
                "Column names in dataset %d do not match dataset 1.", i))
    }

    if (!is.null(soft_power)) {
        if (!is.numeric(soft_power) || soft_power < 0)
            stop("soft_power must be a non-negative number or NULL.")
    }
    if (!is.numeric(ME_diss_thres) || ME_diss_thres <= 0 || ME_diss_thres >= 1)
        stop("ME_diss_thres must be between 0 and 1 (exclusive).")
    if (!is.numeric(min_cluster_size) || min_cluster_size < 2)
        stop("min_cluster_size must be >= 2.")

    cores <- as.integer(cores)
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
    # ==================================================================
    cor_use <- "everything"
    for (i in seq_len(n_ds)) {
        if (anyNA(datan_list[[i]])) {
            cor_use <- "pairwise.complete.obs"
            break
        }
    }
    if (cor_use == "everything")
        msg("  No NAs - using fast BLAS correlation path")

    # --- Filter mask cache ---
    ds_sig <- paste0(vapply(datan_list, function(d)
        paste(nrow(d), ncol(d), sep = "x"), character(1)), collapse = "_")
    cache_key <- sprintf("mask_%s_p%d_%s_step%d",
                          ds_sig, p, filter_method, stepno)
    cache_dir <- if (!is.null(mask_cache_dir)) mask_cache_dir else tempdir()
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    cached <- NULL
    if (file.exists(cache_file)) {
        msg(sprintf("  Loading cached filter masks: %s", cache_file))
        cached <- readRDS(cache_file)
    }

    msg("\nStep 1-2: Computing filtered adjacencies...")
    adj_raw <- list()
    filter_masks <- list()
    filter_per_ds <- list()
    active_genes <- seq_len(p)

    for (i in seq_len(n_ds)) {
        n_i <- nrow(datan_list[[i]])

        # Progressive filtering: for min consensus, skip inactive genes
        use_subset <- (consensus_method == "min" && i > 1 &&
                       length(active_genes) < p)
        if (use_subset) {
            p_i <- length(active_genes)
            dat_i <- datan_list[[i]][, active_genes, drop = FALSE]
            msg(sprintf(
                "  DS%d: progressive - %d/%d active genes (%.0f%% saved)",
                i, p_i, p, 100 * (1 - (p_i / p)^2)))
        } else {
            p_i <- p
            dat_i <- datan_list[[i]]
        }

        if (filter_method != "skip") {
            if (!is.null(cached) && length(cached) >= i) {
                msg(sprintf("  DS%d: using cached %s mask",
                            i, filter_method))
                mask_full <- cached[[i]]
                filter_per_ds[[i]] <- which(
                    mask_full & upper.tri(mask_full), arr.ind = TRUE)
                if (use_subset) {
                    mask <- mask_full[active_genes, active_genes]
                } else {
                    mask <- mask_full
                }
                rm(mask_full)
            } else {
                msg(sprintf("  DS%d: nb_filter (method=%s, stepno=%d)...",
                            i, filter_method, stepno))
                nb_edges <- nb_filter(
                    as.data.frame(dat_i),
                    filter_method = filter_method,
                    stepno = stepno,
                    until  = 0L,
                    mode   = 2L,
                    cores  = cores,
                    verbose = verbose
                )
                filter_per_ds[[i]] <- nb_edges
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
            msg(sprintf("  DS%d: %s (no filter)...", i, method))
            adj_sub <- WGCNA::cor(dat_i, method = method, use = cor_use)
            diag(adj_sub) <- 0
            filter_per_ds[[i]] <- matrix(integer(0), ncol = 2,
                dimnames = list(NULL, c("row", "col")))
        }

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
        fstats$n_after_filter[i] <- sum(
            adj_raw[[i]][upper.tri(adj_raw[[i]])] != 0)
        msg(sprintf("  DS%d: %d edges after filter",
                    i, fstats$n_after_filter[i]))

        # Update active genes for progressive filtering
        if (consensus_method == "min") {
            has_edge <- logical(p)
            for (s in seq(1, p, BS)) {
                e <- min(s + BS - 1, p)
                has_edge[s:e] <- rowSums(
                    adj_raw[[i]][s:e, , drop = FALSE] != 0) > 0
            }
            active_genes <- intersect(active_genes, which(has_edge))
            msg(sprintf("  Active genes after DS%d: %d/%d",
                        i, length(active_genes), p))
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
        in_2   <- n_nz >= 2

        if (filter_dir) {
            ss <- Reduce("+", lapply(blocks, sign))
            same_dir <- abs(ss) == n_nz
            n_dir_filt <- n_dir_filt + sum(in_2 & !same_dir)
            valid <- in_2 & same_dir
        } else {
            valid <- in_2
        }

        if (consensus_method == "min") {
            ab <- lapply(blocks, function(b) {
                a <- abs(b); a[a == 0] <- Inf; a
            })
            ca <- Reduce(pmin, ab)
            ca[is.infinite(ca)] <- 0
        } else {
            ca <- Reduce(pmax, lapply(blocks, abs))
        }

        adj_consensus[s:e, ] <- ca * valid
    }

    diag(adj_consensus) <- 0
    rownames(adj_consensus) <- fnames
    colnames(adj_consensus) <- fnames
    fstats$n_dir_filtered <- as.integer(n_dir_filt / 2)
    fstats$n_consensus    <- sum(adj_consensus[upper.tri(adj_consensus)] != 0)

    # Collect example edges for scatter plots (before freeing adj_raw)
    example_edges <- .sample_example_edges(
        adj_raw, adj_consensus, fnames, n_examples = 6)

    rm(adj_raw); gc()

    msg(sprintf("  Direction-filtered: %d", fstats$n_dir_filtered))
    msg(sprintf("  Consensus edges:    %d", fstats$n_consensus))

    if (fstats$n_consensus == 0)
        stop(paste("No consensus edges.",
                   "Try a different filter_method or disable direction filter."))

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
    # STEP 6: TOM
    # ==================================================================
    msg("\nStep 6: TOM...")
    TOM <- WGCNA::TOMsimilarity(abs(adj_consensus),
                                 TOMType = "unsigned", verbose = 0)
    dissTOM <- 1 - TOM
    rm(TOM, adj_consensus); gc()
    rownames(dissTOM) <- fnames
    colnames(dissTOM) <- fnames

    # ==================================================================
    # STEP 7: Clustering
    # ==================================================================
    msg("\nStep 7: Clustering...")
    geneTree    <- hclust(as.dist(dissTOM), method = "average")
    dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                  deepSplit = deepSplit,
                                  pamRespectsDendro = FALSE,
                                  minClusterSize = min_cluster_size)
    names(dynamicMods) <- fnames
    rm(dissTOM); gc()

    n_initial <- length(unique(dynamicMods[dynamicMods > 0]))
    msg(sprintf("  Initial modules: %d", n_initial))
    if (n_initial > 0) {
        tab <- sort(table(dynamicMods[dynamicMods > 0]), decreasing = TRUE)
        msg(sprintf("  Sizes: %s", paste(tab, collapse = ", ")))
    }
    if (n_initial == 0) {
        message("No modules found.")
        return(invisible(list(
            colors = setNames(WGCNA::labels2colors(dynamicMods), fnames),
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

    # Align MEs to mean expression of module genes
    msg("  Aligning MEs to mean expression...")
    mods_align <- sort(unique(dynamicMods[dynamicMods > 0]))
    for (ds in seq_len(n_ds)) {
        dat <- as.matrix(datan_list[[ds]])
        for (mod in mods_align) {
            me_col <- paste0("ME", mod)
            if (!(me_col %in% colnames(ME_per_ds[[ds]]))) next
            genes <- names(dynamicMods[dynamicMods == mod])
            mean_expr <- rowMeans(dat[, genes, drop = FALSE], na.rm = TRUE)
            if (cor(ME_per_ds[[ds]][, me_col], mean_expr,
                    use = "complete.obs") < 0)
                ME_per_ds[[ds]][, me_col] <- -ME_per_ds[[ds]][, me_col]
        }
    }

    # ==================================================================
    # STEP 9: Cross-dataset module merging
    # ==================================================================
    me_cols <- setdiff(colnames(ME_per_ds[[1]]), "ME0")
    ME_cor_ds <- NULL; ME_tree <- NULL

    if (length(me_cols) > 1) {
        msg(sprintf("\nStep 9: Merging (thres=%.2f, criterion=%s)...",
                    ME_diss_thres, merge_criterion))

        ME_cor_ds <- lapply(ME_per_ds, function(MEs)
            abs(WGCNA::cor(MEs[, me_cols, drop = FALSE],
                           use = "pairwise.complete.obs")))

        ME_cor_cross <- if (merge_criterion == "all")
            Reduce(pmin, ME_cor_ds) else Reduce(pmax, ME_cor_ds)
        diag(ME_cor_cross) <- 1

        ME_tree  <- hclust(as.dist(1 - ME_cor_cross), method = "average")
        ME_clust <- cutree(ME_tree, h = ME_diss_thres)

        mod_nums <- as.integer(sub("^ME", "", me_cols))
        names(ME_clust) <- mod_nums

        merged_colors <- dynamicMods
        for (cl in unique(ME_clust)) {
            mods <- mod_nums[ME_clust == cl]
            if (length(mods) > 1) {
                rep_mod <- min(mods)
                for (m in mods) merged_colors[dynamicMods == m] <- rep_mod
                msg(sprintf("  Merged %s -> %d",
                            paste(mods, collapse = ", "), rep_mod))
            }
        }
    } else {
        msg("\nStep 9: Only 1 module, no merge.")
        merged_colors <- dynamicMods
    }

    n_merged <- length(unique(merged_colors[merged_colors > 0]))
    msg(sprintf("  After merging: %d modules", n_merged))

    # ==================================================================
    # STEP 10: Module retention (size + cross-dataset PC1 cor > 0.8)
    # ==================================================================
    msg(sprintf("\nStep 10: Retention (%s, PC1 cor > 0.8)...",
                module_retention))
    all_mods <- sort(unique(merged_colors[merged_colors > 0]))
    retain   <- integer(0)
    for (mod in all_mods) {
        sz <- sum(merged_colors == mod)
        if (sz < min_cluster_size) {
            msg(sprintf("  Dropping module %d (size %d)", mod, sz))
            next
        }

        # Cross-dataset PC1 loading correlation
        genes <- names(merged_colors[merged_colors == mod])
        pc1_list <- list()
        for (ds in seq_len(n_ds)) {
            ex <- datan_list[[ds]][, genes, drop = FALSE]
            if (ncol(ex) >= 2) {
                pc <- prcomp(ex, center = TRUE, scale. = scale)
                w  <- pc$rotation[, 1]
                if (cor(pc$x[, 1], rowMeans(ex, na.rm = TRUE),
                        use = "complete.obs") < 0)
                    w <- -w
                pc1_list[[ds]] <- w
            }
        }

        cors <- c()
        for (a in seq_len(n_ds - 1)) for (b in (a + 1):n_ds) {
            if (!is.null(pc1_list[[a]]) && !is.null(pc1_list[[b]])) {
                v <- !is.na(pc1_list[[a]]) & !is.na(pc1_list[[b]])
                if (sum(v) >= 2)
                    cors <- c(cors, cor(pc1_list[[a]][v], pc1_list[[b]][v]))
            }
        }
        me_cor <- if (length(cors) > 0) mean(cors) else NA_real_

        if (!is.na(me_cor) && me_cor >= 0.8) {
            retain <- c(retain, mod)
            msg(sprintf("  Retaining module %d (size %d, PC1 cor = %.3f)",
                        mod, sz, me_cor))
        } else {
            msg(sprintf("  Dropping module %d (size %d, PC1 cor = %.3f)",
                        mod, sz, if (is.na(me_cor)) NA else me_cor))
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
    for (ds in seq_len(n_ds)) {
        dat <- as.matrix(datan_list[[ds]])
        for (mod in sort(unique(final_colors[final_colors > 0]))) {
            me_col <- paste0("ME", mod)
            if (!(me_col %in% colnames(ME_final[[ds]]))) next
            genes <- names(final_colors[final_colors == mod])
            mean_expr <- rowMeans(dat[, genes, drop = FALSE], na.rm = TRUE)
            if (cor(ME_final[[ds]][, me_col], mean_expr,
                    use = "complete.obs") < 0)
                ME_final[[ds]][, me_col] <- -ME_final[[ds]][, me_col]
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
        w_list <- list()
        rot_list_ds <- list()
        for (i in seq_len(n_ds)) {
            ex <- datan_list[[i]][, genes, drop = FALSE]
            if (ncol(ex) >= 2) {
                pc <- prcomp(ex, center = TRUE, scale. = scale)
                n_avail <- min(n_pc, ncol(pc$rotation))
                w  <- pc$rotation[, 1]; s <- pc$x[, 1]
                if (cor(s, rowMeans(ex, na.rm = TRUE),
                        use = "complete.obs") < 0)
                    w <- -w
                rot_ds <- pc$rotation[, seq_len(n_avail), drop = FALSE]
                if (cor(pc$x[, 1], rowMeans(ex, na.rm = TRUE),
                        use = "complete.obs") < 0)
                    rot_ds[, 1] <- -rot_ds[, 1]
            } else {
                w <- setNames(rep(NA_real_, length(genes)), genes)
                rot_ds <- matrix(NA_real_, length(genes), n_pc,
                    dimnames = list(genes, paste0("PC", seq_len(n_pc))))
            }
            w_list[[i]] <- w
            rot_list_ds[[i]] <- rot_ds
        }

        df <- data.frame(gene_id = genes, module = mod,
                         stringsAsFactors = FALSE)
        for (i in seq_len(n_ds))
            df[[paste0("PC1_dataset", i)]] <- w_list[[i]][genes]
        pc1c <- grep("^PC1_dataset", names(df))
        df$PC1_mean <- rowMeans(df[, pc1c, drop = FALSE], na.rm = TRUE)

        cors <- c()
        for (a in seq_len(n_ds - 1)) for (b in (a + 1):n_ds) {
            v <- !is.na(w_list[[a]]) & !is.na(w_list[[b]])
            if (sum(v) >= 2)
                cors <- c(cors, cor(w_list[[a]][v], w_list[[b]][v]))
        }
        df$PC1_cor <- if (length(cors) > 0) mean(cors) else NA_real_
        pc1_all[[as.character(mod)]] <- df

        # Mean rotation across datasets for this module
        mean_rot <- Reduce("+", rot_list_ds) / length(rot_list_ds)
        rownames(mean_rot) <- genes
        rotation_list[[as.character(mod)]] <- mean_rot
    }
    if (length(pc1_all) > 0) {
        all_weights <- do.call(rbind, pc1_all)
        rownames(all_weights) <- NULL
    } else {
        all_weights <- data.frame(gene_id = character(0),
                                   module = integer(0),
                                   PC1_mean = numeric(0),
                                   PC1_cor = numeric(0))
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
    # Convert numeric module IDs to WGCNA color labels
    # ==================================================================
    max_mod <- max(c(dynamicMods, merged_colors, final_colors), na.rm = TRUE)
    num2color <- setNames(WGCNA::labels2colors(0:max_mod),
                          as.character(0:max_mod))

    final_colors  <- setNames(num2color[as.character(final_colors)], fnames)
    dynamicMods   <- setNames(num2color[as.character(dynamicMods)], fnames)
    merged_colors <- setNames(num2color[as.character(merged_colors)], fnames)
    retain        <- unname(num2color[as.character(retain)])

    # Rename ME columns: ME1 -> MEturquoise, etc.
    .rename_ME <- function(nms) {
        vapply(nms, function(nm) {
            paste0("ME", num2color[sub("^ME", "", nm)])
        }, character(1), USE.NAMES = FALSE)
    }

    for (ds in seq_len(n_ds))
        colnames(ME_final[[ds]]) <- .rename_ME(colnames(ME_final[[ds]]))

    if (!is.null(ME_cor_ds)) {
        for (ds in seq_along(ME_cor_ds)) {
            rownames(ME_cor_ds[[ds]]) <- .rename_ME(
                rownames(ME_cor_ds[[ds]]))
            colnames(ME_cor_ds[[ds]]) <- .rename_ME(
                colnames(ME_cor_ds[[ds]]))
        }
    }

    if (!is.null(ME_tree))
        ME_tree$labels <- .rename_ME(ME_tree$labels)

    if (nrow(all_weights) > 0)
        all_weights$module <- unname(
            num2color[as.character(all_weights$module)])

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
            module_retention = module_retention, deepSplit = deepSplit,
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
                sum(final_colors != "grey"),
                100 * sum(final_colors != "grey") / p))
    msg(sprintf("  Unassigned:  %d", sum(final_colors == "grey")))
    msg(sprintf("  Modules: %d -> %d -> %d",
                n_initial, n_merged, n_final))
    if (n_final > 0) {
        sz <- table(final_colors[final_colors != "grey"])
        msg(sprintf("  Sizes: %d - %d (mean %.1f)",
                    min(sz), max(sz), mean(sz)))
    }
    if (nrow(all_weights) > 0 && any(!is.na(all_weights$PC1_cor))) {
        pc1_per_mod <- unique(all_weights[, c("module", "PC1_cor")])
        n_above <- sum(pc1_per_mod$PC1_cor >= 0.8, na.rm = TRUE)
        n_total_mod <- nrow(pc1_per_mod)
        msg(sprintf("  PC1 cor: mean=%.3f [%.3f, %.3f]",
                    mean(pc1_per_mod$PC1_cor, na.rm = TRUE),
                    min(pc1_per_mod$PC1_cor, na.rm = TRUE),
                    max(pc1_per_mod$PC1_cor, na.rm = TRUE)))
        msg(sprintf("  PC1 concordance >= 0.8: %d / %d modules",
                    n_above, n_total_mod))
    }
    msg(paste(rep("=", 60), collapse = ""))
    invisible(results)
}

# =============================================================================
# TRANSFER FUNCTION
# =============================================================================

#' Transfer consensus module assignments to new data
#'
#' Projects consensus module eigengenes onto new samples using the rotation
#' matrix from \code{\link{nb_consensus}}.
#'
#' @param nb_consensus_result Result from \code{\link{nb_consensus}}.
#' @param new_data Data frame or matrix (rows = samples, cols = features).
#' @param scale Logical. Scale and center \code{new_data} before projection?
#' @param robust_PCs Logical. Use rank-based PCA (Spearman)?
#' @return Data frame of projected module eigengenes (samples x modules).
#'
#' @examples
#' data('tcga_aml_meth_rna_chr18', package='netboost')
#' d1 <- tcga_aml_meth_rna_chr18[1:40, 1:50]
#' d2 <- tcga_aml_meth_rna_chr18[41:80, 1:50]
#' res <- nb_consensus(datan_list = list(d1, d2),
#'     filter_method = "spearman", stepno = 20L,
#'     min_cluster_size = 5L, ME_diss_thres = 0.25,
#'     verbose = FALSE)
#' new_dat <- tcga_aml_meth_rna_chr18[1:10, 1:50]
#' MEs_new <- nb_consensus_transfer(res, new_dat)
#'
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
        stop(paste("No overlapping feature names between",
                   "new_data and consensus result."))
    if (length(common) < length(fnames))
        message(sprintf("  nb_consensus_transfer: %d / %d features matched",
                        length(common), length(fnames)))

    dat <- new_data[, common, drop = FALSE]
    if (scale) dat <- scale(dat)
    if (robust_PCs) dat <- apply(dat, 2, rank)

    mods <- sort(unique(colors[colors != "grey"]))
    MEs  <- matrix(NA_real_, nrow = nrow(dat), ncol = length(mods),
                    dimnames = list(rownames(dat), paste0("ME", mods)))

    for (j in seq_along(mods)) {
        mod   <- mods[j]
        genes <- names(colors[colors == mod])
        genes <- intersect(genes, common)
        if (length(genes) < 2) next

        # Project using consensus rotation (PC1)
        w <- rot[genes, 1]
        w <- w[!is.na(w)]
        genes_use <- intersect(names(w), colnames(dat))
        if (length(genes_use) < 2) next

        scores <- dat[, genes_use, drop = FALSE] %*% w[genes_use]
        # Align sign: positive correlation with module mean expression
        mu <- rowMeans(dat[, genes_use, drop = FALSE], na.rm = TRUE)
        if (cor(scores, mu, use = "complete.obs") < 0) scores <- -scores
        MEs[, j] <- scores
    }

    as.data.frame(MEs)
}

# =============================================================================
# Helper: sample example edges for report scatter plots
# =============================================================================

.sample_example_edges <- function(adj_raw, adj_consensus, fnames,
                                   n_examples = 6) {
    p <- length(fnames)
    ut <- which(upper.tri(adj_consensus) & adj_consensus != 0, arr.ind = TRUE)
    if (nrow(ut) == 0) return(NULL)

    vals <- abs(adj_consensus[ut])
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
    bp <- barplot(edge_vals, names.arg = edge_labs,
                  col = c(rep("#4575b4", n_ds), "#91cf60"),
                  main = "Edge Counts", ylab = "Edges", las = 2,
                  cex.names = 0.8, ylim = c(0, ymax))
    graphics::text(bp, edge_vals,
                   labels = format(edge_vals, big.mark = ","),
                   pos = 3, cex = 0.7)

    # 1b: Edge density per dataset
    pct <- fs$n_after_filter / fs$n_total * 100
    ymax_pct <- max(pct) * 1.25
    bp2 <- barplot(pct, names.arg = paste("DS", seq_len(n_ds)),
            col = c("#2166ac", "#b2182b", "#1b7837", "#762a83")[seq_len(n_ds)],
            main = "Edge Density per Dataset",
            ylab = sprintf("%%  of total possible (%s)",
                           format(fs$n_total, big.mark = ",")),
            ylim = c(0, ymax_pct), cex.names = 0.8)
    graphics::text(bp2, pct, labels = sprintf("%.1f%%", pct),
                   pos = 3, cex = 0.8)

    # 1c: Edge fate pie
    cons   <- fs$n_consensus
    dir_f  <- fs$n_dir_filtered
    no_sig <- fs$n_total - max(fs$n_after_filter)
    other  <- fs$n_total - cons - dir_f - no_sig
    vals   <- c(cons, dir_f, max(other, 0), no_sig)
    labs   <- c(sprintf("Consensus\n(%s)", format(cons, big.mark = ",")),
                sprintf("Dir-filtered\n(%s)", format(dir_f, big.mark = ",")),
                sprintf("In 1 DS only\n(%s)",
                        format(max(other, 0), big.mark = ",")),
                sprintf("Not signif.\n(%s)", format(no_sig, big.mark = ",")))
    par(mar = c(4, 2, 4, 2))
    pie(vals, labels = labs,
        col = c("#91cf60", "#fc8d59", "#fee08b", "#d9d9d9"),
        main = "Edge Fate", cex = 0.8)

    # 1d: Parameters text
    par(mar = c(2, 2, 3, 2))
    plot.new()
    graphics::title("Parameters", line = 0.5)
    txt <- c(
        sprintf("Filter: %s (nb_filter, stepno=%d)",
                pa$filter_method, pa$stepno),
        sprintf("Network type: %s", pa$network_type),
        sprintf("Adjacency: %s", pa$method),
        sprintf("Direction filter: %s", pa$filter_dir),
        sprintf("Consensus: %s", pa$consensus_method),
        sprintf("Soft power: %s",
                ifelse(is.null(pa$soft_power), "none", pa$soft_power)),
        sprintf("ME merge thres: %.2f (cor >= %.2f, %s)",
                pa$ME_diss_thres, 1 - pa$ME_diss_thres, pa$merge_criterion),
        sprintf("Retention: %s | dS: %d | minSz: %d",
                pa$module_retention, pa$deepSplit, pa$min_cluster_size),
        "",
        sprintf("Features: %d", results$n_features),
        sprintf("In modules: %d (%.1f%%)",
                sum(fc != "grey"), 100 * sum(fc != "grey") / results$n_features),
        sprintf("Modules: %d -> %d -> %d",
                results$n_modules_initial, results$n_modules_merged,
                results$n_modules_final),
        {
            pw <- results$pc1_weights
            if (!is.null(pw) && nrow(pw) > 0 &&
                any(!is.na(pw$PC1_cor))) {
                pc1_mod <- unique(pw[, c("module", "PC1_cor")])
                n_ab <- sum(pc1_mod$PC1_cor >= 0.8, na.rm = TRUE)
                sprintf("PC1 cor > 0.8: %d / %d modules",
                        n_ab, nrow(pc1_mod))
            } else ""
        }
    )
    for (k in seq_along(txt))
        graphics::text(0.05, 0.95 - k * 0.08, txt[k], adj = 0,
                       cex = 0.85, family = "mono")

    # ---- Page 2: ME dendrogram ----
    if (!is.null(results$ME_tree)) {
        par(mfrow = c(1, 1), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
        tryCatch({
            plot(results$ME_tree,
                 main = "Module Eigengene Dendrogram (pre-merge)",
                 xlab = "", ylab = "Dissimilarity (1 - min|cor|)", sub = "",
                 cex = 0.8)
            abline(h = pa$ME_diss_thres, col = "red", lty = 2, lwd = 2)
            graphics::mtext(
                sprintf("merge threshold = %.2f  (cor = %.2f)",
                        pa$ME_diss_thres, 1 - pa$ME_diss_thres),
                side = 3, line = -1.5, col = "red", cex = 0.9)
        }, error = function(e) {
            plot.new()
            graphics::text(0.5, 0.5,
                sprintf("Dendrogram not available\n(%d modules)",
                        results$n_modules_initial), cex = 1.2)
        })
    }

    # ---- Page 3: ME correlation heatmaps per dataset ----
    if (!is.null(results$ME_cor_per_dataset) &&
        length(results$ME_cor_per_dataset) > 0) {
        par(mfrow = c(1, n_ds), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
        for (ds in seq_len(n_ds)) {
            mc <- results$ME_cor_per_dataset[[ds]]
            lbl <- sub("^ME", "", rownames(mc))
            graphics::image(seq_len(nrow(mc)), seq_len(ncol(mc)), mc,
                  col = colorRampPalette(c("white", "orange", "red"))(50),
                  zlim = c(0, 1), axes = FALSE,
                  main = sprintf("DS%d: |ME cor|", ds))
            graphics::axis(1, seq_len(nrow(mc)), lbl,
                           las = 2, cex.axis = 0.5)
            graphics::axis(2, seq_len(ncol(mc)), lbl,
                           las = 2, cex.axis = 0.5)
        }
    }

    # ---- Page 4: Module sizes ----
    par(mfrow = c(1, 1), mar = c(6, 5, 4, 3), oma = c(0, 0, 0, 0))
    tab <- sort(table(fc[fc != "grey"]), decreasing = TRUE)
    if (length(tab) > 0) {
        ymax_tab <- max(tab) * 1.1
        bp <- barplot(tab, main = "Final Module Sizes", xlab = "Module",
                      ylab = "Genes", col = "steelblue", las = 2,
                      ylim = c(0, ymax_tab), cex.names = 0.8)
        graphics::text(bp, tab, labels = tab, pos = 3, cex = 0.65)
    } else {
        plot.new()
        graphics::text(0.5, 0.5, "No modules retained", cex = 1.5)
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
        graphics::mtext("Top Consensus Edges: Gene-Gene Scatters",
                        outer = TRUE, side = 3, line = 0.3, cex = 1.0)
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
        graphics::mtext("PC1 Weights Across Datasets", outer = TRUE,
                        side = 3, line = 0.3, cex = 1.0)
    }

    message(sprintf("  Report: %s (%d pages)", pdf_path,
                    2 + !is.null(results$ME_cor_per_dataset) +
                    1 + !is.null(results$example_edges) + (n_ds == 2)))
}
