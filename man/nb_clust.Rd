% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netboost.R
\name{nb_clust}
\alias{nb_clust}
\title{Netboost clustering step}
\usage{
nb_clust(
  filter,
  dist,
  datan,
  max_singleton = dim(datan)[2],
  min_cluster_size = 2L,
  ME_diss_thres = 0.25,
  cores = getOption("mc.cores", 2L),
  qc_plot = TRUE,
  n_pc = 1,
  robust_PCs = FALSE,
  nb_min_varExpl = 0.5,
  method = c("pearson", "kendall", "spearman")
)
}
\arguments{
\item{filter}{Filter-Matrix as generated by the nb_filter function.}

\item{dist}{Distance-Matrix as generated by the nb_dist function.}

\item{datan}{Data frame were rows correspond to samples and columns to
features.}

\item{max_singleton}{Integer. The maximal singleton in the clustering.
Usually equals the number of features.}

\item{min_cluster_size}{Integer. The minimum number of features in one 
module.}

\item{ME_diss_thres}{Numeric. Module Eigengene Dissimilarity Threshold for
merging close modules.}

\item{cores}{Integer. Amount of CPU cores used (<=1 : sequential)}

\item{qc_plot}{Logical. Create plot.}

\item{n_pc}{Number of principal components and variance explained
entries to be calculated. The number of returned variance explained entries
is currently ‘min(n_pc,10)’. If given ‘n_pc’ is greater than 10, a warning 
is issued.}

\item{robust_PCs}{Should PCA be calculated on ranked data (Spearman PCA)?
Rotations will not correspond to original data if this is applied.}

\item{nb_min_varExpl}{Minimum proportion of variance explained for
returned module eigengenes. The number of PCs is capped at n_pc.}

\item{method}{A character string specifying the method to be used for
correlation coefficients.}
}
\value{
List
}
\description{
Netboost clustering step
}
\examples{
 data('tcga_aml_meth_rna_chr18',  package='netboost')
 cores <- as.integer(getOption('mc.cores', 2))
 datan <- as.data.frame(scale(tcga_aml_meth_rna_chr18, center=TRUE, 
 scale=TRUE))
 filter <- nb_filter(datan=datan, stepno=20L, until=0L, progress=1000L,
 cores=cores,mode=2L)
 dist <- nb_dist(datan=datan, filter=filter, soft_power=3L, cores=cores)
 max_singleton = dim(tcga_aml_meth_rna_chr18)[2]
 pdf("test.pdf",width=30)
 sum_res <- nb_clust(filter=filter, dist=dist, datan=datan,
 max_singleton=max_singleton, min_cluster_size=10L, ME_diss_thres=0.25,
 cores=cores, qc_plot=TRUE, n_pc=2L, nb_min_varExpl=0.5)
 dev.off()
}
