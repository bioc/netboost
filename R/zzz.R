## Note: here also the setup NAMESPACE directives are set (as Roxygen is
## creating the namespace, those are given as Roxygen attributes). The
## importFrom statements are required to instantly load the linked libraries
## Rcpp and RcppParallel, else loading of the own shared lib would fail (the
## imported functions do not matter, but for each package at least one import
## must be present).

## Those are required for CRAN checks but bug out on BiocCheck
#' Package startup: used to fetch installation path of the own package,
#' as required for executing binary programs delivered with it.
#' 
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel setThreadOptions
#' @importFrom parallel mclapply
#' @importFrom colorspace rainbow_hcl
#' @importFrom grDevices dev.off gray pdf colorRampPalette adjustcolor
#' @importFrom graphics abline layout par plot barplot text pie image axis mtext
#'   title plot.new
#' @importFrom stats as.dendrogram as.dist cor cov prcomp hclust
#'   order.dendrogram pt cutree lm setNames
#' @importFrom dynamicTreeCut cutreeDynamic indentSpaces printFlush
#' @importFrom impute impute.knn
#' @importFrom WGCNA allowWGCNAThreads mergeCloseModules plotDendroAndColors
#' @importFrom WGCNA moduleColor.getMEprefix pickSoftThreshold TOMsimilarity
#'   labels2colors
#' @importFrom utils data packageDescription
#' @importFrom methods is
#'
#' @useDynLib netboost
#'
#' @examples 
#' \dontrun{nb_example()}
#' @return none
#' @param libname Path to R installation (base package dir)
#' @param pkgname Package name (should be "netboost")
.onAttach <- function(libname, pkgname) {
    desc <- packageDescription(pkgname)

    # If no default core count given, detect.  
    if (is.null(getOption("mc.cores")) || !is.integer(getOption("mc.cores"))) {
        # logical = FALSE is not working correctly if CPU has logical cores, which
        # are disabled (at least Linux).
        # Means: if CPU has logical cores, core count should be set manually.
        # cores <- parallel::detectCores() # Bioconductor does not like this
        cores <- NA

        if (is.na(cores)) cores <- 1

        options("mc.cores" = cores)
    }
    
    ## Optional startup message, mainly for development.
    packageStartupMessage(
        paste(pkgname,
              desc$Version,
              "loaded"),
        paste(
            "Default CPU cores:",
            getOption("mc.cores"),
            "\n"),
        appendLF = TRUE
    )
    #                              "Loaded from:", libname),

    # Create temp subfolder in tempdir()
    netboostTmpCleanup()
}

## #' If package detached, clean up temporary folders.
## #' @return none
## #' @param libpath Library path (unused)
##.onDetach <- function(libpath) {
##    print("kthnxbye")
##}
