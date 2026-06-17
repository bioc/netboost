#' Returns the absolute path to "exec" folder in the package.
#'
#' @return Absolute path of installed package
netboostPackagePath <- function() {
    return(system.file(package="netboost"))
}

#' Returns the absolute path to temporary folder of the package.
#' To change temporary path, use normal R variables (TEMPDIR etc).
#'
#' @return Absolute path for temporary folder
netboostTmpPath <- function() {
    return(file.path(tempdir(), "netboost"))
}

#' Cleans the netboost temporary folder. This can be useful during the session
#' as mcupgma creates vast directory structures (for iterations).
#' Creates the own folder (all netboost temporary data is stored in
#' netboostTmpPath(), which is equal to tempdir()/netboost).
#' Also used for first time setup of folder.
#'
#' @param verbose Flag verbose
#' @return none
netboostTmpCleanup <- function(verbose = FALSE) {
    folder <- netboostTmpPath()

    if (dir.exists(folder)) {
        if (verbose)
            message(paste("Netboost: cleaning temporary folder:", folder))

        ## Delete and recreate more convenient than globbing through the folders
        unlink(folder, recursive = TRUE)
    }

    # Create or recreate folder.
    dir.create(folder)
}