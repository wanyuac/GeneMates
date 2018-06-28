#' @title Estimating correlations between clades and projections of samples
#'
#' @description This function estimates point-biserial correlation coefficients between each clade and projections of
#' samples on each axis defined by eigenvectors of positive eigenvalues. The output links projections (based on singular-value
#' decomposition of the cgSNP matrix) to population structure. This function ignores the "root clade", which is the clade based
#' on the root when the tree is rooted, because it does not show any variation in the presence/absence status of samples (all
#' samples must be present in this super clade). P-values do not apply in this scenario because we are only describing the correlation
#' between observations but not going to make any predictions.
#'
#' @param clades A presence/absence matrix of samples in every clade. It can be generated with the function tree2Clades.
#' @param projections A matrix of projections of samples on eigenvectors of positive eigenvalues. Rows are sample names and columns
#' are axes. It can be obtained from the element C in the output of the function projectSamples.
#' @param clade.sizes (Optional) A vector of named integers for sizes of every clade in C. It can be obtained from the element
#' "sizes" in the output of the function tree2Clades. Providing this argument makes the function to run faster, although the
#' function calculates clade sizes when this vector is absent. It is recommended to have the order of clade names in this vector
#' identical to column names of the clades matrix, although the function makes a correction for the order when it mismatches.
#' @param n.cores An integer determining the number of cores used for parallel computing. Valid values are the same as those for
#' the findPhysLink function.
#'
#' Dependency: parallel, data.table
#'
#' @examples
#' C <- projectSamples(...)
#' clades <- tree2Clades(...)
#' r <- corCladeProj(clades = clades[["pam"]], projections = C[["C"]], clade.sizes = clades[["sizes"]], n.cores = 8)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 19 July 2017, the lastest edition: 1 Aug 2017

corCladeProj <- function(clades = NULL, projections = NULL, clade.sizes = NULL, n.cores = -1) {
    require(parallel)
    require(data.table)
    print(paste0(Sys.time(), ": estimating correlation coefficients between each clade and projections of samples on each axis."))

    # match rows of matrices clades and projections
    samples <- rownames(projections)
    if (any(rownames(clades) != samples)) {
        print("Warning: rearranging rows of the clade matrix to match samples in the projection matrix.")
        clades <- clades[samples, ]  # An error arises when two matrices involve different samples.
    }

    # calculate clade sizes when they are not provided
    clade.ids <- colnames(clades)
    if (class(clade.sizes) != "integer") {
        print("Calculating clade sizes.")
        clade.sizes <- apply(clades, 2, function(x) as.integer(sum(x)))
    } else {
        # match clade names in the size vector to those on the columns of the matrix clades
        if (any(names(clade.sizes) != clade.ids)) {
            print("Warning: rearranging clade sizes to have their names matched to column names of the clade matrix.")
            clade.sizes <- clade.sizes[clade.ids]
        }
    }

    # remove the "root clade"
    n <- nrow(projections)  # sample size
    axes <- colnames(projections)
    where.root <- which(as.integer(clade.sizes) == n)  # which column corresponds to the "root clade"
    if (length(where.root) > 0) {
        print(paste0("Removing the \"root clade\"", clade.ids[where.root], "from the clade matrix.", sep = " "))
        clades <- clades[, -where.root]  # remove the column corresponding to the root clade
        clade.ids <- colnames(clades)  # Do not need to update the vector clade.sizes as we will not use it anymore.
    }

    # estimate correlation coefficients
    cl <- makeCluster(.setCoreNum(n.cores = n.cores, cores.avai = detectCores()))  # cores.avai >= 1; n.cores >= 1
    clusterExport(cl = cl, varlist = list("clades", "projections", "axes"), envir = environment())
    print(paste0(Sys.time(), ": Starting estimating point-biserial correlation coefficients between each clade and projections on every axis."))
    r <- parLapply(cl, clade.ids, .calcCorCladeProj, clades, projections, axes)  # wait until all jobs finish
    stopCluster(cl)
    print(paste0(Sys.time(), ": all correlation coefficients have been estimated."))

    return(rbindlist(r))
}

# Estimating the correlation coefficient for a given clade with projections on every axis.
# Notice projections on different axes are linearly independent.
.calcCorCladeProj <- function(i, clades, projections, axes) {
    c <- clades[, i]  # take a specific clade for constructing the model

    # A point-biserial correlation is a Pearson product-moment correlation when one variable is dichotomous
    # and the other is continuous.
    r <- sapply(axes, function(j) round(cor(x = c, y = projections[, j], method = "pearson"), digits = 8))
    r <- data.frame(clade = rep(i, times = length(r)), axis = names(r), cor = as.numeric(r), stringsAsFactors = FALSE)
    r <- r[order(r$cor, decreasing = TRUE), ]

    return(r)
}
