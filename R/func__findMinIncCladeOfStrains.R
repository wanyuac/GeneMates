#' @title Finding out and summarising the minimal inclusive clade containing all given strains
#'
#' @description This function works in a similar way as the function findMinIncClade, except it works on a
#' given set of strains rather than a pair of alleles.
#'
#' @param strains a vector of strain names
#' @param clade.pam a presence/absence matrix of samples in each clade of an input tree. It can be obtained
#' using the function tree2Clades of phylix.
#' @param clade.sizes (optional) A named vector of integers for the number of samples in each clade. It can
#' be obtained from the element "sizes" in the outputs of the function tree2Clades. Optional.
#'
#' @examples
#' mic <- findMinIncCladeOfStrains(strains = c("strain1", "strain2"), clade.pam = assoc[["struc"]][["clades"]][["pam"]],
#' clade.sizes = assoc[["struc"]][["clades"]][["sizes"]])
#'
#' @author Yu Wan (\email{wanyuac@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition and the lastest edition: 27 Sep 2017

findMinIncCladeOfStrains <- function(strains = NULL, clade.pam = NULL, clade.sizes = NULL) {
    # check and prepare arguments
    if (is.null(strains) | is.null(clade.pam)) {
        stop("Argument error: strains or clades.pam must not be NULL.")
    }
    if (is.null(clade.sizes)) {  # make a named vector for clade sizes when it is not provided
        clade.sizes <- colSums(x = clade.pam)
    }

    # look for the clade in a similar way to the function .searchMinCladePerPair.
    strain.presence <- rownames(clade.pam) %in% strains
    n <- length(strains)
    candidate.clades <- names(clade.sizes)[which(as.integer(clade.sizes) >= n)]  # Each candiate clade must have at least n samples.
    clade.pam <- clade.pam[, candidate.clades]  # to reduces the number of comparisons
    is.inc <- as.logical(apply(clade.pam, 2, function(c) sum(as.logical(c) | strain.presence) == sum(c)))  # sum(c) always >= n; equality holds when input strains form a subset of a clade.
    inc.clades <- candidate.clades[is.inc]  # IDs of all inclusive clades
    min.clade <- inc.clades[which.min(as.integer(clade.sizes[inc.clades]))]  # The usage of which.min is correct for a bipartition tree.
    min.clade.size <- as.integer(clade.sizes[[min.clade]])

    return(list(clade = min.clade, size = min.clade.size, freq = round(n / min.clade.size * 100, digits = 4)))
}
