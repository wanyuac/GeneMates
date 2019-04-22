#' @title Retrive allelic presence-absence information given a vector of allele names
#'
#' @description This function is to be run after the findPhysLink pipeline.
#'
#' @param targets A character vector of allele names whose co-occurrence information is to be retrieved
#' @param assoc The element "assoc" in the output of findPhysLink.
#' @param allelic.pam The element "assoc$A" in the output of findPhysLink.
#' @param min.count (>= 1) The minimal number of alleles that an isolate must carry to be included in the extracted information
#' @param phylo.dists A matrix of phylogenetic distances between tips of a tree. It can be generated using the function cophenetic.phylo(tree).
#' @param dist.method The parameter for the method argument of the base function dist
#' @param clust.method The parameter for the method argument of the base function hclust
#'
#' @examples
#' retrieveAlleleSetInfo(targets = c("aadA2", "dfrA12", "sul1"), assoc = a[["assoc"]], allelic.pam = a[["alleles"]][["A"]], min.count = 2, phylo.dists = cophenetic.phylo(tr))
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Dependency: ape
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  Last edit: 2 June 2017

retrieveAlleleSetInfo <- function(targets, assoc, allelic.pam, min.count = 2, phylo.dists = NULL,
                                  dist.method = "binary", clust.method = "single") {
    require(ape)

    # check argument
    if (min.count < 1) {
        print("Warning: min.count is reset to one as it cannot be smaller than one.")
        min.count <- 1
    }

    # extract association information of these alleles
    clr <- subset(assoc, (y %in% targets) & (x %in% targets))
    clr <- clr[order(clr$pair, decreasing = FALSE), ]  # order rows by pair IDs

    # extract allelic presence/absence status
    A <- allelic.pam[, targets]
    A <- A[as.integer(rowSums(A)) >= min.count, ]  # remove rows that are empty or only have a small number of alleles

    # cluster alleles on columns
    d <- dist(t(A), method = dist.method)
    hc <- hclust(d, method = clust.method)
    A <- A[, hc$order]

    # cluster alleles on rows
    d <- dist(A, method = dist.method)
    hc <- hclust(d, method = clust.method)
    A <- A[hc$order, ]

    # find out isolates that have all target alleles and extract their pairwise phylogenetic distances
    co <- rownames(A)[as.integer(rowSums(A)) == ncol(A)]
    pd <- NULL  # default value of phylo.dists
    if (length(co) == 0) {
        co <- NULL
    } else if (!is.null(phylo.dists)) {
        pd <- phylo.dists[co, co]  # a subset of the distance matrix for isolates having all target alleles
    }  # else, do nothing (pd = NULL, co != "")

    return(list(targets = targets, clr = clr, pam = A, co = co, pd = pd))
}
