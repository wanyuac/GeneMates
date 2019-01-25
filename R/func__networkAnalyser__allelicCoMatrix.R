#' @title Create a co-occurrence matrix from allelic PAM
#'
#' @param clusters A named or unnamed list of allele name vectors that define
#' alleles of each cluster. Every vector must contain at least two alleles.
#' @param pam An allelic presence-absence matrix, where 1 denotes presence of an
#' allele and 0 denotes absence of an allele.
#' @param hclust.method Method for clustering columns in the output co-occurrence
#' matrix. Default: NULL (no clustering)
#' @param hclust.dist Distance metric used for column-wise clustering. Default:
#' binary.
#'
#' @note hclust.method and hclust.dist are passed to R base functions hclust and
#' dist as arguments.
#'
#' @return A list of elements: m, a co-occurrence matrix; c, clusters, recordered
#' if hclust.method != NULL; g: a vector (method, dist) for arguments of
#' column-wise clustering (grouping).
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2018 Yu Wan <wanyuac@gmail.com>
#  Licensed under the Apache License, Version 2.0
#  First and the latest edition: 25 Jan 2019

allelicCoMatrix <- function(clusters = NULL, pam = NULL, hclust.method = NULL,
                            hclust.dist = "binary") {
    if (is.null(clusters) || is.null(pam)) {
        stop("Argument error: both clusters and pam must be provided.")
    }

    if (is.null(names(clusters))) {
        print("Arbitrarily assign cluster IDs.")
        names(clusters) <- paste0("C", as.character(1 : length(clusters)))  # C1, C2, ...
    }

    M <- matrix(0, nrow = nrow(pam), ncol = length(clusters),
                dimnames = list(rownames(pam), names(clusters)))  # a co-occurrence matrix
    for (c in names(clusters)) {
        pam_c <- pam[, clusters[[c]]]  # A subset of
        co <- as.logical(pam_c[, 1])
        for (i in 2 : ncol(pam_c)) {
            co <- co & pam_c[, i]
        }
        M[, c] <- as.integer(co)  # co-occurrence of all alleles in the current cluster
    }

    if (!is.null(hclust.method)) {
        print("Cluster columns of the co-occurrence network.")
        h <- hclust(d = dist(x = t(M), method = hclust.dist), method = hclust.method)
        M <- M[, h$order]
        clusters <- clusters[colnames(M)]
    }

    return(list(m = M, c = clusters, g = c(hclust.method, hclust.dist)))
}
