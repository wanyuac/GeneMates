#' @title Compare similarity between edges' occurrence
#'
#' @description Compute a matrix of Simpson's similarity coefficients for the
#' presence of edges.
#'
#' @param nwk A Graph object created by the function mkNetwork for a single network.
#' It has three attributes: id, E and V. Required columns of E are "y" and "x",
#' which are vertex names.
#' @param apam Allelic presence-absence matrix. Rows: strain names; columns:
#' allele names.
#' @param directed A logical argument specifying whether nwk stores a directed network.
#' Default: TRUE.
#'
#' @return A list of three elements: nwk, an undirected network converted from
#' the input network; ep, a data frame of edge pairs; sim, the similarity matrix
#' for edges, which is a symmetric matrix.
#'
#' @examples compEdgeOccur
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First version: 6 Aug 2018; the latest edition: 11 Aug 2018

compEdgeOccur <- function(nwk, apam, directed = TRUE) {
    # Convert nwk into an undirected network
    if (directed) {  # Vertices will not change in this process.
        if (!any(names(nwk@E) == "pair")) {
            nwk@E <- assignPairID(lmms = nwk@E, from = 1, paired.rows = FALSE)  # A "pair" column is appended.
        }
        nwk@E <- nwk@E[!duplicated(nwk@E$pair), ]  # only preserve one edge per pair of edges
    }

    # Prepare the data frame of edge pairs
    ep <- nwk@E[, c("pair", "y", "x")]
    ep <- ep[order(ep$pair, decreasing = FALSE), ]
    ep$pair <- paste0("e", ep$pair)  # 1, 2, ... => e1, e2, ...
    names(ep)[1] <- "edge_ID"  # replace the first column name
    ep$edge_name <- mapply(function(a1, a2) paste(a1, a2, sep = "&"), ep$y, ep$x)
    eids <- ep$edge_ID
    n <- length(eids)  # number of undirected edges

    # Calculate a Simpson similarity for each pair of edges
    sim <- matrix(1.0, nrow = n, ncol = n, dimnames = list(eids, eids))
    for (i in 1 : (n - 1)) {
        ei <- ep[i, ]
        ei_id <- ei$edge_ID
        ei_x <- ei$x
        ei_y <- ei$y
        co_i <- as.logical(apam[, ei_x]) & as.logical(apam[, ei_y])  # co-occurrence of x, y alleles of the i-th edge
        ei_c <- sum(co_i)  # occurrence count of the i-th edge
        for (j in (i + 1) : n) {
            ej <- ep[j, ]
            ej_id <- ej$edge_ID
            ej_x <- ej$x
            ej_y <- ej$y
            co_j <- as.logical(apam[, ej_x]) & as.logical(apam[, ej_y])  # co-occurrence of x, y alleles of the j-th edge
            ej_c <- sum(co_j)  # occurrence count of the j-th edge
            s <- round(sum(co_i & co_j) / min(ei_c, ej_c), digits = 8)
            sim[ei_id, ej_id] <- s
            sim[ej_id, ei_id] <- s  # This is a symmetric matrix.
        }
    }

    return(list(nwk = nwk, ep = ep, sim = sim))
}
