#' @title Attach new columns summarising distance measurements to LMM results (lmms)
#'
#' @description This function produces summary statistics of distance measurements
#' per pair of alleles.
#'
#' @param lmms a data frame that must contain two columns x and y for allele names.
#' @param lmms.ds a data frame sharing the columns x and y of lmms and all of
#' their distance measurements.
#' @param tree A phylogenetic tree (an object of the phylo class defined in the
#' ape package) accompanying the clades argument.
#' @param clades A PAM of samples in every clade. It can be obtained using the
#' function tree2Clade.
#' @param allele.pam an allelic PAM of bacterial genes. Usually it is the element
#' A in the output list of the function importAllelicPAM.
#' @param max.nodes an integer specifying the maximal number of nodes to which
#' we can trust paths for distance measurements. Set max.nodes = -1 to skip this
#' fiter when the accuracy of measurements does not matter. Notice node number
#' equals zero for overlapping genes.
#' @param max.dist An inclusive upper bound for filterring distance measurements.
#' Measurements above this threshold will be ignored. Set it to Inf to turn off
#' this filter.
#' @param colname.co a column name in the data frame lmms for the number of co-occurrence
#' @param colname.pair a column name in the data frame lmms for indices of a pair
#' of alleles
#' @param qs quantile probabilities (between 0 and 1, inclusive) of physical
#' distances. The minimum, first quantile (25th percentile), median, third quantile
#' (75th percentile) and the maximum will always be calculated.
#' @param n.cores Number of computational cores that will be used in parallel for
#' this function.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2017-2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First editions: 5 July 2017 - 28 July 2017; the latest edition: 1 Apr 2018

summariseDist <- function(lmms, lmms.ds, tree, clades, allele.pam,
                          max.nodes = -1, max.dist = Inf,
                          colname.pair = "pair", colname.co = "n_xy",
                          qs = c(0, 0.25, 0.5, 0.75, 1),
                          n.cores = n.cores) {
    # Initialisation
    print(paste0(Sys.time(), ": summarising distance measurements for LMM outputs."))
    groups <- .checkGroupIDs(lmms, lmms.ds)  # Keys of both lists must be the same.
    qs <- .checkQuantileSetting(qs)  # The outcome is sorted.
    qs.n <- length(qs)  # no. of quantiles

    # Construct the result list
    summary.stats <- vector(mode = "list", length = length(groups))
    names(summary.stats) <- groups  # Each element is a NULL initially.

    # Obtain summary statistics of distance measurements for allele pairs from different groups in the LMM outputs under H1
    for (g in groups) {
        # Get the LMM results
        h1 <- lmms[[g]][["h1"]]  # Association under the alternative hypothesis

        # Rename columns for the convenience of programming
        index.pair <- 0
        index.co <- 0
        if (colname.pair != "pair") {
            index.pair <- which(names(h1) == colname.pair)
            names(h1)[index.pair] <- "pair"
        }
        if (colname.co != "n_xy") {
            index.co <- which(names(h1) == colname.co)
            names(h1)[index.co] <- "n_xy"
        }
        .checkLmms(lmms = h1, group = g)  # stop if failed
        h1.names <- names(h1)  # all column names of the data frame h1
        h1 <- h1[, c("pair", "y", "x", "y_pat", "x_pat", "n_xy",
                     setdiff(h1.names, c("pair", "y", "x", "y_pat", "x_pat", "n_xy")))]  # rearrange columns

        # Get the distance measurements
        ds <- lmms.ds[[g]]  # physical distances between alleles tested for associations
        ds.n <- nrow(ds)  # ds links to h1 via pair IDs under each g
        if (ds.n == 0) {  # Normally, this does not happen.
            # The downstream function .calcSummaryStats can handle this situation.
            print(paste0("[summariseDist] warning: lmms.ds[[", g, "]] is empty."))
        } else {  # Apply a series of filters to the measurements
            if (max.nodes >= 0) {
                ds <- subset(ds, node_number <= max.nodes)  # remove distances measured in paths that have too many nodes
                ds.n.flt <- nrow(ds)
                ds.rm <- ds.n - ds.n.flt  # number of rows removed from ds
                print(paste(as.character(ds.rm), "distance measurements are removed from",
                            g, "as they involve more than", max.nodes, "nodes each.",
                            sep = " "))
                ds.n <- ds.n.flt  # refresh the row count
                if (ds.n.flt == 0) {
                    print(paste("Warning: no rows of", g, "left after filtering for a maximum node number.",
                                sep = " "))
                }
            }
            if (ds.n > 0 & max.dist != Inf & max.dist >= 0) {
                ds <- subset(ds, distance <= max.dist)  # remove measurements that are way too large
                ds.n.flt <- nrow(ds)
                ds.rm <- ds.n - ds.n.flt
                print(paste(as.character(ds.rm), "distance measurements are removed from",
                            g, "as they are greater than", max.dist, "bp.", sep = " "))
                if (ds.n.flt == 0) {
                    print(paste("Warning: no rows of", g, "left after filtering for a maximum distance.",
                                sep = " "))
                }
            }  # The filters finish here.
        }

        # Generate column names for quantiles of the distance measurements
        qnames <- paste("d_P", as.character(qs * 100), sep = "")  # d_P0, d_P25, ...
        qnames[1] <- "d_min"  # particularly rename the first, middle and the last quantiles
        qnames[which(qnames == "d_P50")] <- "d_median"
        qnames[qs.n] <- "d_max"

        # Calculate summary statistics for every pair of rows in h1
        # Notice ds may equal a data frame of no rows, although this is an
        # abnormal condition.
        h1 <- .getDistQuantiles(lmms = h1, lmms.ds = ds, qs = qs, qnames = qnames,
                                pam = allele.pam, tree = tree, clades = clades,
                                n.cores = n.cores)

        # Restore column names
        if (index.co > 0) {
            names(h1)[which(names(h1) == "n_xy")] <- colname.co
        }
        if (index.pair > 0) {
            names(h1)[which(names(h1) == "pair")] <- colname.pair
        }

        # Save results for the current association group
        summary.stats[[g]] <- h1
    }

    return(summary.stats)  # This list may contain NULL elements.
}
