#' @title Summarise allelic physical distances for edges in an association network
#'
#' @description This function summarises allelic physical distances for a table
#' of edges. The result shows how does the incorporation of physical distances
#' affect the edge weights (weighted distance score: w_d).
#'
#' @param E An edge list having columns in the order: allele_1, allele_2, co-occurrence
#' count and the distance score s_d. An additional pair column can be included.
#' @param ds A data frame of allelic physical distances imported from the output of
#' the pipeline physDist.
#' @param d.max Maximum distance to be considered as accruate.
#' @param n.max Maximum node number for distances that are considered as accruate.
#' @param source.graph Name for assembly graphs as a source of distance measurements.
#' @param source.contig Name for contigs as a source of distance measurements.
#' @param source.complete Name for finished-grade genomes as a source of distance measurements.
#' @param sort.output Keep it TRUE to enable sorting of the output data frame in
#' a descending order of the measurability (Mr) and count of reliable distances.
#'
#' @return A data frame of the following columns:
#'   Allele_1, Allele_2: names of associated alleles, ordered alphabetically;
#'   S_d, Co: distance score s_d and co-occurrence count;
#'   M: overall measurability of physical distances based on the co-occurrence count;
#'   Mr: measurability of reliable distances;
#'   N: overall count of physical distances;
#'   Nr: number of all reliable distances;
#'   Ng: overall number of distances from assembly graphs;
#'   Ng_r: number of reliable distances from assembly graphs;
#'   Nc_r: number of reliable distances from contigs;
#'   Nf: number of distances from finished-grade genomes.
#'
#' @note Since the distance measurements may be prioritised according to their
#' sources, Ng and Ng_r may not be accurate when Nc_r or Nf > 0; Nc_r may not be
#' accurate when Nf > 0.
#'
#' @examples
#'   assoc_lmm <- findPhysLink(...)
#'   a_lmm_dif <- subset(assoc_lmm$assoc, beta > 0 & p_adj <= 0.05)
#'   ds_stats <- summariseDistsForEdges(E = a_lmm_dif[, c("pair", "y", "x", "n_xy", "s_d")],
#'                                      ds = assoc_lmm$ds, d.max = 250e3, n.max = 2,
#'                                      source.graph = "graph", source.contig = "contig",
#'                                      source.complete = NA, sort.output = TRUE)
#'   ds_stats <- ds_stats[, c("Allele_1", "Allele_2", "Co", "S_d", "M", "Mr", "N", "Nr", "Nc_r")]  # For prioritised distances, Nc_r and Ng are mutually exclusive.
#'
#'
#' @author Yu Wan, \email{wanyuac@@gmail.com}
#' @export summariseDistsForEdges
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 10 Oct 2017, the lastest edition: 21 Apr 2019

summariseDistsForEdges <- function(E, ds, d.max = 250e3, n.max = 2, source.graph = "graph",
                                   source.contig = "contig", source.complete = "complete",
                                   sort.output = TRUE) {
    require(data.table)

    # E: an edge list of allele names, beta, co-occurrence count, distance score s_d (c * m_in) and pair IDs.
    # ds: a data frame of allelic physical distances.
    # 1. Convert the edge list into an undirected network
    if (is.na(source.graph) && is.na(source.contig && is.na(source.complete))) {
        stop("Argument error: source names must not be all NA's.")
    }
    if (! "pair" %in% names(E)) {
        E <- E[, 1 : 4]
        names(E) <- c("a1", "a2", "n_co", "s_d")
        E <- assignPairID(lmms = E, from = 1, paired.rows = FALSE)  # append a column pair to this data frame
        E <- E[, c("pair", "a1", "a2", "n_co", "s_d")]
    } else {
        names(E)[2 : 5] <- c("a1", "a2", "n_co", "s_d")  # assuming the first column is "pair"
    }
    E <- E[!duplicated(E$pair), ]
    E <- E[, c("a1", "a2", "n_co", "s_d")]

    # 2. Summarise physical distances
    out_list <- mapply(.distsOfEdge, E$a1, E$a2, E$n_co, E$s_d,
                       MoreArgs = list(ds = ds, d.max = d.max, n.max = n.max,
                                       source.graph = "graph", source.contig = "contig",
                                       source.complete = "complete"),
                       SIMPLIFY = FALSE, USE.NAMES = TRUE)  # SIMPLIFY = TRUE: does not return a data frame for rbindlist
    out <- as.data.frame(rbindlist(out_list))

    # 3. Remove uninformative columns
    if (is.na(source.graph)) {
        out <- out[, -which(names(out) %in% c("Ng", "Ng_r"))]
    }
    if (is.na(source.contig)) {
        out <- out[, -which(names(out) == "N_cr")]
    }
    if (is.na(source.complete)) {
        out <- out[, -which(names(out) == "Nf")]
    }

    # 4. Sort the data frame according to the measurability of reliable distances
    if (sort.output) {
        out <- out[order(out$Mr, out$Nr, decreasing = TRUE), ]
    }

    return(out)
}

.distsOfEdge <- function(a1, a2, co, s_d,
                         ds, d.max, n.max, source.graph,
                         source.contig, source.complete) {
    alleles <- sort(c(a1, a2), decreasing = FALSE)  # sort allele names alphabetically
    a1 <- alleles[1]
    a2 <- alleles[2]
    ds <- getRowsXY(df = ds, k1 = a1, k2 = a2, c1 = "query1", c2 = "query2")
    n <- nrow(ds)
    if (n > 0) {
        m <- round(n / co, digits = 4)  # conditional measurability of all distance measurements given the co-occurrence count
        if (!is.na(source.graph)) {
            ng <- sum(ds$source == source.graph)  # number of distances measured in assembly graphs
        } else {
            ng <- 0
        }

        # get reliable measurements
        if (is.na(source.complete)) {  # no distances from complete genomes
            ds_reli <- subset(ds, distance <= d.max & node_number <= n.max)  # reliable distance measurements
        } else {
            sel <- ds$source == source.complete
            ds1 <- ds[sel, ]  # distances from complete genomes
            ds2 <- ds[!sel, ]  # distances from draft genomes
            ds_reli <- rbind.data.frame(ds1, subset(ds2, distance <= d.max & node_number <= n.max),
                                        stringsAsFactors = FALSE)
        }

        # summarise reliable distances
        nr <- nrow(ds_reli)  # number of reliable distances
        if (nr > 0) {
            mr <- round(nr / co, digits = 4)  # conditional measurability of reliable distances given the co-occurrence count
            source_tab <- table(ds_reli$source)
            source_list <- names(source_tab)
            if (!is.na(source.graph) && source.graph %in% source_list) {
                ngr <- source_tab[[source.graph]]
            } else {
                ngr <- 0
            }

            if (!is.na(source.contig) && source.contig %in% source_list) {
                ncr <- source_tab[[source.contig]]
            } else {
                ncr <- 0
            }

            if (!is.na(source.complete) && source.complete %in% source_list) {
                nf <- source_tab[[source.complete]]
            } else {
                nf <- 0
            }
        } else {
            mr <- ngr <- ncr <- nf <- 0
        }
    } else {
        m <- nr <- mr <- ng <- ngr <- ncr <- nf <- 0
    }

    out <- data.frame(Allele_1 = a1, Allele_2 = a2, Co = co, S_d = s_d,
                      M = m, Mr = mr, N = n, Nr = nr, Ng = ng, Ng_r = ngr,
                      Nc_r = ncr, Nf = nf, stringsAsFactors = FALSE)

    return(out)
}
