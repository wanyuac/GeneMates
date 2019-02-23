#' @title Summarise allele content per clique
#'
#' @description This function extract allele names from cliques and return a
#' summary data frame.
#'
#' @param q A list produced by the function max_cliques of the package igraph.
#' @param s An optional list produced by the function compEdgeOccur for edge similarities.
#'
#' @return A data frame of three columns: Clique, Size and Alleles.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export summariseCliques
#'
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the latest edition: 28 Oct 2018

summariseCliques <- function(q, s = NULL) {
    require(data.table)

    d <- as.data.frame(rbindlist(lapply(q, .makeLines, sim = s)))
    d <- d[order(d$Size, decreasing = TRUE), ]  # By default, the cliques functions orders cliques in an ascending order of clique sizes.
    d <- cbind.data.frame(Clique = 1 : nrow(d), d)  # assign arbitrary, continuous clique IDs

    return(d)
}

.makeLines <- function(i, sim = NULL) {  # A subordinate function of summariseCliques
    # Basic summary
    alleles <- i$name
    d <- data.frame(Size = length(alleles), Alleles = paste(alleles, collapse = ","),
                    stringsAsFactors = FALSE)

    # Summarise Edge similarities
    if (is.list(sim)) {  # assuming elements "es" and "sim" are included
        es <- subset(sim[["es"]], y %in% alleles & x %in% alleles) # extract undirected, non-duplicated edges for the current clique
        eids <- es$edge_ID
        eid_n <- length(eids)
        if (eid_n > 2) {  # eid_n >= 3 when length(alleles) >= 3.
            sim_mat <- sim[["sim"]][eids, eids]  # extract the relevant similarity matrix
            sim_coes <- numeric(0)
            for (i in 1 : (eid_n - 1)) {  # sim_coe has a minimum row number of three.
                for (j in (i + 1) : eid_n) {
                    sim_coes <- append(sim_coes, sim_mat[eids[i], eids[j]])
                }
            }
            # Actually, for any Uni_edge > 0, S_num = Uni_edge * (Uni_edge - 1) / 2.
            # Similarly, Uni_edges = Alleles * (Alleles - 1) / 2.
            # However, I put Uni_edges and S_num here for users to check their outputs.
            d <- cbind.data.frame(d, data.frame(Uni_edges = eid_n, S_num = length(sim_coes),
                                                S_min = min(sim_coes), S_med = median(sim_coes),
                                                S_max = max(sim_coes), stringsAsFactors = FALSE))  # Uni_edges: number of unidirected edges
        } else {  # eid_n = 1 when length(alleles) == 2, for which no edge similarity coefficient is applicable.
            d <- cbind.data.frame(d, data.frame(Uni_edges = eid_n, S_num = 0,
                                                S_min = NA, S_med = NA, S_max = NA,
                                                stringsAsFactors = FALSE))
        }
    }

    return(d)  # returns a single-line data frame
}
