#' @title Count neighbours per node
#'
#' @description Count how many neighbours does each node have in a network. The
#' function treats the network as undirected.
#'
#' @param net A Graph object, which can be created by the function mkNetwork.
#' @param groups A data frame mapping alleles to groups, such as genes. The first
#' column stores allele names and the second column stores group information.
#' @param na A character argument specifying the default group in the output.
#'
#' @return When groups != NULL, the function returns a list of two data frames
#' for allele-level and gene-level neighbour counts, respectively. Otherwise, a
#' single data frame of allele-level neighbour counts is returned.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 9 Sep 2018, the lastest edition: 29 Oct 2018

countNeighbours <- function(net, groups = NULL, na = "Mix") {
    require(igraph)  # Notice it overrides V, E methods in the package GeneMates.

    # Deduplicate rows when treating them as undirected
    E <- GeneMates::E(net)
    if (! "pair" %in% names(E)) {
        E <- assignPairID(lmms = E, from = 1, paired.rows = FALSE)
    }
    E <- E[!duplicated(E$pair), ]

    # Allele-level counts
    a <- graph_from_data_frame(d = E, directed = FALSE, vertices = GeneMates::V(net))
    a_deg <- degree(graph = a, mode = "all")  # for an undirected network, all = in = out degrees
    a_deg <- data.frame(Allele = names(a_deg), Neighbours = as.integer(a_deg), stringsAsFactors = FALSE)
    a_deg <- a_deg[order(a_deg$Neighbours, decreasing = TRUE), ]

    # Gene-level counts
    if (!is.null(groups)) {
        a_deg$Group <- groups[, 2][match(a_deg$Allele, groups[, 1])]
        a_deg$Group[is.na(a_deg$Group)] <- na
        a_deg <- a_deg[, c("Group", "Allele", "Neighbours")]
        groups <- sort(unique(a_deg$Group), decreasing = FALSE)
        g_deg <- data.frame(Gene = character(0), Neighbours = integer(0), stringsAsFactors = FALSE)
        for (g in groups) {
            a_deg_sub <- subset(a_deg, Group == g)
            g_deg <- rbind.data.frame(g_deg,
                                      data.frame(Group = g, Neighbours = sum(a_deg_sub$Neighbours),
                                                 stringsAsFactors = FALSE))  # may have duplicated counts for the same neighbour gene
        }
        g_deg <- g_deg[order(g_deg$Neighbours, decreasing = TRUE), ]
        out <- list(a = a_deg, g = g_deg)
    } else {
        out <- a_deg
    }

    return(out)
}
