#' @title Compile graphs in a GraphSet object into a single network while retaining
#' separation between individual graphs.
#'
#' @description This function pools graphs in a Graph object into a single network
#' while keeping individual graphs separated. It is particularly useful when a
#' user wants to display cliques or other kinds of clusters in a single view (such
#' as in Cytoscape). The Graph object is produced by the function extractSubgraphs.
#'
#' This function does not behave like the Cytoscape utility "merge networks".
#' Specifically, it keeps the separation between graphs through appending a cluster
#' ID to each node name so that the node shared by different graphs become distinguishable.
#' It adds a label column to save the node names. In Cytoscape, for example, the
#' user can choose this label column as node labels. Notice the first two columns,
#' namely, node1 and node2, correspond to allele names in the third and fourth
#' columns, respectively. Therefore, to make a directed network in Cytoscape, the
#' source and target columns are specified based on the third and fourth columns.
#'
#' @param graphs A Graph object created using the function extractSubgraphs.
#' @param edge.attr A character vector of column names for edge attributes.
#' @param node.attr A character vector of column names for node attributes.
#'
#' @examples g <- compileGraphs(graphs = clusters, edge.attr = c("beta", "score"),
#' node.attr = c("n", "class"))
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First version: 14 Mar 2018; the lastest edition: 10 Aug 2018
#  In memorial to the physicist Stephen Hawking, 8/1/1942 - 14/3/2018 (aged 76).

compileGraphs <- function(graphs, edge.attr = NULL, node.attr = NULL) {
    # Sanity check
    if (as.character(class(graphs)) != "GraphSet") {
        stop("Argument error: graphs must be an object of the class GraphSet.")
    }

    # Create an edge table and a node table
    V <- NULL
    E <- NULL
    cids <- names(graphs@E)  # cluster IDs
    prev.keys <- names(graphs@E[[1]])[c(1, 2)]  # Column names must be the same for all data frames in the list graphs@E.
    prev.node.col <- names(graphs@V)[1]  # the column for node names
    names(graphs@V)[1] <- "label"
    for (cid in cids) {
        Ec <- graphs@E[[cid]]  # Ec is a data frame (under the current cid).
        names(Ec)[c(1, 2)] <- c("label1", "label2")  # temporarily rename the keys to simplify the script
        Ec <- Ec[, c("label1", "label2", edge.attr)]  # I assume Ec always has at least two rows.
        Vc <- subset(graphs@V, label %in% union(Ec$label1, Ec$label2))
        Vc <- Vc[, c("label", node.attr)]

        # Make the edge table
        Ec <- cbind.data.frame(node1 = paste(Ec$label1, cid, sep = "__"),
                               node2 = paste(Ec$label2, cid, sep = "__"),
                               Ec,
                               cid = rep(cid, times = nrow(Ec)),
                               stringsAsFactors = FALSE)
        names(Ec)[c(3, 4)] <- prev.keys  # restore the orignial column names
        E <- rbind.data.frame(E, Ec, stringsAsFactors = FALSE)

        # Make the node table
        Vc <- cbind.data.frame(node = paste(Vc$label, cid, sep = "__"), Vc,
                               cid = rep(cid, times = nrow(Vc)),
                               stringsAsFactors = FALSE)
        names(Vc)[2] <- prev.node.col
        V <- rbind.data.frame(V, Vc, stringsAsFactors = FALSE)
    }

    # return a Graph object
    return(new("Graph", id = NULL, E = E, V = V))
}
