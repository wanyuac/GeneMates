#' @title Creating objects of the class "Graph" for subgraphs of a large graph
#'
#' @description This is a generator function used for creating objects of subgraphs. Particularly in this package, a
#' subgraph refers to a connected graph component in a linkage network.
#'
#' @param V A data frame where the first column contains vertex names and the other columns (optional) are vertex attributes.
#' @param E A data frame consisting of all edges of an overall graph. The first and second columns are considered
#' as "from" and "to" vertices of directed edges. Other columns are treated as weights of edges.
#' @param clusters A data frame of three columns for cluster IDs (the first column), "from" vertices (the second column) and
#' "to" vertices (the third column). Recommend to name columns of clusters by c("cluster", "from", "to") for convenience of
#' understanding, although column names are not used in this function. Notice for overlapping clusters, some edges may belong
#' to multiple clusters. Cluster IDs will be treated as characters in this function.
#'
#' @note vertices in V are not necessarily the same as those in E. The idea to set a slot V is to keep vertex attributes
#' vertices as well as to deal with singleton.
#'
#' @examples g <- extractSubgraphs(V = nodes, E = edges, clusters = c)
#'
#' @return An object of the class "Graph", which has two slots:
#' V: A data frame for vertices and vertex attributes.
#' E: A data frame for edges and edge attributes.
#' Method getV(object, cluster.id, sorted = TRUE), which extracts vertex names of the cluster 'cluster.id' from a Graph object.
#' The names are sorted in an ascending order alphabetically if sorted = TRUE.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First version: 3 Oct 2017; the lastest edition: 10 Aug 2018

extractSubgraphs <- function(V, E, clusters) {
    # sanity check
    if (!is.data.frame(E)) {
        stop("Error: E must be a data frame.")
    } else if (ncol(E) < 2) {
        stop("Error: E must have at least two columns for 'from' and 'to' vertices.")
    } else {
        prev.names <- names(E)[c(1, 2)]  # store column names given by users
        names(E)[c(1, 2)] <- c("from", "to")
    }

    if (is.character(V)) {
        V <- data.frame(vertex = V, stringsAsFactors = FALSE)  # no node attributes
    } else if (!is.data.frame(V)) {
        stop("Error: V must be either a single character vector or a data frame.")
    }

    if (!is.data.frame(clusters)) {
        stop("Error: 'clusters' must be a data frame.")
    } else if (ncol(clusters) < 3) {
        stop("Error: 'clusters' must have at least three columns for 'cluster', 'from' and 'to' values.")
    } else {
        names(clusters) <- c("cluster", "from", "to")
    }

    # make a named list for the slot E and a vector of vertex names
    cluster.ids <- unique(as.character(clusters[, 1]))
    cluster.edges <- vector(mode = "list", length = length(cluster.ids))
    names(cluster.edges) <- cluster.ids
    vertices <- NULL
    for (c in cluster.ids) {
        edges <- subset(clusters, cluster == c)
        vertices <- union(vertices, union(edges$from, edges$to))
        edges <- merge(x = edges, y = E, by = c("from", "to"), all.x = TRUE, all.y = FALSE, sort = FALSE)
        names(edges)[c(1, 2)] <- prev.names  # restore column names
        cluster.edges[[c]] <- edges
    }

    # pack data elements into G, a GraphSet object
    # see lib__networkAnalyser.R for the definition of .createGraphSet and the class
    # GraphSet.
    G <- .createGraphSet(V = V[V[, 1] %in% vertices, ], E = cluster.edges)

    return(G)
}
