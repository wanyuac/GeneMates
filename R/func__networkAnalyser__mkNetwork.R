#' @title Create network data from association results.
#'
#' @description This function generates edge and node tables for a network. A user
#' may print its output as an edge file and a node file for network construction.
#'
#' @param assoc A data frame from which the network is created. It may be pre-filtered
#' for some criteria by users.
#' @param name.x Column name (character) for the explanatory variable X.
#' @param name.y Column name (character) for the response variable Y.
#' @param edge.attr A character vector of column names for edge attributes.
#' @param node.x.attr A character vector of column names for node attributes of X.
#' @param node.y.attr A character vector of column names for node attributes of Y.
#' Vectors node.x.attr and node.y.attr must have the same length.
#' @param node.attr.names A character vector for column names of the final data frame
#' of nodes. It should have one more name (the "name" column for node names) than
#' do the vectors node.x.attr and node.y.attr.
#' @param id (optional) Network name
#'
#' @return G A Graph object having three slots: id, E and V. E: a data frame for
#' edges and their attributes; V, a data frame for nodes and their attributes; id,
#' an optional parameter for the network name.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' nwk <- mkNetwork(assoc = assoc[["assoc"]], name.x = "x", name.y = "y",
#' edge.attr = c("beta", "score", "m", "s_d"), node.y.attr = "n_y",
#' node.x.attr = "n_x", node.attr.names = c("allele", "n"), id = "linkage")
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the lastest edition: 13 Mar 2018

mkNetwork <- function(assoc = NULL, name.x = "x", name.y = "y", edge.attr = NULL,
                      node.x.attr = NULL, node.y.attr = NULL,
                      node.attr.names = "name", id = NULL) {
    # Sanity check
    if (length(node.x.attr) != length(node.y.attr)) {
        stop("Error: the vectors node.x.attr and node.y.attr do not have the same length.")
    }
    names.all <- c(name.x, name.y, edge.attr, node.x.attr, node.y.attr)  # NULL values are dropped out automatically.
    i <- which(!(names.all %in% names(assoc)))
    if (length(i) > 0) {
        stop(paste("Error: column names", paste(names.all[i], collapse = ","),
                   "are not present in the assoc data frame.", sep = " "))
    }

    # Create the data frame for edges
    E <- assoc[, c(name.y, name.x, edge.attr)]

    # Create the data frame for nodes
    vx <- assoc[, c(name.x, node.x.attr)]
    vy <- assoc[, c(name.y, node.y.attr)]
    names(vx) <- node.attr.names
    names(vy) <- node.attr.names
    V <- rbind.data.frame(vx, vy, stringsAsFactors = FALSE)
    V <- unique(V)  # de-duplicate rows

    return(new("Graph", id = id, E = E, V = V))
}
