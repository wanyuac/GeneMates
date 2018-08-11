# Defining data types and functions for the networkAnalyser module
# Export the following names to NAMESPACE, otherwise the methods cannot be called
# by users even though they are methods of classes.
#
# Example:
#     g <- new("Graph", id = "graph1", V = nwk$V, E = nwk$E)
#     gs <- new("GraphSet", V = g@V, E = list("a" = g@E))
#     getV_id(gs, id = "a")
#     g1 <- getGraph(gs, id = "a")
#     getV(g1)
#
# Copyright 2017-2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 10 Aug 2018; the latest update: 11 Aug 2018
#
#' @export Graph
#' @export GraphSet
#' @export getV
#' @export countEdges
#' @export getV_id
#' @export countEdges_id
#' @export getGraph

# Class Graph ###############
# Store a single graph (including a network)
# Slots: V (vertices) and E (edges).
Graph <- setClass(
    Class = "Graph",
    slots = list(id = "character", V = "data.frame", E = "data.frame"),
    prototype = list(id = NULL, V = NULL, E= NULL),

    validity = function(object) {
        v <- is.data.frame(object@V) && is.data.frame(object@E)
        if (v) {
            ve <- union(object@E[, 1], object@E[, 2])
            vv <- object@V[, 1]

            # E must not contain any vertices that are not present in V, and vice versa.
            unmatch <- any(! ve %in% vv) || any(! vv %in% ve)
            v <- v && !unmatch
        }

        return(v)
    }
)

# Define methods of Graph ===============
setGeneric(name = "getV",
           # Do not use "signature = Graph" here.
           def = function(object) {  # Declaring the default value is necessary here, otherwise, the default value does not work.
               standardGeneric(f = "getV")
           }
)

setMethod(f = "getV",
          signature = "Graph",
          definition = function(object) {  # return vertex names
              return(object@V[, 1])
          }
)

setGeneric(name = "countEdges",
           def = function(object) {
               standardGeneric(f = "countEdges")
           }
)

setMethod(f = "countEdges",
          signature = "Graph",
          definition = function(object) {

              return(nrow(object@E))
          }
)


# Class GraphSet ###############
# This class is used for storing all graphs that share nodes from the same set.
# It is useful for saving all cliques or communities identified in a network.
#
# The class has two slots:
#   V: a data frame for attributes of nodes in all member graphs
#   E: a list (named or unnamed) of edge tables. Each table defines a member
#      graph in this stack.
#
# Define a generator function for objects of the class GraphSet  ===============
GraphSet <- setClass(
    # class name
    Class = "GraphSet",

    # Data slots of every object. E can be a named list.
    slots = list(V = "data.frame", E = "list"),

    # default values of slots
    # By far, this class does not provide a function to check if there are
    # redundant nodes (namely, not present in any member graphs) or missing nodes
    # in the V slot in accordance with the E slot. So it is the user's responsibility
    # to ensure the nodes in V matches those in graphs.
    prototype = list(V = NULL, E = NULL)
)

# Check if arguments for initialisation are correct ===============
.checkGraphObject <- function(object) {
    if ((is.null(object@V) | is.data.frame(object@V) | is.character(object@V)) & (is.null(object@E) | is.list(object@E))) {
        valid = TRUE
    } else {
        valid = "Error: V can only be a data frame, a character vector or NULL; E can only be a list or NULL."
    }

    return(valid)
}

# Register the function .checkGraphObject as a validation function of the class Graph
setValidity(Class = "GraphSet", method = .checkGraphObject)

# Define methods for the class GraphSet ===============
# reference: http://www.cyclismo.org/tutorial/R/s4Classes.html
# Get vertex names ---------------
#
# 1. reserve a name for a generic function
setGeneric(name = "getV_id",
           def = function(object, id, sorted = TRUE) {  # Declaring the default value is necessary here, otherwise, the default value does not work.
               standardGeneric(f = "getV_id")
           }
)

# 2. define the actual method using the reserved name
# Retrieve vertices of a specific member graph (given an ID)
setMethod(f = "getV_id",
          definition = function(object, id, sorted = TRUE) {  # cluster.id: an integer or a character value for the cluster index or name, respectively.
              edges <- object@E[[id]]
              vertices <- union(edges[, 1], edges[, 2])  # discard duplicated values
              if (sorted) {
                  vertices <- sort(vertices, decreasing = FALSE)
              }

              return(vertices)
          }
)

# Count the number of edges in a given member graph ---------------
setGeneric(name = "countEdges_id",
           def = function(object, id) {
               standardGeneric(f = "countEdges_id")
           }
)

setMethod(f = "countEdges_id",
          signature = "GraphSet",
          definition = function(object, id) {

              return(nrow(object@E[[id]]))
          }
)

# Extract a member graph from the set given its ID ---------------
# Return a Graph object
setGeneric(name = "getGraph",
           def = function(object, id) {
               standardGeneric(f = "getGraph")
           }
)

setMethod(f = "getGraph",
          signature = "GraphSet",
          definition = function(object, id) {
              e <- object@E[[id]]
              v <- union(e[, 1], e[, 2])
              v_tab <- object@V[object@V[, 1] %in% v, ]

              return(new("Graph", id = id, V = v_tab, E = e))
          }
)
