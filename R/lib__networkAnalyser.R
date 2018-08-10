# Defining data types and functions for the networkAnalyser module
# Export the following names to NAMESPACE, otherwise the methods cannot be called
# by users even though they are methods of classes.
#
# Copyright 2017-2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 10 Aug 2018; the latest update: 10 Aug 2018
#
#' @export getV
#' @export countEdges

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
.createGraphSet <- setClass(
    # class name
    Class = "GraphSet",

    # data slots of every object
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
setGeneric(name = "getV",
           def = function(object, cluster.id, sorted = TRUE) {  # Declaring the default value is necessary here, otherwise, the default value does not work.
               standardGeneric(f = "getV")
           }
)

# 2. define the actual method using the reserved name
# Retrieve vertices of a specific member graph
setMethod(f = "getV",
          signature = "GraphSet",
          definition = function(object, cluster.id, sorted = TRUE) {  # cluster.id: an integer or a character value for the cluster index or name, respectively.
              edges <- object@E[[cluster.id]]
              vertices <- union(edges[, 1], edges[, 2])  # discard duplicated values
              if (sorted) {
                  vertices <- sort(vertices, decreasing = FALSE)
              }

              return(vertices)
          }
)

# Count the number of edges in a given member graph ---------------
setGeneric(name = "countEdges",
           def = function(object, cluster.id) {
               standardGeneric(f = "countEdges")
           }
)

setMethod(f = "countEdges",
          signature = "GraphSet",
          definition = function(object, cluster.id) {

              return(nrow(object@E[[cluster.id]]))
          }
)
