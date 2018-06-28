#' @title Mapping a vector of numerics to vertex sizes through a linear transformation
#'
#' @description This function generates vertex sizes for plotting networks based on a vertex attribute.
#' It works in the same way as the Cytoscape in mapping continuous values to an attribute. This is a
#' "smoother" way than simply multiply values with a constant.
#'
#' @param x Input numeric vector to be transformed.
#' @param x.min Minimum of x or the variable from which the x is drawn.
#' @param x.max Maximum of x or the variable from which the x is drawn.
#' @param size.min Minimal vertex size.
#' @param size.max Maximal vertex size.
#'
#' @return A vector of sizes with a precision of two decimals each.
#'
#' @examples
#' library(igraph)
#' ...  # network construction
#' assoc <- findPhysLink(...)
#' V(net)$size <- vertexAttr2Size(x = V(net)$n, x.min = 2, x.max = max(assoc[["mapping"]]$count),
#' size.min = 15, size.max = 60)
#'
#' @author Yu Wan (\email{wanyuac@gmail.com})
#' @export vertexAttr2Size
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition and the lastest edition: 10 Oct 2017

vertexAttr2Size <- function(x, x.min = min(x), x.max = max(x), size.min = 15, size.max = 60) {
    size <- size.min + (x - x.min) / (x.max - x.min) * (size.max - size.min)

    return(round(size, digits = 2))
}
