#' @title Summarise allele content per clique
#'
#' @description This function extract allele names from cliques and return a
#' summary data frame.
#'
#' @param q A list produced by the function max_cliques of the package igraph.
#'
#' @return A data frame of three columns: Clique, Size and Alleles.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export summariseCliques
#'
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the latest edition: 27 Oct 2018

summariseCliques <- function(q) {
    require(data.table)

    d <- as.data.frame(rbindlist(lapply(q, .makeLines)))
    d <- d[order(d$Size, decreasing = TRUE), ]  # By default, the cliques functions orders cliques in an ascending order of clique sizes.
    d$Clique <- 1 : nrow(d)
    d <- d[, c("Clique", "Size", "Alleles")]

    return(d)
}

.makeLines <- function(i) {
    # a subordinate function of summariseCliques
    alleles <- i$name
    d <- data.frame(Alleles = paste(alleles, collapse = ","), Size = length(alleles),
                    stringsAsFactors = FALSE)

    return(d)
}
