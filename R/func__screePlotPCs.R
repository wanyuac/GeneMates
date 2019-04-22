#' @title Calculating proportion of genetic variation captured by each principal
#' component
#'
#' @description Since eigenvalues is proportional to the amount of genetic variation
#' captured by the corresponding principal component (PC), this function computes
#' the proportion for every PC so that users can draw a scree plot for PCs.
#'
#' @param ev A vector of eigenvalues, which can be obtained from the element
#' struc$C$ev in the output of the functions lmm or findPhysLink.
#' @param ev.unsorted A logical argument telling the function that eigenvalues
#' have been sorted in a descending order. Default: FALSE.
#'
#' @examples pr <- screePlotPCs(ev = assoc_lmm$struc$C$ev)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the lastest edition: 29 May 2018

screePlotPCs <- function(ev, ev.unsorted = FALSE) {
    # ev: all eigenvalues, pre-sorted in an ascending order
    # npcs: number of principal components to be considered
    # First, calculate how much variantion does each projection axis captures.
    if (ev.unsorted) {
        ev <- sort(ev, decreasing = TRUE)
    }
    n <- length(ev)
    axes <- 1 : n
    pr <- data.frame(Axis = axes, Ev = ev, Pr = numeric(n), Pr_cum = numeric(n))
    ev.total <- sum(ev)
    pr$Pr <- round(pr$Ev / ev.total * 100, digits = 4)  # proportion of variation captures on each axis of sample projections

    # Then calculate the cumulative proportion of variation captured
    for (i in axes) {
        pr$Pr_cum[i] <- round(sum(pr$Ev[1 : i]) / ev.total * 100, digits = 4)
    }

    return(pr)
}
