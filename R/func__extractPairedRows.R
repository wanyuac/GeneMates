#' @title Extract paired and unpaired rows from a data frame.
#'
#' @description Assuming a data frame with a column storing pair IDs, this function extracts rows sharing the same pair
#' IDs from the data frame. It can be used to isolate symmetric associations from asymmetric associations.
#' Specifically, a pair of random variables X and Y is considered as in a symmetric association if both
#' models Y ~ X and X ~ Y display a significant effect of the explanatory variable, respectively. On the
#' contrary, X and Y is in an asymmetric association if either model shows insignificance.
#'
#' @param x A data frame containing pair IDs.
#' @param pair.col A character or integer specifying which column is providing pair IDs.
#'
#' @examples
#' a <- rbind.data.frame(subset(assoc[["assoc"]], score >= 1.5 | (s_d == 1 & dif == 0)), stringsAsFactors = FALSE)
#' ap <- extractPairedRows(x = a, pair.id = "pair")
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the lastest edition: 16 August 2017

extractPairedRows <- function(x, pair.col = "pair") {
    pair.ids <- x[, pair.col]
    pair.count <- table(pair.ids)
    rows.paired <- as.integer(names(pair.count)[as.integer(pair.count) > 1])  # pair IDs occurring multiple times
    rows.single <- as.integer(names(pair.count)[as.integer(pair.count) == 1])  # pair IDs occurring only once in the data frame x
    is.paired <- pair.ids %in% rows.paired

    return(list(paired = x[is.paired, ], single = x[!is.paired, ]))
}
