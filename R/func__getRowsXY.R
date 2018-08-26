#' @title Extract rows matching to a pair of keys k1 and k2.
#'
#' @description This is a handy small tool to help investigating connections between two objects. A typical usage is to
#' extract or display results relating to alleles x and y.
#'
#' @param df a data frame whose rows are to be extracted.
#' @param k1 a character, numeric or integer type key being searched against the column c1 (see below).
#' @param k2 the same kind of key being searched against the column c2 (see below).
#' @param c1 name or index of the first value column (variable) in df. Default: "y".
#' @param c2 name or index of the second key column in df. Default: "x".
#'
#' @examples
#' assoc <- findPhysLink(...)
#' z <- getRowsXY(df = assoc[["assoc"]], k1 = "IntI1_1.1", k2 = "SulI_1616", c1 = "y", c2 = "x")
#' z <- getRowsXY(df = assoc[["assoc"]], k1 = "IntI1_1.1", k2 = "SulI_1616", c1 = 1, c2 = 2)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 20 August 2017, the lastest edition: 22 August 2017

getRowsXY <- function(df, k1 = NULL, k2 = NULL, c1 = "y", c2 = "x") {
    if (any(c(nrow(df) == 0, is.null(k1), is.null(k2)))) {
        stop("Error: arguments df, q1 or q2 are not completely specified.")
    } else {
        v1 <- df[, c1]
        v2 <- df[, c2]
        z <- df[(v1 == k1 & v2 == k2) | (v1 == k2 & v2 == k1), ]
    }

    return(z)
}
