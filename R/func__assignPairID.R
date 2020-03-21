#' @title Assign an ID to each pair of alleles or patterns.
#'
#' @description This is a generic function that assigns an ID to each pair of object names. For example, a single ID is
#' assigned to two rows in a data frame for associations between patterns (pat_1, pat_2) and (pat_2, pat_1).
#'
#' @param lmms A data frame whose first two columns (named y and x respectively) are used to assign pair IDs.
#' @param from The initial pair ID.
#' @param paired.rows: whether there are only two rows per pair of alleles, each in an inverse direction.
#' Set paired.rows = FALSE to process allele pairs in the data frame lmms.ds[["dif"/"idd"]].
#'
#' @examples
#' a0 <- assoc[["lmms.pat"]][["dif"]][["h1"]]
#' a0 <- assignPairID(lmms = a0, from = 1, paired.rows = TRUE)
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: March 2017, lastest edition: 15 August 2017

assignPairID <- function(lmms, from = 1, paired.rows = TRUE) {
    # check if y and x are column names
    if (any(!(c("y", "x") %in% names(lmms)))) {
        rename <- TRUE
        names.prev <- names(lmms)[1 : 2]
        names(lmms)[1 : 2] <- c("y", "x")
    } else {
        rename <- FALSE
    }

    # assign pair IDs
    new.df <- NULL
    k <- from
    n <- nrow(lmms)
    while (n > 1) {  # The loop executes until n = 0 or 1.
        ay <- lmms$y[1]  # allele y
        ax <- lmms$x[1]  # allele x
        marks <- rep(FALSE, times = n)
        if (paired.rows) {  # There are always only one/two rows present: (x, y) and/or (y, x)
            marks[1] <- TRUE  # The first row of the current lmms is always chosen.
            marks <- marks | (lmms$x == ay & lmms$y == ax)  # The other one may or may not be found.
        } else {  # Rows may be (x, y), (x, y), (y, x), (x, y), etc.
            marks <- (lmms$x == ax & lmms$y == ay) | (lmms$x == ay & lmms$y == ax)
        }
        allele.pair <- lmms[marks, ]  # extract rows from the data frame lmms: either one or two rows
        m <- nrow(allele.pair)  # m = 1 or 2. Comparing to matrices, lmms[1, ] remains a data frame.
        allele.pair$pair <- rep(k, times = m)
        new.df <- rbind.data.frame(new.df, allele.pair, stringsAsFactors = FALSE)  # push selected into the stack
        if (m == n) {  # m = n = 2, no rows remain: reaches the bottom of iterations
            n <- 0
        } else {
            lmms <- lmms[!marks, ]  # chop off assessed rows from the data frame
            n <- nrow(lmms)
            k <- k + 1
        }
    }
    if (n == 1) {  # nrow(lmms) = 1: reaches the bottom of the above loops
        lmms$pair <- k
        new.df <- rbind.data.frame(new.df, lmms, stringsAsFactors = FALSE)  # push the last row into the stack
    }

    # restore names of the first two columns where rename = TRUE
    if (rename) {
        names(new.df)[1 : 2] <- names.prev
    }

    return(new.df)
}
