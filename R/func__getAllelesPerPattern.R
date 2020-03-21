#' @title Converting an allele-pattern data frame into a pattern-alleles data frame
#'
#' @description This functions reformats the data frame "alle.pat" or mapping[, c("allele", "pattern")] in outputs of the findPhysLink function.
#' It is useful when a user wants to list all alleles under each particular pattern.
#'
#' For example, given an input data frame:
#'     allele | pattern; a1 | 1; a2 | 1
#'
#' The function returns a data frame containing a row:
#'     pattern | alleles; 1 | a1,a2
#'
#' @param mapping the data frame alle.pat or mapping in outputs of the findPhysLink function. It can be other data frames, but must
#' contain names "allele" and "pattern".
#' @param sep the delimiter for alleles that belong to the same pattern in the output. It is a comma by default.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' ap1 <- getAllelesPerPattern(assoc[["mapping"]], sep = ",")
#' ap2 <- getAllelesPerPattern(assoc[["alleles"]][["alle.pat"]], sep = ",")  # alternatively
#'
#' @author Yu Wan (\email{wanyuac@126.com})
#' @export
#
# Copyright 2017 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# First and the lastest edit: 19 Sep 2017

getAllelesPerPattern <- function(mapping, sep = ",") {
    # sanity check
    if (any(!(c("allele", "pattern") %in% names(mapping)))) {
        stop("Column name error: the data frame 'mapping' must contain two columns named 'allele' and 'pattern' respectively.")
    }

    # pool alleles under the same pattern
    pats <- sort(unique(mapping$pattern), decreasing = FALSE)
    d <- data.frame(pattern = pats,
                    alleles = sapply(pats, function(p) paste(mapping$allele[which(mapping$pattern == p)], collapse = sep)),
                    stringsAsFactors = FALSE)

    return(d)
}
