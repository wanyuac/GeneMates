#' @title Evaluate evidence of physical linkage for positive associations with summary statistics of distance measurements
#'
#' @description This function must be used after running through the findPhysLink function.
#'
#' @param lmms a list produced by .summariseDist, which contains LMM parameters
#' and distance measurements
#' @param min.beta only associations with beta's >= min.beta will be considered
#' as showing evidence of physical linkage.
#' @param max.p only associations with P values <= max.p will be considered as signficant.
#' @param max.range An integer specifying the maximum range (in bp) of in-group
#' allelic physical distances allowed to call consistent.
#' @param min.pIBD (optional) Minimum probability of the root of a minimum
#' inclusive clade to display a positive binary trait, such as having a pair of
#' alleles co-occurring or a specific allelic physical distance. Default: 0.9 (90\%).
#'
#' @author Yu Wan (\email{wanyuac@gmail.com})
#' @export
#
#  Copyright 2017-2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 1 June 2017; latest edition: 2 April 2018

evalPL <- function(lmms, min.beta = 0, max.p = 0.05, max.range = 2000, min.pIBD = 0.9) {
    print(paste0(Sys.time(), ": Evaluating evidence of physical linkage."))

    # Sanity check
    if (class(lmms) != "list") {  # The input must be a list produced by the function summariseDist.
        stop("Argument error: \"lmms\" must be a list produced by the function summariseDist.")
    } else if (!("dif" %in% names(lmms))) {  # two requisite elements
        stop("Structural error: the list \"lmms\" must contain an elements called \"dif\".")
    }

    # Concatenate data frame elements in the list lmms into a single one
    assoc <- NULL
    for (g in names(lmms)) {  # c("dif", "idd") or "dif" only
        d <- lmms[[g]]
        if (!is.null(d)) {  # assign groups of association results
            if (g == "dif") {
                d$dif <- rep(1, times = nrow(d))  # Are alleles differently distributed? Yes
            } else {
                d$dif <- rep(0, times = nrow(d))  # identically distributed; # beta = 1 for all idd. alleles and all associations are considered as significant.
            }
        }
        assoc <- rbind.data.frame(assoc, d, stringsAsFactors = FALSE)
    }

    # Score evidence for physical linkage
    assoc$s_a <- as.integer(apply(as.matrix(assoc[, c("beta", "p_adj")]), 1,
                                  .scoreAssocEvidence, min.beta, max.p))
    assoc$s_d <- as.integer(apply(assoc[, c("d_in_n", "d_in_range", "pIBD_in")], 1,
                                  .scoreDistEvidence, max.range, min.pIBD))
    assoc$w_d <- round(assoc$m_in * assoc$s_d, digits = 6)  # weighted s_d. Measurability has six decimals.
    assoc$score <- assoc$s_a + assoc$w_d  # The score belongs to [-2, 2] or equals NA.

    return(assoc)
}

# Scoring evidence for the presence of physical linkage based on associations
# This is a subordinate function for evalPL.
# scores for physical linkage based on the support of physical distances
#   1: supports the presence of linkage
#   0: cannot determine or lack of evidence
#   -1: evidence against linkage
# Notice p-values are NA's for intra-pattern LMM results.
.scoreAssocEvidence <- function(x, min.beta, max.p) {  # x: a single row from the matrix of statistics
    # direction of the current association based on the size of this fixed effect
    # The ifelse function returns NA when beta = NA, which is a desirable behaviour.
    s <- ifelse(x[["beta"]] > min.beta, 1, -1)  # s = 1, -1 or NA

    # significance of the current association
    # significant: beta != 0, insignificant: cannot reject the hypothesis that beta = 0.
    # reset the score to zero as the beta may actually equal zero (cannot reject the null hypothesis)
    # Notice NA values are kept as NA * n = NA and NA * NA = NA.
    s <- s * as.integer(x[["p_adj"]] <= max.p)  # s = 0 when x$p_adj > max.p and s = s otherwise

    return(s)  # return 1, -1 or NA
}

# Scoring evidence for the presence of physical linkage based on physical distances
# This is also a subordinate function for evalPL.
# scores for physical linkage based on the support of physical distances
#   1: supports the presence of linkage
#   0: cannot determine or lack of evidence
#   -1: evidence against linkage
.scoreDistEvidence <- function(x, max.range, min.pIBD) {
    n.d <- x[["d_in_n"]]
    not.IBD <- x[["pIBD_in"]] < min.pIBD
    if (n.d > 1) {
        if (x[["d_in_range"]] <= max.range) {  # The distances are similar to each other.
            # For the next logical statement, if it returns TRUE, we know the
            # most-recent common ancestor of the minimum inclusive clade of the
            # strains where the in-group distances are measured may not display
            # a similar distance, that is, the genomic structure is unlikely to
            # be identical by descendent for these strains; if it returns FALSE,
            # it is unsurprising to see similar distances present in closely
            # related strains (hence s = 0, indicating insufficient evident and
            # the necessity of having more observations).
            s <- ifelse(not.IBD, 1, 0)
        } else {  # In-group distances are too divergent.
            # Unlikely to be IBD: not surprising to see diverse distances;
            # IBD: evidence against physical linkage (hence s = -1).
            s <- ifelse(not.IBD, 0, -1)
        }
    } else {  # n.d = 1
        s <- 0  # insufficient evidence to make an inference
    }

    return(s)
}
