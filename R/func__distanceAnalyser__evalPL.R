#' @title Evaluate evidence of physical linkage for positive associations with summary statistics of distance measurements
#'
#' @description This function evaluates evidence of physical linkage using association
#' status and consistency in allelic physical distances.
#'
#' @param lmms This argument can be a list produced by the function summariseDist,
#' which contains LMM parameters and distance measurements. It can also be the
#' data frame assoc in the output of findPhysLink or lmms in the output of function
#' lmm. When the latter two kinds of input are used, users can rescore the evidence
#' using new criteria.
#' @param min.beta only associations with beta's >= min.beta will be considered
#' as showing evidence of physical linkage.
#' @param max.p only associations with P values <= max.p will be considered as signficant.
#' @param max.range An integer specifying the maximum range (in bp) of in-group
#' allelic physical distances allowed to call consistent.
#' @param min.pIBD (optional) Minimum probability of the root of a minimum
#' inclusive clade to display a positive binary trait, such as having a pair of
#' alleles co-occurring or a specific allelic physical distance. Default: 0.9 (90\%).
#' @param score.dist (optional) A logical parameter specifying whether allelic
#' physical distances should be taken into account for scoring edges. Default: FALSE.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2017-2019 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 1 June 2017; latest edition: 20 April 2019

evalPL <- function(lmms, min.beta = 0, max.p = 0.05, max.range = 2000, min.pIBD = 0.9,
                   score.dist = FALSE) {
    print(paste0(Sys.time(), ": Evaluating evidence of physical linkage."))

    # Sanity check
    arg_cls <- class(lmms)
    if (arg_cls == "list") {  # lmms = summariseDist(...)
        if (!("dif" %in% names(lmms))) {
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
    } else if (arg_cls == "data.frame") {  # element "assoc" from findPhysLink
        assoc <- lmms
    } else {
        stop("Argument error: \"lmms\" must be a list produced by the function summariseDist.")
    }

    # Score evidence for physical linkage
    assoc$s_a <- as.integer(apply(as.matrix(assoc[, c("beta", "p_adj")]), 1,
                                  .scoreAssocEvidence, min.beta, max.p))
    if (score.dist) {
        assoc$c <- as.integer(apply(assoc[, c("d_in_n", "d_in_range", "pIBD_in")], 1,
                                    .scoreDistEvidence, max.range, min.pIBD))  # consistency score
        assoc$s_d <- round(assoc$m_in * assoc$c, digits = 6)  # distance score (measurability-weighted consistency score). The measurability has six decimals.
        assoc$s <- assoc$s_a + assoc$s_d  # The score belongs to [-2, 2] or equals NA.
    } else {
        assoc$m <- rep(NA, times = nrow(assoc))  # I choose NA rather than zero in this circumstance since we do not know the measurability when the distance measurement is not carried out.
        assoc$m_in <- assoc$m
        assoc$s <- assoc$s_a
        assoc$s_d <- assoc$m_in
        assoc$c <- assoc$s_d
    }

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
