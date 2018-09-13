#' @title Calculate sign congruence of fixed effects (betas) between linear mixed
#' models (LMMs) and penalised logistic models (PLMs)
#'
#' @description This function calculates the percentage of matched fixed effects
#' having the same sign. The outcome can be visualised as a scatter plot or so.
#'
#' @param lmms A data frame from the function lmm for association status between
#' patterns. Assuming assoc = lmm(...), then lmms = assoc$lmms.pat$dif$h1.
#' @param plms A data frame from the function plr for association status between
#' patterns. Assuming assoc = plr(...), then plms = assoc$pat. Pattern orders for
#' association tests must match to those in lmms. Otherwise, both data frames
#' cannot merge correctly.
#' @param p.max Upper bound of raw p-values for congruence assessment.
#' @param step Step size of betas for calculating the congruence percentage.
#' @param nevLog A logical argument determining whether -log10 will be applied to
#' p-values in the result.
#'
#' @examples sg <- betaSignCongruence(lmms = assoc_lmm$lmms.pat$dif$h1,
#' plms = assoc_lgr$pat, p.max = 1, step = 0.005/nrow(assoc_lmm$lmms.pat$dif$h1),
#' nevLog = TRUE)
#'
#' @return A data frame of percentages as measurements of sign congruence.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 13 Sep 2018, the lastest edition: 13 Sep 2018

betaSignCongruence <- function(lmms, plms, p.max = 1, step = 0.0005, nevLog = TRUE) {
    # Prepare the target data frame
    b <- merge(x = lmms[, c("y_pat", "x_pat", "beta", "p_wald")],
               y = plms[, c("y_pat", "x_pat", "beta", "p_chisq")],
               by = c("y_pat", "x_pat"), all = TRUE, sort = FALSE)
    names(b)[c(3, 5)] <- c("beta_lmm", "beta_plm")
    b <- b[, c("y_pat", "x_pat", "beta_lmm", "beta_plm", "p_wald", "p_chisq")]

    # Calculate proportion of betas sharing the same signs between PLMs and LMMs
    if (step <= 0) {
        step <- 0.0005
    }
    p0_set <- seq(from = step, to = p.max, by = step)  # upper bound of p-values for significance based on Bonferroni correction

    # Define the output data frame
    # P0: upper bound of p-values for significance. For instance, P0 = 0.05.
    # N: number of significant associations shared between PLMs and LMMs
    # Nc: number of betas showing the same signs
    # Perc: percentage of rows showing sign congruence
    c <- data.frame(P0 = numeric(0), N = integer(0), Nc = integer(0),
                    Perc = numeric(0))

    for (p0 in p0_set) {
        bp <- subset(b, p_wald <= p0 & p_chisq <= p0)  # shared associations in both kinds of models
        m <- nrow(bp)
        if (m > 0) {
            nc <- sum((sign(bp$beta_lmm) * sign(bp$beta_plm)) == 1)  # number of rows showing identical signs of the two betas
            c <- rbind.data.frame(c, data.frame(P0 = p0, N = m, Nc = nc,
                                                Perc = round(nc / m * 100, digits = 4)))
        }  # else, do nothing
    }

    if (nevLog) {
        c$P0_nevLg <- round(-log10(c$P0), digits = 4)
        c <- c[, c("P0", "P0_nevLg", "N", "Nc", "Perc")]
    }

    return(c)
}
