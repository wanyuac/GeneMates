#' @title Firth's penalised logistic regression
#'
#' @description Fit pairwise allelic presence/absence data with Firth's penalised
#' logistic regression without a control for bacterial population structure.
#' Because this function is not the focus of GeneMates, it takes as inputs the
#' outputs of the lmm or findPhysLink function. As a result, users must run
#' either function before calling this one.
#'
#' Dependency: parallel, logistf and data.table
#'
#' @param pat An uncentred matrix of patterns of the allelic presence/absence status
#'      Expect every cell of this matrix to be a dichotomous variable.
#' @param tests.pat A data frame specifying pairs of patterns whose association
#' will be tested for with the logistic regression.
#' @param tests.allele A data frame recording pairs of alleles whose association
#' have been tested for by the findPhysLink function.
#' @param n.cores Number of cores used to run GEMMA in parallel where possible.
#'      -1: automatically detect the number of available cores N, but use N - 1
#' cores (recommended)
#'      0: automatically detect the number of available cores and use all of them.
#' Be careful when the current R session is not running through SLURM.
#'      >= 1: use the number of cores as specified. n.cores is reset to the maximum
#' number of available cores N when n.cores > N.
#' @param p.adj.method Method for correcting p-values. Refer to the base function
#' p.adjust for legitimate values of this argument.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' flr <- penalisedLogisticReg(pat = assoc$alleles$B,
#'        tests.pat = assoc$lmms.pat$dif$h1[, c("y_pat", "x_pat")],
#'        tests.allele = assoc$tests$dif$tests,
#'        n.cores = 8, p.adj.method = "bonferroni")
#'
#' @author (Copyright) Yu Wan (\email{wanyuac@gmail.com})
#' @export
#
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 17 June 2017, last edition: 23/5/2018

penalisedLogisticReg <- function(pat, tests.pat, tests.allele, n.cores = -1,
                                 p.adj.method = "bonferroni") {
    require(parallel)
    require(data.table)  # rbindlist function

    # perform pairwise Firth's logistic regression
    nc <- .setCoreNum(n.cores = n.cores, cores.avai = detectCores())
    print(paste("Use", nc, "cores to fit logistic models.", sep = " "))
    print(paste0(Sys.time(), ": Starting penalised logistic regression."))
    cl <- makeCluster(nc)
    clusterExport(cl = cl,
                  varlist = list("pat", "tests.pat"),
                  envir = environment())  # make variables accessible to different cores
    flms <- parLapply(cl, 1 : nrow(tests.pat), .flr)
    stopCluster(cl)
    print(paste0(Sys.time(), ": Regression analysis finished."))

    # concatenate results and adjust p-values for multiple tests
    print("Concatenating results and adjusting p-values.")
    flms.pat <- as.data.frame(rbindlist(flms))  # pattern-level results
    flms.pat$p_chisq_adj <- p.adjust(p = flms.pat$p_chisq, method = p.adj.method)

    # utilise the existing function .patternToAlleles to restore allele names from patterns
    flms.a <- .patternToAlleles(lmms = flms.pat,
                                tests = list(tests = tests.allele,
                                             y.pats = unique(tests.pat$y_pat)),
                                h1 = TRUE, mapping = NULL)

    # assign pair IDs
    flms.a <- assignPairID(lmms = flms.a)
    print(paste0(Sys.time(), ": This job is finished successfully."))

    return(list(pat = flms.pat, allele = flms.a))  # pattern-level and allele-level results
}

# This is a subordinate function of penalisedLogisticReg.
# Requires data frames pat, tests and mapping in its parental environment.
# i: row number in the data frame mapping
.flr <- function(i) {
    require(logistf)

    t <- tests.pat[i, ]  # take a single row from the data frame
    y <- pat[, paste0("pat_", t[["y_pat"]])]
    x <- pat[, paste0("pat_", t[["x_pat"]])]

    # fit a Firth simple logistic model: logit(Y) ~ X
    r <- logistf(y ~ 1 + x)

    # obtain estimates of parameters and summary statistics; cf. the logistf manual for parsing results of the logistf function
    p.chisq <- r$prob[["x"]]
    d <- data.frame(y_pat = t[["y_pat"]], x_pat = t[["x_pat"]],
                    n_y = sum(y), n_x = sum(x), n_xy = sum(as.logical(y) & as.logical(x)),
                    beta = r$coefficients[["x"]],
                    se = sqrt(r$var[2, 2]),
                    lower_95 = r$ci.lower[["x"]],
                    upper_95 = r$ci.upper[["x"]],
                    chisq = qchisq(p = p.chisq, df = r$df, lower.tail = FALSE),
                    p_chisq = p.chisq, stringsAsFactors = FALSE)

    return(d)
}
