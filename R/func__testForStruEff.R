#' @title Testing for random structural effects of sample projections on each pattern of allele distributions
#'
#' @description This function tests for effects of bacterial subpopulations on each pattern of allele
#' distributions across the whole population. This kind of effects constitutes a part (the group specific
#' effect vg) of random effects of a linear mixed model (the other part is the environmental effect
#' ve). More specifically, it carries out chi-square tests on the vector of structural effects gamma
#' given observations (a pattern y and a relateness matrix K). Since the phylix package aims to
#' identify physically linked bacterial genes, this function only tests for the structural effects
#' of patterns that have been used for linear mixed models in the function findPhysLink.
#'
#' This function is designed under the assumption that users have run functions findPhysLink and
#' projectSamples.
#'
#' Dependent package: parallel
#'
#' @param pat.h0 A single data frame or a list of data frames of parameter estimates under
#' null linear mixed models. It expects to take as input a concatenation of elements
#' lmms.pat[["dif"/"idd"]]$h0 in outputs of the function findPhysLink.
#' @param Y A column-wisely zero-centred pattern matrix for all haploid samples. An ideal input
#' is the element alleles[["Y"]] on the output list of findPhysLink.
#' @param C A projection matrix of samples on eigenvectors corresponding to positive eigenvalues.
#' It is expected to be the element C on the output list of the function projectSamples. Sample
#' names (row names) of C must be the same as those of Y (row names) and K (row and column names).
#' @param K A relatedness matrix with row and column names. Importantly, row and column names,
#' that is, sample names, must follow the same order as sample names in the Y matrix. The function
#' findPhysLink ensures this match between samples in both matrices. This matrix can be acquired
#' from the element K in the output of the function projectSamples.
#' @param L The number of biallelic cgSNPs used to perform singular-value decomposition of the SNP matrix.
#' Assuming assoc <- findPhysLink(...), L equals ncol(assoc[["snps"]][["G"]]).
#' @param n.cores Number of computational cores that will be used in parallel for this function.
#' It follows the same convention defined in the function findPhysLink. For simplicity, set it to
#' zero to automatically detect and use all available cores; set it to -1 to leave one core out.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' C <- projectSamples(...)
#' s <- testForStruEff(pat.h0 = list(assoc[["lmms.pat"]][["dif"]][["h0"]], assoc[["lmms.pat"]][["idd"]][["h0"]]),
#' Y = assoc[["alleles"]][["Y"]], C = C[["C"]], K = C[["K"]], L = ncol(assoc[["snps"]][["G"]]) n.cores = 8)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 6 July 2017, latest edition: 7 July 2017

testForStruEff <- function(pat.h0, Y, C, K, L, n.cores = 0) {
    #===== check if samples are the same between Y, C and K =====
    samples.y <- rownames(Y)
    samples.C <- rownames(C)
    samples.K <- rownames(K)

    # match samples when they are not matched to each other
    # An error of subscript overflow arises when there are any samples unique to a matrix.
    if (any(samples.C != samples.y)) {
        print("Warning: rearranging rows of C to match Y.")
        C <- C[samples.y, ]  # match rows of C to Y
    }
    if (any(samples.K != samples.y)) {
        print("Warning: rearranging rows and columns of K to match Y.")
        K <- K[samples.y, samples.y]
    }

    #===== concatenate data frames when they come as a list =====
    if (class(pat.h0) == "list") {
        pat.h0 <- .catH0(pat.h0)
    }

    #===== launch a test for each pattern using parallel computing =====
    require(parallel)
    require(data.table)

    # apply a common test to every pattern
    patterns <- sort(pat.h0$y_pat, decreasing = FALSE)
    print(paste0(Sys.time(), ": starting chi-square tests for structural effects underlying each of ",
                 length(patterns), " patterns."))

    cl <- makeCluster(.setCoreNum(n.cores = n.cores, cores.avai = detectCores()))
    clusterExport(cl = cl, varlist = list("pat.h0", "Y", "C", "K", "L"), envir = environment())
    tests <- parLapply(cl, patterns, .testForEffectsPerPattern, pat.h0, Y, C, K, L)
    stopCluster(cl)

    print(paste0(Sys.time(), ": all patterns have been tested."))

    # concatenate results of every pattern
    tests <- rbindlist(tests)  # number of rows = ncol(Y) * ncol(C)
    tests$p_adj <- p.adjust(tests$p, method = "bonferroni")  # correction of p-values for multiple tests

    return(tests)
}

# Concatenate data frames of LMM parameters estimated under H0.
# This is a subordinate function of testStruEff.
# Since parameters are estimated under the null hypothesis, there is only a single variable
# y in every data frame. As such, the following algorithm starts with the first element on
# the list and iteratively append rows from the rest of data frames when the corresponding
# y pattern has not been present in the working data frame yet.
.catH0 <- function(pat.h0) {
    n <- length(pat.h0)
    if (n == 1) {
        pat.h0 <- pat.h0[[1]]
    } else {
        tmp1 <- pat.h0[[1]]
        for (i in 2 : n) {  # go through other data frames
            tmp2 <- pat.h0[[i]]
            m <- nrow(tmp2)
            for (j in 1 : m) {
                pat <- tmp2$y_pat[j]
                if (!(pat %in% tmp1$y_pat)) {
                    tmp1 <- rbind.data.frame(tmp1, tmp2[j, ])
                }
            }
        }
        pat.h0 <- tmp1
    }

    return(pat.h0)
}

# A subordinate function of testForStruEff, testing for structural effects for each pattern
# pat: an integer as a pattern index
.testForEffectsPerPattern <- function(pat, pat.h0, Y, C, K, L) {
    id <- paste0("pat_", as.character(pat))

    # make column vectors
    y <- Y[, id]
    y.len <- length(y)
    y <- matrix(data = y, nrow = y.len, ncol = 1)  # convert y into a column vector
    l <- matrix(data = 1, nrow = y.len, ncol = 1)  # a column vector of 1's

    # retrieve the REML estimate of lambda under the null hypothesis for LMMs
    lambda <- pat.h0$lambda0_remle[which(pat.h0$y_pat == pat)]

    # calculate tau^(-1) from lambda
    H.inv <- solve(lambda * K + diag(y.len))  # H(lambda)^(-1); diag(n) returns a unit matrix of n rows and columns; inv: the inverse of a square matrix;
    tl.Hinv <- crossprod(x = l, y = H.inv)  # t(l) * H.inv
    P <- H.inv - H.inv %*% l %*% solve(tl.Hinv %*% l) %*% tl.Hinv  # In fact, it is the same as H.inv numerically when the intercept term is numerically cancelled (which is the case for this package).
    tau.inv <- as.numeric((crossprod(x = y, y = P) %*% y) / (y.len - 1))  # under the null model

    # calculate posterior mean and variance of (gamma | y, K)
    c.num <- ncol(C)  # the number of informative eigenvectors (of positive eigenvalues)
    M <- solve(crossprod(x = C, y = C) + L / lambda * diag(c.num))
    mu <- as.numeric(M %*% crossprod(x = C, y = y))  # converts a matrix into a numeric vector
    sigma <- tau.inv * M
    gamma.var <- diag(sigma)

    # Chi-square tests for every component for the currrent Y pattern
    chisq <- mu ^ 2 / gamma.var  # construct a chi-square statistic from a normal random variable with its mean and variance
    p <- pchisq(q = chisq, df = 1, lower.tail = FALSE)  # upper-tail p-value of a chi-square test

    # construct an output data frame for the current pair of patterns
    d <- data.frame(y_pat = rep(pat, times = c.num), lambda = rep(lambda, times = c.num), tau_inv = rep(tau.inv, times = c.num),
                    c = 1 : c.num, mu = mu, var = gamma.var, chisq = chisq, p = p)
    d <- d[order(d$p, d$c, decreasing = FALSE), ]  # The most significant association comes first.

    return(d)
}
