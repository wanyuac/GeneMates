#' @title Estimating correlation between an allelic distribution pattern and each
#' projection vector.
#'
#' @description This function computes a point-biserial correlation coefficient
#' between each pattern and each projection vector. It works in a similar way as
#' the function corCladeProj. The point-biserial correlation coefficient is a
#' special case of the well-known Pearson correlation coefficient.
#'
#' @param pat.mat A uncentred pattern matrix, which is the element alleles$B in
#' the output of the function lmm or findPhysLink.
#' @param proj.mat A projection matrix, which is the element struc$C$C in the
#' output of the function lmm or findPhysLink.
#'
#' @examples a <- lmm(...)
#' r <- corPatternProj(pat.mat = a$alleles$B, proj.mat = a$struc$C$C)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the latest edition: 25 May 2018

corPatternProj <- function(pat.mat, proj.mat) {
    # Initialisation
    r0 <- data.frame(pattern = integer(0), axis = integer(0), cor = numeric(0),
                     count = integer(0))
    r <- r0  # the final output
    patterns <- colnames(pat.mat)  # pat_1, pat_2, ...
    projs <- colnames(proj.mat)  # c1, c2, ...
    pattern.ids <- sapply(patterns, function(x) as.integer(substr(x, start = 5,
                                                                  stop = nchar(x))))  # a named vector of integers
    proj.ids <- sapply(projs, function(x) as.integer(substr(x, start = 2,
                                                            stop = nchar(x))))

    # Sample count per pattern
    counts <- colSums(pat.mat)  # a named vector of integers

    # Calculating the correlation coefficients
    for (i in patterns) {
        p <- pat.mat[, i]  # a binary vector
        p.id <- pattern.ids[[i]]
        r.i <- r0  # to store rows of each i
        for (j in projs) {
            c <- proj.mat[, j]  # a continuous vector
            r.ij <- data.frame(pattern = p.id,
                               axis = proj.ids[[j]],
                               cor = round(cor(x = p, y = c, method = "pearson"),
                                           digits = 8),
                               count = counts[[i]])
            r.i <- rbind.data.frame(r.i, r.ij)
        }
        r.i <- r.i[order(abs(r.i$cor), decreasing = TRUE), ]  # sort rows by the absolute values of correlation coefficients under each pattern
        r <- rbind.data.frame(r, r.i)
    }

    return(r)
}
