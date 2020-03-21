#' @title Find out rows corresponding to alleles of a given pair of genes
#'
#' @description This is a little tool function for getting all rows relating to
#' all alleles of two given genes.
#'
#' @param x A data frame to be subset for the result.
#' @param map A data frame to be searched against, which maps allele names to gene
#' information. There must be two columns named "allele" and "gene".
#' @param g1 A string for the first gene name.
#' @param g2 A string for the second gene name.
#' @param col1 Name or index for the first column of alleles in x. Default: y.
#' @param col2 Name or index for the second column of alleles in x. Default: x.
#'
#' @examples getAssocGenePair(x = assoc$mapping, g1 = "StrA_AGly", g2 = "DfrA8_Tmt", col1 = "y", col2 = "x")
#'
#' @return A data frame from x.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#'
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the lastest edition: 27 August 2018

getAssocGenePair <- function(x, map, g1, g2, col1 = "y", col2 = "x") {
    a_x <- map$allele[map$gene == g1]
    a_y <- map$allele[map$gene == g2]
    a <- c(a_x, a_y)
    out <- x[(x[, col1] %in% a) & (x[, col2] %in% a), ]

    return(out)
}
