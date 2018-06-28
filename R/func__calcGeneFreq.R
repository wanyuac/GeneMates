#' @title Calculate gene frequency
#'
#' @description This function takes as input the output data frame of countAlleles (more specifically, only uses its
#' second and third columns) and computes the frequency per gene given a two-column data frame of gene IDs and counts
#' across a bacterial population.
#'
#' @param x a two-column data frame where the first column stores gene names and the second column stores counts.
#'     Columns must not be factor variables.
#' @param n total number of bacterial isolates/strains in a population
#' @param ord a logical parameter specifying whether the output data frame should be sorted by the gene frequency
#'
#' @examples
#' mapping <- countAlleles(...)
#' gf <- calcGeneFreq(x = mapping[, c("gene", "count")], n = 1125)  # assuming there are 1125 strains
#'
#' @author Yu Wan (\email{wanyuac@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 26 May 2017; the latest edition: 24 Jan 2018


calcGeneFreq <- function(x, n, ord = TRUE) {
    names(x) <- c("gene", "count")
    genes <- unique(x$gene)
    y <- data.frame(gene = character(0), count = integer(0), freq = numeric(0), n_a = integer(0), stringsAsFactors = FALSE)  # n_a: number of alleles
    for (g in genes) {
        r <- subset(x, gene == g)  # There is at least one row being selected.
        c <- sum(r$count)  # total number of times observing the current gene
        y.new <- data.frame(gene = g, count = c, freq = round(c / n * 100, digits = 2), n_a = nrow(r), stringsAsFactors = FALSE)
        y <- rbind.data.frame(y, y.new, stringsAsFactors = FALSE)
    }
    if (ord) {
        y <- y[order(y$freq, decreasing = TRUE), ]
    }
	rownames(y) <- NULL

    return(y)
}
