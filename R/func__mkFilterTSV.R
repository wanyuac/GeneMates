#' @title Making a guidance TSV file as an input for the physDist pipeline
#'
#' @description The function prints allele names of each strain when it has at least two alleles.
#'
#' @param allele.mat A matrix in the output format of SRST2, which shows an allele per gene for each sample.
#' This argument overrides pam.a. Row names of this matrix are sample names and column names are gene names.
#' There must not be a column of sample names.
#' @param pam.a A presence-absence matrix of alleles. It is an alternative to the allele.mat argument and it
#' must not be a data frame. Row names of this matrix are sample names and column names are allele names. Again,
#' there must not be a column of sample names.
#' @param output Name of the output file
#'
#' Either allele.mat or pam.a is expected to be used as an input for this function.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# First edition: 23 Janury 2017, 21 May 2017; the latest edition: 18 July 2018.
# Copyright 2017-2018 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0

mkFilterTSV <- function(allele.mat = NULL, pam.a = NULL, output = "targeted_isolates_alleles.tsv") {
    # get input
    if (!is.null(allele.mat)) {
        if (is.data.frame(allele.mat)) {
            allele.mat <- as.matrix(allele.mat, stringsAsFactors = FALSE)
        }
        sample.alleles <- apply(allele.mat, 1, function(r) r[r != "-"])  # a list of character vectors
    } else if (!is.null(pam.a)) {
        if (is.data.frame(pam.a)) {  # assumes a data frame with row names and no sample name column
            pam.a <- as.matrix(pam.a)  # converts it into a matrix of zero and one while keeping its row names
        }
        alleles <- colnames(pam.a)  # Allele names are retrieved from column names.
        sample.alleles <- apply(pam.a, 1, function(r) alleles[as.logical(r)])  # returns a list of character vectors per sample
    } else {
        stop("Error: alleles and pam.a must not be both NULL.")
    }

    # remove samples that do not have any allele calls
    len <- as.logical(sapply(sample.alleles, length))  # len[i] = FALSE when length(sample.alleles[[i]]) = 0
    sample.alleles <- sample.alleles[len]

    # write the output file line-by-line
    if (length(sample.alleles) > 0) {  # when there are some elements left
        cat(NULL, file = output)  # erase previous outputs where present
        for (sample in names(sample.alleles)) {  # go through every sample
            allele.calls <- sample.alleles[[sample]]
            if (length(allele.calls) > 1) {  # To measure distances, we need at least two alleles present in each strain.
                write(paste(sample, paste(allele.calls, collapse = ","), sep = "\t"),
                      file = output, append = TRUE)
            }  # otherwise, skip this sample
        }
    }

    return(sample.alleles)
}
