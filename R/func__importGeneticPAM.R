#' @title Import the genetic presence/absence (p/a) matrix
#'
#' @description This function reads the genetic presence-absence matrix and returns a list of relevant matrices.
#' It automatically matches isolates to those in the SNP matrix.
#'
#' @param pam Expect pam to be "modified_allele_matrix.txt" produced by the cdhitFormatter pipeline.
#'     First column of pam consists of isolate names. Other columns: presence/absence profile of every gene
#'     The genetic PAM may contain more alleles than the allelic PAM does, because the latter applies a more stringent criterion to filter data.
#' @param pam.delim a single character for the delimiter in the text file of the genetic PAM
#' @param outliers a vector of isolate names to be excluded from rows of the PAM
#' @param min.count the minimal number of times that a gene occurs in all isolates excluding outliers
#'     By default, only genes that do not occur at all are removed (min.count = 1).
#'     When min.count = 2, isolate-specific genes are removed as well.
#' @param genes.rm a vector of gene names to be excluded from columns of the genetic PAM.
#' @param sample.order a vector of isolate names for the PAM to be matched with.
#' It can be used to filter the genetic PAM for in-group isolates.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2017 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# An early version: 21 May 2017; the latest edition: 22 Aug 2018

importGeneticPAM <- function(pam, pam.delim = "\t", outliers = NULL, min.count = 1,
                             genes.rm = NULL, sample.order = NULL) {
    pamData.class <- class(pam)
    if (pamData.class == "character") {
        print(paste0("Reading the genetic presence-absence matrix from the file ", pam))
        pam <- read.delim(pam, header = TRUE, sep = pam.delim, check.names = FALSE, stringsAsFactors = FALSE)
        isolates <- pam[, 1]
        pam <- as.matrix(pam[, -1])  # chop out the first column, which contains isolate names, and make a genuine genetic PAM
        rownames(pam) <- isolates

        # filter the genetic PAM for ingroup and outlier isolates
        if (!is.null(outliers)) {
            pam <- pam[!(isolates %in% outliers), ]  # exclude profiles of outlier isolates (the "which" function does not work well)
            isolates <- rownames(pam)  # refresh the list of isolate names
        }

        # match the order of isolates to that of the core-genome SNP table
        # sum(isolates %in% sample.order) < length(isolates): some isolates in PAM are not found in the SNP matrix
        # sum(sample.order %in% isolates) < length(sample.order): some isolates in the SNP matrix are not found in the PAM
        # Correct condition: all isolates in the PAM are present in the SNP matrix, and vice versa.
        if (class(sample.order) == "character") {
            # The above condition is a weak check for this argument. So please do not challenge it through providing unexpected values.
            # Assuming a character vector of isolate names are provided.
            if (any(! sample.order %in% isolates)) {
                # Abnormal situation: the required sample set contains samples that are absent in the allelic PAM
                stop("Error: some samples in the required sample set are not found in the genetic PAM.")
            } else if (any(! isolates %in% sample.order)) {
                print("Warning: sample names in the genetic matrix is a subset of the required sample set.")  # Two sets of names must be the same.
                print("Therefore, the genetic matrix will be trimmed for the required sample set.")
            }

            pam <- pam[sample.order, ]  # An error arises when sample.order has extra strains comparing to pam.
        }

        # first, remove empty columns (with all values equalling "-") after excluding outliers
        gene.counts <- as.integer(apply(pam, 2, function(x) sum(x != "-")))  # number of isolates harbouring each gene
        genes.incl <- gene.counts > 0
        pam <- pam[, genes.incl]
        genes <- colnames(pam)  # genes having at least one allele
        gene.counts <- gene.counts[genes.incl]
        pam.all <- pam  # copy the PAM with empty columns removed, which may be useful for users

        # then, remove genes that do not occur at enough times when specified
        if (min.count <= 0) {
            min.count <- 1  # turn off the filter for gene frequency
        }
        if (min.count > 1) {
            genes.incl <- gene.counts >= min.count  # observing each gene (through any alleles of it) at least twice among in-group isolates (no signs of horizontal gene transfer if the gene occurs only once)
            if (sum(genes.incl) < length(genes)) {  # when there are some genes to exclude
                genes.excl <- genes[!genes.incl]  # alleles to be omitted due to insufficient abundance
                pam <- pam[, genes.incl]  # drop some columns
                genes <- colnames(pam)  # refresh the gene list
                print(paste0(length(genes.excl), " genes are excluded from the presence-absence matrix due to insufficient occurrence counts: "))
                print(paste(genes.excl, collapse = ","))
            } else {
                print("No gene is excluded due to insufficient occurrence.")
            }
        }

        # remove genes from the PAM in accordance with genes.rm (a user-specified list of genes to be excluded)
        if (is.character(genes.rm)) {
            print(paste0("Removing genes according to the list: ", paste(genes.rm, collapse = ","), "."))
            flag <- genes %in% genes.rm  # flag rows where gene is present in the set genes.excl
            if (sum(flag) > 0) {  # when some genes are on the removal list
                pam <- pam[, !flag]  # filter out alleles of these genes. This is a subset of the original PAM.
            }
        }

        # construct the return list
        pam <- list(pam = pam, all = pam.all, alleles.inc = .getAlleleNames(pam))  # pam[["pam"]] is an isolate-by-gene matrix.
    } else if (pamData.class == "list") {
        if (sum(names(pam) %in% c("pam", "all", "alleles.inc")) == 3) {
            print("Skip importing the genetic presence-absence data as they are already provided.")
        } else {
            stop("Input error of pam: it is a list but does not have all elements.")
        }
    } else {
        stop("Input error: requiring either a list or a path to the raw matrix file.")
    }

    return(pam)
}
