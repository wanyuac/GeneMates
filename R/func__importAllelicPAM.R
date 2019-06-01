#' @title Import the allelic presence-absence matrix (PAM) of bacterial genes.
#'
#' @description This function reads the allelic p/a matrix and returns a list of relevant matrices.
#'
#' @param pam either the final PAM or the file name of a raw PAM created by the PAMmaker pipeline.
#'     First column: Sample; other columns: allele identifiers.
#' @param pam.delim a single character for the delimiter in the text file of the allelic PAM.
#' @param outliers a vector of isolate names to be excluded from rows of the allelic PAM.
#' @param min.count the minimal number of times that an allele occurs in all isolates excluding outliers.
#'     By default, only alleles that do not occur at all are removed (min.count = 1).
#'     When min.count = 2, isolate-specific alleles are removed as well.
#' @param alleles.inc a vector of allele names to be included for columns of the allelic PAM.
#' @param output.y path to the output file for the centred pattern matrix Y.
#'     Keep output.y = "" to prevent the function from printing the matrix.
#'     In this case, only a list is returned to the user.
#' @param sample.order a vector of isolate names for the PAM to be matched with.
#' It can be used to filter the allelic PAM for in-group isolates.
#' @param skip wheter to skip overwriting existant output files.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# An early version: 21 May 2017; the latest edition: 22 Aug 2018

importAllelicPAM <- function(pam, pam.delim = "\t", outliers = NULL, min.count = 1,
                             alleles.inc = NULL, output.y = "", sample.order = NULL,
                             skip = TRUE) {
    pamData.class <- class(pam)
    if (pamData.class == "character") {
        print(paste0("Reading the allelic presence-absence matrix from the file ", pam))
        pam <- read.delim(pam, header = TRUE, sep = pam.delim, check.names = FALSE, stringsAsFactors = FALSE)
        isolates <- pam[, 1]
        pam <- as.matrix(pam[, -1])  # chop out the first column, which contains isolate names, and make a genuine p/a matrix
        rownames(pam) <- isolates

        # exclude profiles of outlier isolates
        if (!is.null(outliers)) {
            pam <- pam[!(isolates %in% outliers), ]
            isolates <- rownames(pam)  # refresh the list of isolate names
        }

        # match isolates in the SNP matrix and the allelic PAM
        if (class(sample.order) == "character") {  # Otherwise, isolates will not be selected and rearranged in pam.
            # The above condition is a weak check for this argument. So please do not challenge it through providing unexpected values.
            # Assuming a character vector of isolate names are provided.
            if (any(! sample.order %in% isolates)) {
                # Abnormal situation: the required sample set contains samples that are absent in the allelic PAM
                stop("Error: some samples in the required sample set are not found in the allelic PAM.")
            } else if (any(! isolates %in% sample.order)) {
                print("Warning: sample names in the allelic matrix is a subset of the required sample set.")  # Two sets of names must be the same.
                print("Therefore, the allelic matrix will be trimmed for the required sample set.")
            }

            pam <- pam[sample.order, ]  # match the order of isolates to that of the core-genome SNP table
        }

        # remove empty columns
        allele.counts <- as.integer(apply(pam, 2, sum))
        inc <- allele.counts > 0  # identify empty columns and flag them
        pam <- pam[, inc]
        allele.counts <- allele.counts[inc]
        alleles <- colnames(pam)
        pam.all <- pam  # copy the original PAM with empty column removed; not filtered for allele counts yet.

        # remove alleles that do not occur at a minimal times
        if (min.count <= 0) {
            min.count <- 1  # only remove empty columns (with all values equalling zero)
        }
        if (min.count > 1) {
            inc <- as.integer(allele.counts) >= min.count  # mark alleles of sufficient abundance
            if (sum(inc) < length(alleles)) {  # Some alleles do not have sufficient abundance.
                alleles.excl <- alleles[!inc]  # alleles to be omitted due to insufficient abundance
                pam <- pam[, inc]  # Every allele must present in at least two isolates.
                alleles <- alleles[inc]
                print(paste0(length(alleles.excl), " alleles are excluded from the presence-absence matrix due to insufficient counts:"))
                print(paste(alleles.excl, collapse = ","))
            } else {
                print("No gene is excluded due to an insufficient count.")
            }
        }

        # remove alleles of genes that have been excluded in accordance with the alleles.inc argument
        if (is.character(alleles.inc)) {
            print("Filtering alleles according to the argument alleles.inc.")
            inc <- colnames(pam) %in% alleles.inc
            alleles.excl <- alleles[!inc]
            pam <- pam[, inc]
            print(paste0(length(alleles.excl), " alleles are removed according to the removal list alleles.inc:"))
            print(paste(alleles.excl, collapse = ","))
        }

        # prepare genotype data for GEMMA: both X and Y are pattern matrices
        a.pat <- .getPatterns(pam)  # compress pam into an uncentred pattern matrix
        Y <- scale(a.pat[["B"]], center = TRUE, scale = FALSE)  # a centred pattern matrix subjecting to association analysis
        X <- t(Y)  # samples go to columns to make it compatible with the BIMBAM format for genotypes
        rownames(X) <- NULL  # A bimbam-formatted genotype matrix does not take row names into account.

        # print the BIMBAM-format "phenotype" file (In fact, Y refers to genotypes)
        # write a centred pattern matrix into a tab-delimited file
        if (output.y != "") {
            .writeData(x = Y, output = output.y, message = "[Matrix Y]", skip = skip)
        }

        # wrap results into a list; B: an uncentred pattern matrix as a reduction of columns of the allelic p/a matrix A
        # B is the raw pattern matrix. dim(B) = dim(Y), dim(B) = dim(t(X)). Y is the column-centred version of B. X = t(Y).
        pam <- list(Y = Y, X = X, A = pam, B = a.pat[["B"]], alle.pat = a.pat[["alle.pat"]],
                    pat.sizes = a.pat[["pat.sizes"]], all = pam.all)
    } else if (pamData.class == "list") {
        # This branch works when this function is called within the findPhysLink function and
        # a user has imported the allelic PAM outside of the findPhysLink function already.
        if (sum(names(pam) %in% c("Y", "X", "A", "B", "alle.pat", "pat.sizes", "all")) == 7) {  # check for completeness
            print("Skip importing the allelic presence-absence data as they are already provided.")
            if (output.y != "") {  # print the BIMBAM-compatible "phenotype" file (In fact, Y refers to genotypes)
                .writeData(x = pam[["Y"]], output = output.y, message = "[Matrix Y]", skip = skip)
            }
        } else {
            stop("Input error of pam: it is a list but does not have all elements.")
        }
    } else {
        stop("Input error: requiring either a list or a path to the raw matrix file as an input.")
    }

    return(pam)
}
