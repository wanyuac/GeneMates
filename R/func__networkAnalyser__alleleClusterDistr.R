#' @title Compute a matrix for presence of allele clusters based on allelic presence-absence status
#'
#' @description This function generates a matrix for co-occurrence of all alleles
#' of each cluster in a give collection of bacterial samples. This function is
#' particularly useful in investigating the distribution of allele clusters (such
#' as maximal cliques) in bacterial samples. It can be considered as a simplified
#' version of function getClusterMemberCoocurrence.
#'
#' @param cls A data frame of two columns for cluster IDs (first column) and allele
#' contents (the second column, comma-delimited). This data frame can be generated
#' from the output of function summariseCliques. Each cluster must consist of at
#' least two alleles.
#' @param pam An uncentred presence-absence matrix for alleles in the data frame cls.
#' This is a binary matrix, where 1 denotes presence and 0 denotes absence.
#'
#' @return An n-by-m matrix, where n is the sample size and m is the cluster number.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2019 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First and the latest edition: 23 Feb 2019

alleleClusterDistr <- function(cls, pam) {
    names(cls) <- c("ID", "Alleles")
    C <- matrix(0, nrow = nrow(pam), ncol = nrow(cls),
                dimnames = list(rownames(pam), cls$ID))  # the output: a co-occurrence network
    for (i in 1 : nrow(cls)) {
        r <- cls[i, ]
        id <- r[["ID"]]
        alleles <- strsplit(x = r[["Alleles"]], split = ",", fixed = TRUE)[[1]]  # assume length(alleles) > 1
        if (length(alleles) > 1) {  # valid situation
            A <- pam[, alleles]
            C[, id] <- !apply(A, MARGIN = 1, function(r) any(r == 0))  # 1: presence of the cluster; 0, otherwise
        } else {
            print(paste("Error in cluster definition: cluster", id, "only has one allele.",
                        sep = " "))
            C[, id] <- rep(NA, times = nrow(C))
        }
    }

    return(C)
}
