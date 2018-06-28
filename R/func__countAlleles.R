#' @title Mapping alleles to genes
#'
#' @description Make a data frame to summarise genes, alleles and allele frequencies.
#' @param pam.genes a compiled matrix showing allele calls across isolates, where extended allele IDs have been assigned
#' @param pam.alleles a presence-absence matrix of alleles, which is derived from the matrix pam.genes.
#' @param patterns (optional) a data frame consisting of two columns: allele and pattern ID
#' @param pat.sizes (optional) a two-column data frame (pattern and size) for the number of alleles per pattern
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First version: 26 May 2017; the latest edition: 24 May 2018

countAlleles <- function(pam.genes, pam.alleles, patterns = NULL, pat.sizes = NULL) {
    print("Summarising the allelic presence/absence data.")
    genes <- colnames(pam.genes)  # genes whose alleles (any) are present in at least two isolates
    alleles <- colnames(pam.alleles)  # alleles present in at least two isolates

    # initialise the result
    mapping <- data.frame(allele = character(0), gene = character(0), class= character(0),
                          count = integer(0), stringsAsFactors = FALSE)

    # populate the list
    for (gene in genes) {  # for each gene
        # get all alleles of it on the genetic PAM
        gene.alleles <- pam.genes[, gene]  # take a column in the pam.genes matrix
        gene.alleles <- unique(gene.alleles[gene.alleles != "-"])  # remove empty values, getting all alleles of the current gene
        gene.class <- getGeneClass(ids = gene)

        # However, a few alleles may appear only once and thus are not included in the allelic PAM when the filter for allelic counts has been applied outside of this function (such as in importAllelicPAM).
        for (a in gene.alleles) {  # So we need to check every allele if it is in the allelic PAM.
            if (a %in% alleles) {  # This is a valid allele.
                mapping <- rbind(mapping,
                                 data.frame(allele = a, gene = gene, class = gene.class,
                                            count = sum(pam.alleles[, a]), stringsAsFactors = FALSE))
            } else {  # Some alleles may be excluded from the allele matrix for some reasons.
                cat("Warning: the allele", a, "is not present in the allelic matrix.\n", sep = " ")
            }
        }
    }

    mapping$freq <- round(mapping$count / nrow(pam.alleles) * 100, digits = 2)  # relative frequency in the current collection of isolates; notice nrow(pam.alleles) = nrow(pam.genes)

    # attach optional columns to the output
    if (is.data.frame(patterns)) {
        mapping$pattern <- patterns$pattern[match(mapping$allele, patterns$allele)]  # retrieve pattern IDs, sorted by allele frequencies
        mapping <- mapping[, c("allele", "gene", "class", "count", "freq", "pattern")]  # rearrange columns
        if (is.data.frame(pat.sizes)) {  # in addition, if the pattern sizes are provided
            mapping$pat.size <- pat.sizes$size[match(mapping$pattern, pat.sizes$pattern)]  # retrieve pattern sizes
        }
    }

    # sort by allele counts
    mapping <- mapping[order(mapping$count, decreasing = TRUE), ]  # rearrange alleles in accordance with a decreasing order of counts
    rownames(mapping) <- NULL

    return(mapping)
}
