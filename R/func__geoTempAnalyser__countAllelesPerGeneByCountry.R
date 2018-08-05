#' @title Count each gene per country using their alleles
#'
#' @description This function cares any alleles.
#'
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#' @param mapping A data frame in the output of lmm or findPhysLink.
#' @param gm Gene content matrix. Its column names must match the gene names in the mapping data frame.
#' @param nul A value for zero counts. It can be zero definitely, but setting it to NA makes it easier to draw a heat map.
#'
#' @examples cc <- countAllelesPerGeneByCountry(sam = sam, mapping = assoc$mapping, gm = assoc$genes$all, nul = -35)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First version: 2 Aug 2018; the latest edition: 3 Aug 2018

countAllelesPerGeneByCountry <- function(sam, mapping, gm, nul = NA) {
    # Prepare the table of allele counts
    ga <- unique(mapping[, c("class", "gene")])
    ga <- ga[order(ga$class, ga$gene, decreasing = FALSE), ]  # Y axis on the heat map

    # Get country information
    sam <- sam[!is.na(sam$Country), ]  # drop strains whose country information is missing
    countries <- unique(sam$Country)

    # Construct an output matrix
    cc <- matrix(nul, nrow = nrow(ga), ncol = length(countries),
                 dimnames = list(ga$gene, countries))  # allele count through time. Rows: genes sorted by classes; columns: years

    # Populate the matrix
    gm <- gm[, ga$gene]  # an error arises when there are any column names not in ga$gene
    for (c in countries) {  # given each country
        strains <- sam$Strain[sam$Country == c]
        n <- length(strains)  # There must be at least one country in the list.
        if (n > 1) {
            gm_c <- gm[strains, ]
            mk <- as.logical(apply(gm_c, 2, function(x) any(x != "-")))  # markers used for dropping empty columns
            ng <- sum(mk)  # number of genes left
            if (ng > 1) {
                gm_c <- gm_c[, mk]
                for (g in colnames(gm_c)) {
                    cc[g, c] <- sum(as.character(gm_c[, g]) != "-")  # count the total number of whatever alleles
                }
            } else if (ng == 1) {
                g <- colnames(gm_c)[mk]
                gm_c <- gm_c[, mk]  # a named character vector
                cc[g, c] <- sum(as.character(gm_c) != "-")
            } else {
                print(paste("No strain in", c, "had any target gene.", sep = " "))
            }
        } else if (n == 1) {
            gm_c <- gm[strains, ]  # a named character vector
            genes <- names(gm_c)[as.character(gm_c) != "-"]
            if (length(genes) > 0) {
                for (g in genes) {
                    cc[g, c] <- 1
                }
            } else {
                print(paste("The only strain from", c, "did not have any target genes.", sep = " "))
            }
        } else {  # skip this column in ac if no strain was isolated in the current year
            print(paste0("Warning: no strain was found in ", c, "."))
        }
    }

    # Simplify row names by dropping class information from gene names
    ga$gene_sim <- sapply(ga$gene, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1])
    rownames(cc) <- ga$gene_sim[match(rownames(cc), ga$gene)]

    return(list(count = cc, mapping = ga))
}
