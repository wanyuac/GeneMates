#' @title Count each gene per year using their alleles
#'
#' @description Neither function cares any alleles.
#'
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#' @param mapping A data frame in the output of lmm or findPhysLink.
#' @param gm Gene content matrix. Its column names must match the gene names in the mapping data frame.
#' @param nul A value for zero counts. It can be zero definitely, but setting it to NA makes it easier to draw a heat map.
#'
#' @examples
#' ac <- countAllelesPerGeneByYear(sam = sam, mapping = assoc$mapping, gm = assoc$genes$all, nul = -30)
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#'
# Copyright 2018 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# First version: 2 Aug 2018; the latest edition: 3 Aug 2018

countAllelesPerGeneByYear <- function(sam, mapping, gm, nul = NA) {
    # Prepare the table of allele counts
    ga <- unique(mapping[, c("class", "gene")])
    ga <- ga[order(ga$class, ga$gene, decreasing = FALSE), ]  # Y axis on the heat map

    # Get the year information
    ys <- min(sam$Year_low) : max(sam$Year_up)

    # Construct an output matrix
    ac <- matrix(nul, nrow = nrow(ga), ncol = length(ys),
                 dimnames = list(ga$gene, as.character(ys)))  # allele count through time. Rows: genes sorted by classes; columns: years

    # Create a named vector for weights of allelic presence for all strains
    year_range <- sam$Year_up - sam$Year_low + 1
    ws <- round(1 / year_range, digits = 8)  # weight for allelic presence per year: the wider the range, the smaller the weight
    names(ws) <- sam$Strain

    # Populate the matrix
    gm <- gm[, ga$gene]  # an error arises when there are any column names not in ga$gene
    for (y in ys) {  # given each year
        strains <- sam$Strain[(sam$Year_low <= y) & (sam$Year_up >= y)]
        n <- length(strains)
        if (n > 1) {
            gm_y <- gm[strains, ]  # a sub-matrix
            gy <- as.logical(apply(gm_y, 2, function(c) any(c != "-")))  # to drop empty columns
            ng <- sum(gy)  # number of genes left
            if (ng > 1) {  # multiple genes left
                gm_y <- gm_y[, gy]  # get a sub-matrix from gm_y
                genes_y <- colnames(gm_y)
                for (g in genes_y) {
                    strains_g <- strains[as.character(gm_y[, g]) != "-"]  # strains having this gene in this year
                    ac[g, as.character(y)] <- sum(as.numeric(ws[strains_g]))  # weighted allele count for the current gene in this year
                }
            } else if (ng == 1) {  # only one gene left
                g <- colnames(gm_y)[gy]
                gm_y <- gm_y[, gy]  # returns a named character vector across strains
                strains_g <- names(gm_y)[as.character(gm_y) != "-"]
                ac[g, as.character(y)] <- sum(as.numeric(ws[strains_g]))
            } else {
                print(paste("None of strains collected in", y, "had any target genes detected.", sep = " "))
            }
        } else if (n == 1) {
            gm_y <- gm[strains, ]  # a named character vector when there is only a single row selected
            genes_y <- names(gm_y)[as.character(gm_y) != "-"]
            if (length(genes_y) > 0) {
                c <- ws[[strains]]
                for (g in genes_y) {
                    ac[g, as.character(y)] <- c  # There is only a single strain.
                }
            } else {
                print(paste("The only strain collected in", y, "had no target gene detected.", sep = " "))
            }
        } else {  # skip this column in ac if no strain was isolated in the current year
            print(paste0("No strain was collected in ", y, "."))
        }
    }

    # Simplify row names by dropping class information from gene names
    ga$gene_sim <- sapply(ga$gene, function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][1])
    rownames(ac) <- ga$gene_sim[match(rownames(ac), ga$gene)]

    return(list(count = ac, mapping = ga))
}
