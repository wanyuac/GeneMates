#' @title Count alleles per year
#'
#' @param alleles A character vector of allele names.
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#' @param mapping A data frame in the output of lmm or findPhysLink.
#' @param apam Allelic presence-absence matrix.
#' @param nul A value for zero counts. It can be zero definitely, but setting it to NA makes it easier to draw a heat map.
#'
#' @examples
#' a_yr <- countAllelesPerYear(alleles = nwk$V$allele, sam = sam, mapping = assoc$mapping,
#'                             apam = assoc$alleles$A, nul = -30)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First version: 3 Aug 2018; the latest edition: 12 Aug 2018

countAllelesPerYear <- function(alleles = NULL, sam, mapping, apam, nul = NA) {
    # For each allele, count its occurrence in all strains in each year.
    # am: alleleic presence-absence matrix, which only contains binary values

    # Use all alleles in the mapping table when the vector of alleles is unspecified.
    if (!is.null(alleles)) {
        mapping <- subset(mapping, allele %in% alleles)
    }
    mapping <- mapping[order(mapping$class, mapping$gene, mapping$allele, decreasing = FALSE), ]  # The main purpose for using mapping.
    alleles <- mapping$allele  # get a vector of reordered alleles, assuming there is no duplicate in mapping$allele

    # Subset the matrix for these alleles
    apam <- apam[, alleles]
    if (length(alleles) == 1) {
        apam <- matrix(as.integer(apam), ncol = 1, dimnames = list(names(apam), alleles))  # restore the vector back to a column matrix
    }

    # Drop empty columns
    apam <- apam[, as.logical(apply(apam, 2, function(x) any(x > 0)))]
    if (ncol(apam) == 0) {
        stop("Warning: there is no allele present in any strains.")  # an abnormal situation
    }

    # Initialise the outcome
    years <- min(sam$Year_low) : max(sam$Year_up)
    year_range <- sam$Year_up - sam$Year_low + 1
    ws <- round(1 / year_range, digits = 8)  # weight for allelic presence per year: the wider the range, the smaller the weight
    names(ws) <- sam$Strain
    ac <- matrix(nul, nrow = ncol(apam), ncol = length(years),
                 dimnames = list(alleles, as.character(years)))  # allele count through time. Rows: sorted allele names; columns: years

    # For each allele, count its occurrence in the population in each year
    for (yr in years) {
        strains <- sam$Strain[(sam$Year_low <= yr) & (sam$Year_up >= yr)]
        n <- length(strains)
        if (n > 1) {
            apam_yr <- apam[strains, ]
            mk <- as.logical(apply(apam_yr, 2, function(x) any(x > 0)))  # for removal of empty columns
            n_alleles <- sum(mk)
            if (n_alleles > 1) {
                apam_yr <- apam_yr[, mk]
                for (a in colnames(apam_yr)) {
                    strains_a <- strains[as.logical(apam_yr[, a])]
                    ac[a, as.character(yr)] <- sum(as.numeric(ws[strains_a]))
                }
            } else if (n_alleles == 1) {
                a <- colnames(apam_yr)[mk]
                apam_yr <- apam_yr[, mk]  # becomes an integer vector
                strains_a <- strains[as.logical(apam_yr)]
                ac[a, as.character(yr)] <- sum(as.numeric(ws[strains_a]))
            } else {
                print(paste("None of strains collected in", yr, "had any target allele detected.", sep = " "))
            }
        } else if (n == 1) {
            apam_yr <- apam[strains, ]  # a named row vector of integers
            alleles_yr <- names(apam_yr)[as.logical(apam_yr)]
            if (length(alleles_yr) > 0) {
                c <- ws[[strains]]
                for (a in alleles_yr) {
                    ac[a, as.character(yr)] <- c
                }
            } else {
                print(paste("The only strain collected in", yr, "had no target allele detected.", sep = " "))
            }
        } else {
            print(paste0("There was no strain collected in ", yr, "."))
        }
    }

    return(list(count = ac, mapping = mapping))
}
