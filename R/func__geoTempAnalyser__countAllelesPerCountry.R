#' @title Count alleles per country
#'
#' @param alleles A character vector of allele names.
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#' @param mapping A data frame in the output of lmm or findPhysLink.
#' @param apam Allelic presence-absence matrix.
#' @param nul A value for zero counts. It can be zero definitely, but setting it to NA makes it easier to draw a heat map.
#' @param freq A logical argument determining if the allele frequency (in %)
#' rather than allele count is returned.
#'
#' @examples
#' a_nat <- countAllelesPerCountry(alleles = nwk$V$allele, sam = sam, mapping = assoc$mapping, apam = assoc$alleles$A, nul = -30)
#'
#' @return This function returns a list of two elements: count and mapping, when
#' freq = FALSE, else, returns a list of three elements: freq, mapping and sample_size_country.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export countAllelesPerCountry
#'
# Copyright 2018 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# First edition: 16 Aug 2018; the latest update: 1 Nov 2018

countAllelesPerCountry <- function(alleles = NULL, sam, mapping, apam, nul = NA,
                                   freq = FALSE) {
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

    # Get country information
    sam <- sam[!is.na(sam$Country), ]  # drop strains whose country information is missing
    countries <- unique(sam$Country)
    if (freq) {  # count number of strains per country
        sample_size_country <- table(sam$Country)  # returns a named vector of integers
    } else {
        sample_size_country <- NULL
    }

    # Initialise the outcome
    cc <- matrix(nul, nrow = ncol(apam), ncol = length(countries),
                 dimnames = list(alleles, countries))  # allele count through time. Rows: sorted allele names; columns: years

    # For each allele, count its occurrence in the population in each year
    for (c in countries) {  # given each country
        strains <- sam$Strain[sam$Country == c]
        n <- length(strains)  # There must be at least one country in the list.
        if (n > 1) {
            apam_c <- apam[strains, ]
            mk <- as.logical(apply(apam_c, 2, function(x) any(x > 0)))  # markers used for dropping empty columns
            n_alleles <- sum(mk)  # number of alleles left
            if (n_alleles > 1) {
                apam_c <- apam_c[, mk]
                for (a in colnames(apam_c)) {
                    num <- sum(apam_c[, a])  # count the total number of this allele
                    cc[a, c] <- ifelse(freq, round(num / sample_size_country[[c]] * 100, digits = 2), num)
                }
            } else if (n_alleles == 1) {
                a <- colnames(apam_c)[mk]  # save the allele name
                num <- sum(apam_c[, mk])  # sum of elements in a column vector
                cc[a, c] <- ifelse(freq, round(num / sample_size_country[[c]] * 100, digits = 2), num)
            } else {
                print(paste("No strain in", c, "had any target gene.", sep = " "))
            }
        } else if (n == 1) {
            apam_c <- apam[strains, ]  # returns a named character vector of a single row
            alleles <- names(apam_c)[as.logical(apam_c)]
            if (length(alleles) > 0) {
                for (a in alleles) {
                    cc[a, c] <- ifelse(freq, round(n / sample_size_country[[c]] * 100, digits = 2), n)
                }
            } else {
                print(paste("The only strain from", c, "did not have any target genes.", sep = " "))
            }
        } else {  # skip this column in ac if no strain was isolated in the current year
            print(paste0("Warning: no strain was found in ", c, "."))
        }
    }

    # Produce the output
    if (freq) {
        out <- list(freq = cc, mapping = mapping, sample_size = sample_size_country)
    } else {
        out <- list(count = cc, mapping = mapping)
    }

    return(out)
}
