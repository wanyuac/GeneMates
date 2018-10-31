#' @title Allelic diversity per country or any kind of regions
#'
#' @description Calculate Simpson's or Shannon's diversity index for alleles per
#' country or any other kind of geographic regions. Each strain must belong to
#' only a single country.
#'
#' @param sam A data frame whose first and second columns are strain and country
#' names, respectively.
#' @param pam A binary allelic presence-absence matrix, where row names are
#' strain names and column names are allele names.
#' @param alleles A character vector specifying alleles to be analysed. Keep NULL
#' to include all alleles in pam.
#' @param countries A character vector specifying countries to be considered. Keep
#' NULL to include all countries in pam.
#' @param method A character argument of either "simpson" or "shannon".
#' @param shannon_base The base for Shannon's diversity index. It is not used for
#' Simpson's diveristy index. Default: exp(1).
#'
#' @return A list of four elements:
#'   c: a matrix of n(countries) by n(alleles) for allele counts per country;
#'   d: a named numeric vector of diversity indices.
#'   p: the final allelic presence-absence matrix on which the diversity indices
#' are computed.
#'   s: isolation information of the final set of strains.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 31 Oct 2018

calcAllelicDiveristyPerCountry <- function(sam, pam, alleles = NULL, countries = NULL,
                                           method = "simpson", shannon_base = exp(1)) {
    require(vegan)

    # Preparation
    sam <- sam[, c(1, 2)]
    names(sam) <- c("Strain", "Country")
    sam <- subset(sam, !is.na(Country))

    if (!is.null(countries)) {
        sam <- subset(sam, Country %in% countries)  # to reduce the volumn of data
    }
    strains <- intersect(x = rownames(pam), y = sam$Strain)  # the final set of strains. Some strains may be removed from the original pam due to absence of any alleles.
    pam <- pam[strains, ]  # Unicity of strain names is guaranteed by the previous function intersect.
    pam <- pam[, as.logical(colSums(pam))]  # remove empty colums after excluding some strains
    sam <- subset(sam, Strain %in% strains)
    countries <- unique(sam$Country)  # the final set of countries to be analysed

    if (!is.null(alleles)) {
        pam <- pam[, intersect(x = alleles, y = colnames(pam))]
    }
    alleles <- colnames(pam)  # the final set of alleles to be analysed

    # Calculation
    m <- matrix(0, nrow = length(countries), ncol = length(alleles),
                dimnames = list(countries, alleles))
    for (k in countries) {
        strains_k <- sam$Strain[sam$Country == k]
        pam_k <- pam[strains_k, ]
        if (length(strains_k) > 1) {
            m[k, ] <- as.numeric(colSums(pam_k))
        } else {  # The length must equal one, which means a single strain was isolated in the current country.
            m[k, ] <- as.numeric(pam_k)  # pam_k is a vector in this case.
        }
    }

    if (method != "simpson") {
        d <- diversity(x = m, index = "shannon", MARGIN = 1, base = shannon_base)
    } else {
        d <- diversity(x = m, index = "simpson", MARGIN = 1)
    }

    return(list(d = d, c = m, p = pam, s = sam))
}
