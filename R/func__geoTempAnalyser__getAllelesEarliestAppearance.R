#' @title Get earliest appearance of each allele in a vector
#'
#' @description For a group of alleles, this function returns the earliest
#' occurrence time and country of each allele.
#'
#' @param alleles A vector of allele names.
#' @param pam A presence-absence matrix of these alleles. It may come from the element
#' pam in the output list of the function ringPlotPAM.
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#'
#' @examples emg_1 <- getAllelesEarliestAppearance(alleles = targets, pam = rp_1$pam, sam = sam)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#'
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 4 Aug 2018
#'

getAllelesEarliestAppearance <- function(alleles, pam, sam) {
    strains <- rownames(pam)
    emg <- data.frame(Allele = character(0), Strain = character(0),
                      Year_up = integer(0), Country = character(0),
                      stringsAsFactors = FALSE)
    for (a in alleles) {
        sa <- strains[which(pam[, a] > 0)]  # presence of this allele in strains; pam[, a] returns a vector of integers here
        sam_s <- subset(sam, Strain %in% sa)
        year_earliest <- min(sam_s$Year_up)
        i <- which(sam_s$Year_up == year_earliest)
        if (length(i) == 1) {
            emg <- rbind.data.frame(emg, data.frame(Allele = a,
                                                    Strain = sam_s$Strain[i],
                                                    Year_up = year_earliest,
                                                    Country = sam_s$Country[i],
                                                    stringsAsFactors = FALSE))
        } else {
            emg <- rbind.data.frame(emg,
                                    data.frame(Allele = a,
                                               Strain = paste(sam_s$Strain[i], collapse = ","),
                                               Year_up = year_earliest,
                                               Country = paste(sam_s$Country[i], collapse = ","),
                                               stringsAsFactors = FALSE))
        }
    }
    emg <- emg[order(emg$Year_up, decreasing = FALSE), ]

    return(emg)
}
