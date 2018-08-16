#' @title Count alleles per year
#'
#' @param alleles A character vector of allele names.
#' @param sam A data frame of sample information. Required columns: Strain, Country, Year_low and Year_up.
#' @param mapping A data frame in the output of lmm or findPhysLink.
#' @param apam Allelic presence-absence matrix.
#' @param nul_count A value for zero counts. It can be zero definitely, but setting it to NA makes it easier to draw a heat map.
#' @param nul_freq A value for zero frequency.
#' @param combine A named list defining periods when counts of each allele are
#' added up to make a new count. Each element in this list is an integer vector
#' of years. For example, combine = list("1985-2005" = c(1985, 1986, 1988, 2005)).
#' The vector must not contain a single element.
#'
#' @examples
#' a_yr <- countAllelesPerYear(alleles = nwk$V$allele, sam = sam, mapping = assoc$mapping,
#'                             apam = assoc$alleles$A, nul = -30)
#'
#' @return A list of three elements: count, freq, n (number of strains per year),
#' w (weighted counts per strain) and mapping. When combine is a list, another
#' two elements called comb_count, comb_freq and period are included.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export countAllelesPerYear
#'
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First version: 3 Aug 2018; the latest edition: 16 Aug 2018

countAllelesPerYear <- function(alleles = NULL, sam, mapping, apam, nul_count = NA,
                                nul_freq = NA, combine = NULL) {
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
    ws <- round(1 / year_range, digits = 6)  # weight for allelic presence per year: the wider the range, the smaller the weight
    names(ws) <- sam$Strain
    ac <- matrix(nul_count, nrow = ncol(apam), ncol = length(years),
                 dimnames = list(alleles, as.character(years)))  # allele count through time. Rows: sorted allele names; columns: years
    af <- matrix(nul_freq, nrow = ncol(apam), ncol = length(years),
                 dimnames = list(alleles, as.character(years)))
    n_yr <- data.frame(Year = integer(0), Count = integer(0), Count_w = numeric(0))  # weighted sample size per year

    # For each allele, count its occurrence in the population in each year
    for (yr in years) {
        strains <- sam$Strain[(sam$Year_low <= yr) & (sam$Year_up >= yr)]
        n <- length(strains)
        if (n > 1) {
            apam_yr <- apam[strains, ]
            mk <- as.logical(apply(apam_yr, 2, function(x) any(x > 0)))  # for removal of empty columns
            n_alleles <- sum(mk)
            n_yr_w <- sum(as.numeric(ws[strains]))  # weighted strain number of this year
            n_yr <- rbind.data.frame(n_yr, data.frame(Year = yr, Count = n, Count_w = n_yr_w))
            if (n_alleles > 1) {
                apam_yr <- apam_yr[, mk]
                for (a in colnames(apam_yr)) {
                    strains_a <- strains[as.logical(apam_yr[, a])]
                    n_a_yr_w <- sum(as.numeric(ws[strains_a]))  # weighted allele count
                    ac[a, as.character(yr)] <- n_a_yr_w
                    af[a, as.character(yr)] <- round(n_a_yr_w / n_yr_w * 100, digits = 2)
                }
            } else if (n_alleles == 1) {
                a <- colnames(apam_yr)[mk]
                apam_yr <- apam_yr[, mk]  # becomes an integer vector
                strains_a <- strains[as.logical(apam_yr)]
                n_a_yr_w <- sum(as.numeric(ws[strains_a]))  # weighted allele count
                ac[a, as.character(yr)] <- n_a_yr_w
                af[a, as.character(yr)] <- round(n_a_yr_w / n_yr_w * 100, digits = 2)
            } else {
                print(paste("None of strains collected in", yr, "had any target allele detected.", sep = " "))
            }
        } else if (n == 1) {
            apam_yr <- apam[strains, ]  # a named row vector of integers
            alleles_yr <- names(apam_yr)[as.logical(apam_yr)]
            n_yr_w <- as.numeric(ws[[strains]])
            n_yr <- rbind.data.frame(n_yr, data.frame(Year = yr, Count = n, Count_w = n_yr_w))
            if (length(alleles_yr) > 0) {
                for (a in alleles_yr) {
                    ac[a, as.character(yr)] <- n_yr_w
                    af[a, as.character(yr)] <- 100.00
                }
            } else {
                print(paste("The only strain collected in", yr, "had no target allele detected.", sep = " "))
            }
        } else {  # n = 0
            print(paste0("There was no strain collected in ", yr, "."))
            n_yr <- rbind.data.frame(n_yr, data.frame(Year = yr, Count = 0, Count_w = 0))
        }
    }

    # Combine columns
    if (is.list(combine)) {
        mat_comb <- .combineCounts(ac = ac, af = af, comb = combine, n_yr = n_yr,
                                   nul_count = nul_count, nul_freq = nul_freq)
        out <- list(count = ac, freq = af, n = n_yr, mapping = mapping, w = ws,
                    count_comb = mat_comb[["count"]], freq_comb = mat_comb[["freq"]],
                    period = mat_comb[["period"]])
    } else {
        out <- list(count = ac, freq = af, n = n_yr, mapping = mapping, w = ws)
    }

    return(out)
}

.combineCounts <- function(ac, af, comb, n_yr, nul_count = 0, nul_freq = 0) {
    periods <- names(comb)
    cc <- matrix(nul_count, nrow = nrow(ac), ncol = length(periods),
                 dimnames = list(rownames(ac), periods))
    cf <- matrix(nul_freq, nrow = nrow(af), ncol = length(periods),
                 dimnames = list(rownames(af), periods))  # combined frequencies
    p_sum <- data.frame(Period = character(0), Count = integer(0), Count_w = numeric(0),
                        stringsAsFactors = FALSE)
    yrs_excl <- integer(0)  # years to exclude

    for (p in periods) {
        yrs <- comb[[p]]
        yrs_excl <- union(yrs, yrs_excl)
        if (length(yrs) == 1) {
            stop("Argument error: the year vector must not contain only a single element.")
        }
        ac_p <- ac[, as.character(yrs)]
        ac_p <- .replaceVals(m = ac_p, top = 0, new_val = 0)  # Otherwise, the combined counts are wrong if nul_count < 0.
        n_p <- 0
        n_p_w <- 0

        for (yr in yrs) {
            i <- which(n_yr$Year == yr)
            n_p <- n_p + n_yr$Count[i]
            n_p_w <- n_p_w + n_yr$Count_w[i]
        }
        cc[, p] <- rowSums(ac_p)
        cf[, p] <- round(cc[, p] / n_p_w * 100, digits = 2)
        p_sum <- rbind.data.frame(p_sum, data.frame(Period = p, Count = n_p,
                                                    Count_w = n_p_w,
                                                    stringsAsFactors = FALSE),
                                  stringsAsFactors = FALSE)
    }
    cc <- .replaceVals(m = cc, top = 1e-7, new_val = nul_count)  # Six digits are preserved for weights.
    cf <- .replaceVals(m = cf, top = 0.001, new_val = nul_freq)   # Since only two ditigs are preserved, 0.001 is an unrealistic positive number.

    c_ac <- ac[, !(colnames(ac) %in% as.character(yrs_excl))]
    f_ac <- af[, !(colnames(af) %in% as.character(yrs_excl))]
    c_ac <- c_ac[, order(as.integer(colnames(c_ac)), decreasing = FALSE)]
    f_ac <- f_ac[, order(as.integer(colnames(f_ac)), decreasing = FALSE)]

    out <- list(count = cbind.data.frame(cc, c_ac),
                freq = cbind.data.frame(cf, f_ac),
                period = p_sum)
}

.replaceVals <- function(m, top = 0, new_val = 0) {
    for (i in 1 : nrow(m)) {
        for (j in 1 : ncol(m)) {
            x <- m[i, j]
            m[i, j] <- ifelse(x < top, new_val, x)
        }
    }

    return(m)
}
