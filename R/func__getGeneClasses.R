#' @title Extracting gene classes from SRST2-formatted allele IDs.
#'
#' @description SRST2 defines each allele ID using several fields:
#' [sequence cluster ID]__[cluster symbol (the name of a gene, product class or group)]__[gene or allele name]__[sequence ID] [additional annotations separated by ";"].
#' For instance, a valid sequence ID is:
#' 225__SulI_Sul__SulI__1616 no;no;SulI;Sul;U37105;4069-4908;840.
#' Usually, the cluster symbol is the antimicrobial class for an antimicrobial
#' resistance gene. This function aims to extract this symbol from the allele ID.
#' It is a simple function, however, this package includes it because of its
#' frequent usage.
#'
#' @param ids A character vector of allele IDs.
#'
#' @examples
#' assoc <- lmm(...)
#' arg_cls <- data.frame(allele = assoc$mapping$allele, class = getGeneClass(x = assoc$mapping$gene), stringsAsFactors = FALSE)
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export getGeneClass
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 18 Apr 2018; the lastest edition: 23 July 2018

.getClassName <- function(g) {
    # This is a subordinate function of getGeneClass.
    # g: gene ID.
    c <- strsplit(g, "_", fixed = TRUE)[[1]]
    c <- c[length(c)]  # take the last element

    return(c)
}

getGeneClass <- function(ids) {
    cls <- as.character(sapply(ids, .getClassName))

    return(cls)
}
