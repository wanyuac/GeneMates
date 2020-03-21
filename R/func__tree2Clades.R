#' @title Decompose a tree into a presence-absence matrix of samples in each clade of the tree
#'
#' @description This function is an adaptation of the tree2patterns function in the BugWAS package for phylix.
#'
#' @param tr A phylo object of an input tree.
#' @param sample.order A character vector of sample names to reorder tips of the input tree. In practice, the names
#' follow the order in the SNP matrix, the allelic PAM and the projection matrix C.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' clade.pam <- treeToClades(tr = assoc[["C"]][["tr"]], sample.order = rownames(assoc[["C"]][["C"]]))
#'
#' @author Yu Wan (\email{wanyuac@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 18 July 2017, the lastest edition: 27 July 2017

tree2Clades <- function(tr, sample.order) {
    print(paste0(Sys.time(), ": converting the input tree into a PAM of samples in each clade."))

    require(ape)

    sample.n <- length(tr$tip.label)
    nodes.in <- (sample.n + 1) : (sample.n + tr$Nnode)  # indices of internal nodes, which start from the tip number + 1
    clade.ids <-paste0("N", nodes.in)  # Each clade is named with the index of the internal node that it connects to.

    # initialise a PAM for the presence/absence of each sample in each clade
    # tr$Nnode: number of internal nodes. For a rooted tree, the clade of root has all samples. Thus it is not informative.
    clade.pam <- matrix(0, nrow = sample.n, ncol = tr$Nnode)
    rownames(clade.pam) <- sample.order
    colnames(clade.pam) <- clade.ids

    # fill the PAM with presence/absence status of every sample in each clade
    for (i in nodes.in) {
        clade <- extract.clade(phy = tr, node = i)  # extract a subtree rooting on the selected internal node
        samples.in <- match(clade$tip.label, sample.order)  # locate rows corresponding to samples in this subtree
        clade.pam[samples.in, paste0("N", i)] <- 1  # set corresponding cells to one
    }

    # measure the size of every clade and put results into a named vector
    clade.sizes <- apply(clade.pam, 2, function(x) as.integer(sum(x)))

    return(list(pam = clade.pam, sizes = clade.sizes))
}
