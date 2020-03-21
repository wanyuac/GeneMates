#' @title Merging nodes representing identically distributed alleles into a single
#' node in an output network of findPhysLink.
#'
#' @description Merge nodes of identically distributed alleles into a single node to simplify the
#' linkage network. When lmms.only = FALSE, two alleles are mergable only if the physical distances
#' measured between them are well supporting the presence of physical linkage (namely, s_d = 1)
#' and the distributions of alleles are not identical by descent (IBD); when lmms.only = TRUE, two
#' alleles are mergable if they are identically distributed (hence they are represented by the same
#' pattern for linear mixed models).
#'
#' This function takes as input the output of the core function findPhysLink or lmm.
#' Hence it can only be used after running either function. The output is a data frame
#' for network visualisation, which consists of associations between differently distributed
#' alleles and merged clusters. Identically distributed alleles that are not supported by
#' either physical distances or non-IBD status are not included in this data frame. In
#' addition, internal associations (where beta = 1) between identically distributed alleles
#' are not covered by this network either.
#'
#' This function accepts a data frame having unidirected edges as the argument of assoc.
#'
#' @param assoc The data frame named assoc in findPhysLink's output. By design, alleles are all
#' paired in this data frame.
#' @param other.cols Additional columns of the data frame assoc that a user wants to include in the results.
#' Notice variables corresponding to these columns must no be affected by idd. alleles of the same
#' cluster. For example, the beta value does not differ with any alleles in the same cluster. In
#' other words, the estimate of beta for the model y ~ x1 does not differ from y ~ x2 when x1 = x2.
#' By contrast, summary statistics of distance measurements for the distance between y and x1 may
#' differ from those between y and x2. Results can be misleading if this requirement is not satisfied.
#' @param lmms.only A logical argument specifying whether the input data frame "assoc" comes from the
#' function lmm rather than findPhysLink.
#' @param replace.names A logical argument specifying whether to replace cluster names in
#' the assoc data frame with actual allele names. For example, "c1" becomes "allele1&allele2"
#' if replace.names equals TRUE.
#'
#' @examples
#' # Example 1
#' assoc <- findPhysLink(...)
#' a <- mergeIddAlleles(assoc[["assoc"]])
#'
#' # Example 2
#' a <- lmm(...)
#' lmms <- subset(a$dif$h1, p_adj <= 0.05 & beta > 0)
#' lmms$dif <- rep(1, times = nrow(lmms))  # not identically distributed
#' lmms <- rbind.data.frame(lmms,
#'                          cbind.data.frame(a$idd$h1, dif = rep(0, times = nrow(a$idd$h1))),
#'                          stringsAsFactors = FALSE)
#' lmms_merged <- mergeIddAlleles(assoc = lmms, lmms.only = TRUE)
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2017 Yu Wan <wanyuac@126.com>
# Licensed under the Apache License, Version 2.0
# First edition: 30 June 2017; the latest edition: 21 April 2019.

mergeIddAlleles <- function(assoc, other.cols = NULL, lmms.only = FALSE,
                            replace.names = TRUE) {
    # extract linkage information of idd. alleles ==========
    # For identically distributed alleles, n_y = n_x = n_xy.
    if (lmms.only) {
        a <- subset(assoc, dif == 0)[, c("pair", "y", "x")]
        assoc <- subset(assoc, dif == 1)[, c("pair", "y", "x", "y_pat", "x_pat",
                                             "n_y", "n_x", "n_xy", "beta",
                                             "p_adj", "l_remle", other.cols)]
    } else {
        a <- subset(assoc, dif == 0)[, c("pair", "y", "x", "s_d", "c", "d_in_n",
                                         "m_in", "pIBD_in")]  # "Error in 1:nrow(cls) : argument of length 0" arises when this becomes empty.
        assoc <- subset(assoc, dif == 1)[, c("pair", "y", "x", "y_pat", "x_pat",
                                             "n_y", "n_x", "n_xy", "d_in_n", "m_in",
                                             "s", "s_a", "s_d", "c", "beta",
                                             "p_adj", "l_remle", "clade", "size",
                                             "ds_max", "f_xy", "pIBD_in",
                                             other.cols)]
    }
    assoc <- assoc[order(assoc$pair, decreasing = FALSE), ]

    # Identify clusters and merge alleles into nodes
    cls <- .identifyPerfectLinkedCluster(a, lmms.only = lmms.only)  # given dif = 0

    # modify the data frame of associations between differently distributed alleles ==========
    # merge rows under allele names of clusters
    # b$a: a linkage table of differently distributed alleles/clusters
    # b$edges: a list of edges between differently distributed alleles before merging under each cluster IDs
    b <- .modifyAssocTable(assoc, cls[["cls"]], lmms.only, replace.names)  # given dif = 1

    return(list(cls = cls[["cls"]], merged = cls[["merged"]], unmerged = cls[["unmerged"]],
                assoc = b$a, edges = b$edges))
}

# This is a subordinate function of mergeIddAlleles. It only works for identically
# distributed alleles. It starts with an initial vector of allele names and looks
# for all other alleles connecting to these initial alleles. In other words, it
# iteratively finds out all alleles that form a single cluster in perfect linkage.
# This function is able to work on an incomplete association graph for identically
# distributed alleles. Such a graph is produced when there are different but identically
# distributed alleles from the same gene. Since the lmm function does not test for
# associations between alleles of the same gene, an association graph for physically
# linked idd. alleles is not necessarily complete in this scenario.
.identifyPerfectLinkedCluster <- function(a, lmms.only = FALSE) {
    a <- a[order(a$pair, decreasing = FALSE), ]  # sort rows by pair IDs

    # The sum of to.merge must be an even number because idd. alleles form
    # bidirectional associations  ==========
    if (lmms.only) {
        keep.sep <- NULL  # Everything will be merged.
    } else {
        # Both alleles are mergable if their physical linkage is well supported
        # and the distributions of alleles are not due to IBD (c takes this into
        # account).
        to.merge <- a$s_d == 1  # Pairs of idd. alleles that fail this criterion will not appear as a merged node but separate nodes in the output association network.
        keep.sep <- a[!to.merge, ]  # a data frame for un-mergable alleles
        if (nrow(keep.sep) == 0) {
            keep.sep <- NULL
        }
        a <- a[to.merge, ]  # extract mergable edges
    }

    # initialise result variables ==========
    b <- NULL
    cl <- data.frame(ID = character(0), Alleles = character(0), Size = integer(0),
                     stringsAsFactors = FALSE)  # a data frame for cluster-level information
    cl.count <- 0
    cl.id <- ""  # cluster ID
    merged <- NULL  # to be a data frame for mergable alleles

    # Find out as many as alleles until the vector becomes constant ==========
    # Are there other edges from the current alleles to identically distributed alleles?
    while (nrow(a) >= 2) {
        alleles <- a$x[c(1, 2)]  # start from the first two alleles
        delta.n <- 2  # measures the growth of the allele vector; Two is an arbitary initial number for entering the iterations.
        while (delta.n > 0) {  # Extract all idd. alleles of the current pair of alleles using the transitivity of identity
            # identify first neighbours based on outgoing edges
            # This vector (pairs) is non-redundant by nature.
            n <- length(alleles)  # number of previously identified alleles
            pairs <- a$pair[a$x %in% alleles]  # returns at least a single pair a$x[1 : 2]
            b <- subset(a, pair %in% pairs)
            alleles <- unique(b$x)  # extend the vector of alleles; Associations between idd. alleles are always symmetric.
            delta.n <- length(alleles) - n  # becomes zero when all neigbours of the current pair of alleles are extracted
        }
        cl.count <- cl.count + 1
        cl.id <- paste0("c", as.character(cl.count))
        b$cluster <- rep(cl.id, times = nrow(b))
        merged <- rbind.data.frame(merged, b, stringsAsFactors = FALSE)
        cl <- rbind.data.frame(cl, data.frame(ID = cl.id,
                                              Alleles = paste(alleles, collapse = ","),
                                              Size = length(alleles),
                                              stringsAsFactors = FALSE),
                               stringsAsFactors = FALSE)
        a <- subset(a, !(pair %in% pairs))
    }

    if (nrow(cl) > 0) {
        y <- list(cls = cl, merged = merged, unmerged = keep.sep)
    } else {
        y <- list(cls = NULL, merged = NULL, unmerged = keep.sep)
    }

    return(y)
}

# Modify the original table of associations (the element "assoc" in the output
# of findPhysLink) according to identified clusters. This function essentially
# merges vertex properties in a network.
.modifyAssocTable <- function(assoc, cls, lmms.only = FALSE, replace.names = TRUE) {
    # assuming that assoc is a data frame of symmetric associations
    if (!is.null(cls)) {
        edges <- vector(mode = "list", length = 0)
        for (i in 1 : nrow(cls)) {  # for every cluster
            cl.id <- cls$ID[i]  # get the current cluster ID
            alleles <- strsplit(x = cls$Alleles[i], split = ",", fixed = FALSE)[[1]]  # alleles comprising this cluster
            cluster.edges <- (assoc$x %in% alleles) | (assoc$y %in% alleles)  # select edges involving allele clusters

            if (sum(cluster.edges) > 0) {
                a <- subset(assoc, cluster.edges)  # select all edges involving alleles of the current cluster
                b <- subset(assoc, !cluster.edges)  # other edges, which will be processed in the next iteration

                # substitute allele names in a with the current cluster ID
                a$y[a$y %in% alleles] <- cl.id
                a$x[a$x %in% alleles] <- cl.id
                edges[[cl.id]] <- a  # add modified edges into the list

                # replace cluster names with allele names
                if (replace.names) {
                    allele.names <- paste(alleles, collapse = "&")  # paste allele names back into a single string
                    a$y[which(a$y == cl.id)] <- allele.names
                    a$x[which(a$x == cl.id)] <- allele.names
                }

                # merge rows where edges become the same (namely, x1 = x2, y1 = y2)
                # after the allele-name substitution in this iteration
                a <- .mergeIdenticalEdges(a, lmms.only)

                assoc <- rbind.data.frame(a, b, stringsAsFactors = FALSE)
            }  # No modification is performed when this cluster is not connected to any other alleles
        }
        assoc <- assoc[order(assoc$pair, decreasing = FALSE), ]
        print(paste("Reminder: this modified association table does not include",
                    "associations and physical distances between identically",
                    "distributed alleles. For these alleles, plase inspect the",
                    "original association table to infer their co-transfer status.",
                    sep = " "))  # An important message for users.
    } else {  # No edges will be merged.
        print("No edges are merged.")
        edges <- NULL
    }

    return(list(a = assoc, edges = edges))
}

# Aggregate edges (rows) where y1 = y2 and x1 = x2
# This is a subordinate function of .modifyAssocTable, so its behaviour must be
# understood in the context of its parental function.
.mergeIdenticalEdges <- function(a, lmms.only = FALSE) {
    # Merge vertex properties
    merged <- NULL
    while (nrow(a) > 0) {
        x <- a$x[1]
        y <- a$y[1]
        to.merge <- (a$x == x) & (a$y == y)  # a logical vector determining edges of the same direction to be merged
        b <- subset(a, to.merge)  # extract rows corresponding to these pairs; nrow(b) must be even.
        a <- subset(a, !to.merge)  # other rows, which are not processed in this iteration

        # make two rows per group of identical edges (the data frame b)
        # Under the assumption of symmetric directions for each pair of alleles, b must have at least two rows even
        # if there is no duplicated row.
        r <- b[1, ]  # returns a single-row data frame

        if (!lmms.only) {
            # Since idd. alleles of the same cluster are physically linked (s_d = 1),
            # I consider the new link betwwen the allele and the cluster is
            # supported by the distances as long as there is a single link between
            # an allele and any allele of a cluster is supported by physical distances.
            # So I use the max function to get the strongest support for physical
            # linkage based on all c scores relating to this cluster.
            # c = -1, 0 or 1; m_in in [0, 1]
            # To calculate the new measurability, since the cluster of idd. alleles
            # has been merged into a single node, I treat all edges connecting to
            # members of this node as equally contributing to the distance
            # measurements. Therefore, assuming there are n in-bound edges for y ~ x
            # connecting to this cluster, we have n times the co-occurrence events
            # (n_xy is the same for all edges as alleles are idd.), hence the
            # new in-group distance frequency m_in is the sum of original in-group
            # distance frequencies averaged by n.
            r$d_in_n <- sum(b$d_in_n)  # c is a data frame of a single row

            # Notice if a1 and a2 are idd. alleles and a3 ~ a1 or a1 ~ a3 is
            # not significant, then edges a3 ~ a2 or a2 ~ a3 does not exist as well.
            # Therefore, it is impossible to see a3 ~ a1, a3 ~ a2, a2 ~a3 in the
            # data frame b, it must be either {a3 ~ a1, a1 ~ a3, a3 ~ a2, a2 ~a3}
            # or {a3 ~ a1, a3 ~ a2}. As a result, unidirectional associations
            # does not produce differences in the sum of d_in_n in this algorithm.
            r$n_xy <- sum(b$n_xy)  # Will be greater than n_x and n_y.
            r$m_in <- round(r$d_in_n / r$n_xy, digits = 6)  # In fact, the equation is [sum(b$n_d) / 2] / [sum(b$n_xy) / 2] as every two rows are symmetric.
            r$c <- max(b$c)  # A connection to one, then a connection to all.
            r$s_d <- r$m_in * r$c  # s_a must be the same in b as alleles are idd; c = 1; r$m_in must always be one as all the edges are in perfect physical linkage (whose m_in = 1).
            r$s <- r$s_a + r$s_d  # So basically s_d = c here.
            r$pIBD_in <- max(b$pIBD_in)
        }  # Else, do nothing.
        merged <- rbind.data.frame(merged, r, stringsAsFactors = FALSE)
    }
    merged <- merged[order(merged$pair, decreasing = FALSE), ]

    return(merged)
}
