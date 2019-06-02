# A module of GeneMates for analysing physical distances
# Copyright 2017-2018 Yu Wan (wanyuac@gmail.com)
# Licensed under the Apache License, Version 2.0
# Development history: 23 Janury - 27 March, 2017, 26 July 2017, ...
# Latest update: 2 June 2019

# Collapse any pair of directed edges into an undirected edge ###############
# One row of two allele names in an arbitrary order per pair of LMMs.
# edges: a two-column matrix or data frame
.collapseEdges <- function(edges) {
    d <- NULL  # to initialise a data frame of edges, one row of edges per pair of LMMs
    n <- nrow(edges)
    while (n > 1) {  # As we keep reducing the number of rows, the searching time gets shorter and shorter.
        q1 <- edges$query1[1]  # query 1; always start from the first row
        q2 <- edges$query2[1]  # query 2
        mate <- which(edges$query1 == q2 & edges$query2 == q1)  # find out the paired LMM on the inverse direction
        selected <- 1
        if (length(mate) > 0) {  # something is found - should alway be true as we fit LMMs on both directions
             selected <- append(selected, mate)
        }
        d <- rbind(d, edges[1, ])  # only push in a single row
        edges <- edges[-selected, ]  # but delete multiple rows of the same pair of alleles
        n <- nrow(edges)
    }
    if (n == 1) {  # n = 0 when the initial row number of the data frame "edges" is even.
        d <- rbind(d, edges)  # push in the last row
    }

    return(d)
}

# Remove edge redundancy. For example, edges (a, b) and (b, a) provide the same information about physical distances.
# Therefore, we merge edges (a, b) and (b, a) into a single edge (a, b).
# The direction of the collapsed edge does not matter because the distance is an undirected scalar.
.makeEdgeList <- function(edges) {
    # run within the function attachDistances
    if (is.matrix(edges)) {
        edges <- as.data.frame(edges, stringsAsFactors = FALSE)  # Function nrow(edges) returns NULL and 0 when "edges" is a matrix and it only has one row or zero row.
    }

    # only use the first two columns of the input
    if (ncol(edges) > 2) {
        edges <- edges[, 1 : 2]
    }
    colnames(edges) <- c("query1", "query2")  # The order of a and b alleles does not matter as the edge is undirected.

    # collapse bidirectional edges into a combination of alleles
    el <- .collapseEdges(edges)

    return(el)
}

# ap: a pair of alleles, which is a two-element vector of named characters.
# d: a merged data frame of distance measurements
.retriveDistanceMeasurements <- function(ap, d) {
    # implemented within the function attachDistances
    q1 <- ap[["query1"]]  # convert a named character into a pure character
    q2 <- ap[["query2"]]
    rows <- subset(d, (query1 == q1 & query2 == q2) | (query1 == q2 & query2 == q1))
    n <- nrow(rows)
    rows$query1 <- rep(q1, times = n)  # unify the order of query names
    rows$query2 <- rep(q2, times = n)

    return(rows)
}

# Attach distance measurements to LMM results under the alternative hypothesis H1
.attachDistances <- function(lmms = NULL, dists = NULL) {
    # check arguments
    if (is.null(lmms) | is.null(dists)) {
        stop("Error: both arguments (\"lmms\" and \"dists\") must be supplied.")
    }

    # initialise an edge list from the first two columns (assuming to be Y and X alleles) of lmms[[i]]
    for (i in names(lmms)) {  # i = "idd" or "dif"
        lmms.h1 <- lmms[[i]][["h1"]]
        rownames(lmms.h1) <- NULL
        el <- .makeEdgeList(lmms.h1[, 1 : 2])  # the edge list

        # attach distance measurementss onto the edge list
        lmms[[i]][["h1"]] <- do.call("rbind", apply(el, 1, .retriveDistanceMeasurements, dists))
    }

    # construct the return variable
    if (length(lmms) == 2) {
        r <- list(dif = lmms[["dif"]][["h1"]], idd = lmms[["idd"]][["h1"]])
    } else {  # when idd. alleles are absent
        r <- list(dif = lmms[["dif"]][["h1"]])
    }

    return(r)
}

# Find out the pair ID in lmms corresponding to a pair of alleles (query1, query2) in lmms.ds
.retrievePairID <- function(lmms.ds, lmms) {
    for (g in names(lmms.ds)) {  # Lmms and lmms.ds must have the same element names.
        ds <- lmms.ds[[g]]
        h1 <- lmms[[g]][["h1"]]
        ds$pair <- as.integer(apply(ds[, c("query1", "query2")], 1,
                                    function(r) h1$pair[which(h1$y == r[["query1"]] & h1$x == r[["query2"]])]))
        lmms.ds[[g]] <- ds[, c("query1", "query2", "pair", "sample", "distance",
                               "node_number", "source", "orientation", "distance_path")]
    }

    return(lmms.ds)
}

# Sort distance measurements and calculate their differences
.calcDeltaD <- function(d) {
    n <- length(d)
    if (n == 1) {
        dd <- NA
    } else {
        d <- sort(d, decreasing = FALSE)  # in an increasing order
        dd <- integer(n - 1)  # differences of distances: dd[i] = d[i + 1] - d[i]
        for (i in 1 : (n - 1)) {
            dd[i] <- d[i + 1] - d[i]
        }
    }

    return(dd)
}

# Summarising the distances ###############
.checkGroupIDs <- function(lmms, lmms.ds) {
    # A subordinate function of summariseDist.
    # Group IDs in lmms and lmms.ds must be the same.
    groups <- names(lmms)  # c("dif", "idd") or "dif"
    groups.ds <- names(lmms.ds)  # should be the same as groups
    if (sum(groups %in% groups.ds) < length(groups) | sum(groups.ds %in% groups) < length(groups.ds)) {
        stop("Argument error: element names in lmms and lmms.ds must be the same.")
    }

    return(groups)
}

.checkQuantileSetting <- function(qs) {
    # This is a subordinate function of summariseDist. It ensures 0%, Q1 (25th
    # percentile) and Q3 (75th percentile) are included in the vector.
    # qs: a numeric vector of quantiles

    # Check if there is any values smaller than zero.
    if (any(qs < 0) | any(qs > 1)) {
        stop("Argument error: quantile probabilities must not exceed the range [0, 1].")
    }

    # Check probabilities 0, 0.25, 0.75 and 1
    for (p in c(0, 0.25, 0.5, 0.75, 1)) {  # minimum, Q1, median, Q3 and maximum
        if (! p %in% qs) {
            qs <- append(x = qs, values = p)
        }
    }
    qs <- sort(qs, decreasing = FALSE)

    return(qs)
}

.checkLmms <- function(lmms, group) {
    # This is a subordinate function of summariseDist.
    # Check if every pair of alleles link to two rows in h1, which is a
    # prerequisite for this function.
    pairID.count <- unique(as.integer(table(lmms$pair)))
    if (any(pairID.count != 2)) {
        stop(paste0("[summariseDist] error: lmms[[", g, "]] must not contain unpaired rows."))
    }
}

.getDistQuantiles <- function(lmms, lmms.ds, qs, qnames, pam, tree, clades,
                              n.cores = -1) {
    # This is the second subordinate function of summariseDist. It calculate
    # quantiles of allelic physical distances for every pair of alleles.

    require(parallel)
    require(data.table)

    # Sanity check: match sample names of the clades matrix to those in the allelic PAM.
    samples.pam <- rownames(pam)
    if (any(rownames(clades[["pam"]]) != samples.pam)) {
        print("Warning: matching samples in clades to those in the allelic PAM.")
        clades[["pam"]] <- clades[["pam"]][samples.pam, ]
        clades[["sizes"]] <- clades[["sizes"]][samples.pam]  # Matching samples is critical for using the function .summariseMinIncCladeOfDists.
    }

    # Print run information
    pair.ids <- unique(lmms$pair)
    print(paste0("Quantiles to be calculated: ", paste(qs, collapse = ", "),
                 " for ", length(pair.ids), " pairs of alleles."))
    print(paste0("Column names for quantiles: ", paste(qnames, collapse = ", "), "."))

    # Compute summary statistics
    cl <- makeCluster(.setCoreNum(n.cores = n.cores, cores.avai = detectCores()))  # cores.avai >= 1; n.cores >= 1
    clusterExport(cl = cl,
                  varlist = list("lmms", "lmms.ds", "tree", "clades", "pam",
                                 "qs", "qnames"),
                  envir = environment())  # make variables accessible to different cores
    print(paste0(Sys.time(), ": starting computing summary statistics of distance measurements for each pair of alleles."))
    stats <- parLapply(cl, pair.ids, .calcSummaryStats,
                       lmms, lmms.ds, qs, qnames, tree, clades, pam)  # wait until all jobs finish
    stats <- as.data.frame(rbindlist(stats), stringsAsFactors = FALSE)  # concatenate rows into a single data frame; "b" denotes a slope.
    stopCluster(cl)
    stats <- stats[order(stats$pair, stats$y, decreasing = FALSE), ]
    print(paste0(Sys.time(), ": the computation is finished successfully."))

    return(stats)  # stats: summary statistics; ds: distances compared
}

.calcSummaryStats <- function(i, lmms, lmms.ds, qs, qnames, tree, clades, pam) {
    # This is a subordinate function of .getDistQuantiles.
    # Parameters:
    #   i: a pair ID
    #   ds: from lmms.ds
    #   tree: the same as the parameter "tree" specified in the function compareDist

    # Prepare the result data frame
    lmms <- subset(lmms, pair == i)  # It must have two rows (guaranteed by .checkLmms).
    qs.num <- length(qnames)
    qm <- matrix(data = NA, nrow = 2, ncol = qs.num)  # quantiles of the distances
    colnames(qm) <- qnames

    # Append columns for IQR (interquantile range) of the distances
    # m: measurability, which equls n(distances) / n_xy;
    # m_in: n(in-group distances) / n_xy
    # d_n, d_in/out_n: number of all, in-group/outlier distances
    # d_in_min/max: lower/upper whisker of in-group distances
    # pIBD_in: probability of the most common ancestor of strains having the
    # in-group distances to have the same or similar distance.
    # The following statement also defines default values of the columns appended.
    # Set c(0, 0) for d_in_n, d_out_n and m_in for convenience of data visualisation.
    lmms <- cbind.data.frame(lmms,
                             d_n = integer(2),
                             d_in_n = integer(2),
                             d_out_n = integer(2),
                             m = numeric(2),
                             m_in = numeric(2),
                             d_range = c(NA, NA),
                             d_IQR = c(NA, NA),
                             d_in_range = c(NA, NA),
                             d_in_min = c(NA, NA),
                             d_in_max = c(NA, NA),
                             pIBD_in = c(NA, NA),
                             d_in_mic = c(NA, NA),
                             d_in_mic_size = c(NA, NA),
                             d_in_mic_nxy = c(NA, NA),
                             d_in_mic_fxy = c(NA, NA),
                             as.data.frame(qm, stringsAsFactors = FALSE),
                             stringsAsFactors = FALSE)

    # Extract physical distances for the current pair
    n <- nrow(lmms.ds)  # number of distances measured, one in each sample
    if (n > 0) {
        lmms.ds <- subset(lmms.ds, pair == i)
        n <- nrow(lmms.ds)
    }  # Else: keep the default values.

    # Calculate summary statistics
    if (n > 0) {  # when some distance measurements are obtained
        # Get allele names
        a1 <- lmms$y[1]  # "allele 1"
        a2 <- lmms$x[1]  # "allele 2"
        lmms$d_n <- rep(n, times = 2)  # number of distances measured
        lmms$m <- round(lmms$d_n / lmms$n_xy, digits = 6)  # measurability

        # Calculate quantiles of the distances
        # Particularly, all quantiles are the same when there is only a single
        # distance. Therefore, this distance must belong to the in-group.
        Q <- quantile(x = lmms.ds$distance, probs = qs)  # a named numeric vector
        lmms[, qnames] <- rep(as.numeric(Q), each = 2)  # fill the selected block by columns
        Q1 <- Q[["25%"]]
        Q3 <- Q[["75%"]]
        IQR <- Q3 - Q1
        lmms$d_range <- rep(Q[["100%"]] - Q[["0%"]], times = 2)  # maximum - minimum
        lmms$d_IQR <- rep(IQR, times = 2)

        # Find out in-group and outlier distances
        lmms.ds <- .getDistanceGroups(lmms.ds, IQR, Q1, Q3)  # element names: "ingroup" and "outgroup"
        n.ingroup <- nrow(lmms.ds[["ingroup"]])
        n.outgroup <- nrow(lmms.ds[["outgroup"]])
        lmms[, c("d_in_n", "d_out_n")] <- rep(c(n.ingroup, n.outgroup), each = 2)
        lmms$m_in <- round(lmms$d_in_n / lmms$n_xy, digits = 6)
        lmms$d_in_min <- min(lmms.ds[["ingroup"]]$distance)
        lmms$d_in_max <- max(lmms.ds[["ingroup"]]$distance)
        lmms$d_in_range <- lmms$d_in_max - lmms$d_in_min

        # Structural analysis
        samples.in <- unique(lmms.ds[["ingroup"]]$sample)  # samples of in-group distances
        d.mic <- .summariseMinIncCladeOfDists(q1 = a1, q2 = a2, samples = samples.in,
                                              clades = clades, pam = pam)
        lmms$d_in_mic <- rep(d.mic$clade, times = 2)
        lmms$d_in_mic_size <- rep(d.mic$size, times = 2)
        lmms$d_in_mic_nxy <- rep(d.mic$n_xy, times = 2)
        lmms$d_in_mic_fxy <- rep(d.mic$f_xy, times = 2)

        # Determine if all distances are identical (consistent) by descendent
        # using ancestral state reconstruction
        IBD.p <- .estimateIBD(cid = d.mic$clade, samples = samples.in,
                              tr = tree)  # returns a probability
        lmms$pIBD_in <- rep(IBD.p, times = 2)
    }  # Else, keep the default data frame.

    return(lmms)
}

.getDistanceGroups <- function(ds, IQR, Q1, Q3) {
    # This is a subordinate function of .calcSummaryStats. It finds out in-group
    # and outlier distances.
    r <- 1.5 * IQR  # bound for the in-group
    ingroup.marks <- (ds$distance >= Q1 - r) & (ds$distance <= Q3 + r)
    d.groups <- list(ingroup = ds[ingroup.marks, ],
                     outgroup = ds[!ingroup.marks, ])

    return(d.groups)
}

# A subordinate function of .calcSummaryStats
# Parameters:
#   samples: samples where physical distances are collected from
#   clades: a list of elements "pam" and "sizes" for presence-absence of samples
#   in every clade.
#   pam: an allelic presence-absence matrix
# Expect the orders of samples are the same in the clade PAM (clades) and the
# allelic PAM (pam), which is guaranteed in the parental function compareDist.
.summariseMinIncCladeOfDists <- function(q1, q2, samples, clades, pam) {
    # get co-occurrence events of this pair of alleles (q1, q1) in the whole population
    xy <- as.logical(pam[, q1]) & as.logical(pam[, q2])  # mark samples (rows) where q1 and q2 are co-occurring

    s <- rownames(clades[["pam"]]) %in% samples  # convert the vector of samples into a logical vector of presence/absence of these samples
    n.d <- length(samples)  # number of samples having distances measured between this pair of alleles; In extreme cases, n.d = 1.
    if (n.d > 2) {  # size of any clade >= 2 for a bifurcate tree
        candidate.clades <- names(clades[["sizes"]])[which(as.integer(clades[["sizes"]]) >= n.d)]  # Each candiate clade must have at least n samples.
        clades[["pam"]] <- clades[["pam"]][, candidate.clades]  # This step reduces the number of comparisons later on.
    } else {
        candidate.clades <- names(clades[["sizes"]])  # Require no other action when n.d = 1 or 2.
    }

    # Identify the minimum clade that having all distance measurements (the samples parameter)
    # clades.d[i] equals TRUE when the i-th clade is an inclusive clade of samples where the distances are measured.
    clades.d <- apply(clades[["pam"]], 2, function(c) sum(as.logical(c) | s) == sum(c))  # The equation returns TRUE when "samples" form a subset of a clade c.
    inc.clades <- candidate.clades[as.logical(clades.d)]  # IDs of all inclusive clades of samples.

    # get the name of the minimal inclusive clade of input samples
    min.clade.id <- inc.clades[which.min(as.integer(clades[["sizes"]][inc.clades]))]

    # compute the frequency of co-occurrence events
    min.clade.size <- as.integer(clades[["sizes"]][min.clade.id])  # It is more understandable than use clades[["sizes"]][[min.clade]].
    clade.xy <- as.logical(clades[["pam"]][, min.clade.id]) & xy  # co-occurrence events within the min.clade
    n.clade.xy <- sum(clade.xy)  # number of co-occurrence events within this clade
    f.xy <- round(n.clade.xy / min.clade.size, digits = 6)  # intra-clade co-occurrence frequency; f.xy >= n.d / min.clade.size as the distance may not present for some co-occurrence events.

    return(data.frame(clade = min.clade.id, size = min.clade.size,
                      n_xy = n.clade.xy, f_xy = f.xy, stringsAsFactors = FALSE))
}

.estimateIBD <- function(cid, samples, tr) {
    # A subordinate function of .calcSummaryStats.
    # cid: clade ID
	# samples: where the explanatory allele is present
    require(ape)

    nid <- as.integer(substr(cid, start = 2, stop = nchar(cid)))  # e.g. N203 -> 203, which is a node ID in the phylo object
    tr.mic <- ape::extract.clade(phy = tr, node = nid)  # a rooted tree corresponding to the minimum inclusive clade
    if (tr.mic$Nnode > 1) {  # at least two internal nodes
        if (length(tr.mic$tip.label) > length(samples)) {
            # It is an prerequisite of the ace function to have variation in x to
            # estimate parameters of the transition model. Otherwise, when x only
            # contains zeros (not applicable in this function) or ones, ace reports
            # an error that "Error in matexpo(Q * EL[j]): BLAS/LAPACK routine
            # 'DGEEV ' gave error code -13".
            x <- integer(length = length(tr.mic$tip.label))  # an integer vector of 0's
            names(x) <- tr.mic$tip.label
            x[samples] <- 1  # samples showing the "presence" status of a binary trait

            # Ancestral state reconstruction for discrete status under an assumption of
            # equal transition rates.
            asr <- ape::ace(x = x, phy = tr.mic, model = "ER", type = "discrete")

            # Obtain the scaled likelihoods of each ancestral state, also known as the
            # marginal ancestral states and the empirical Bayesian posterior probabilities.
            # See http://www.phytools.org/eqg2015/asr.html for details.
            p <- round(as.numeric(asr$lik.anc[1, "1"]), digits = 8)
        } else {
            # Currently, we arbitrarily treat the ancestral state as the same
            # when all descendents acquired the same genomic structure. Herein we
			# cannot calculate the likelihood because we are not able to estimate the
            # transition rate when there is no variation in x. Notice the branch
			# length is not taken into account because the core-genome history
			# is independent to the acquisition of accessory gene clusters.
            p <- 1  # strongest potential to display vertical gene transfer in these samples
        }
    } else {
        # The internal node count Nnode = 1 happens when there are only two tips
		# in a sub-tree, indicating that both samples are sister descendents of
		# the same most-recent ancestor. In this case, the potential of vertical
		# transfer of the same genomic structure must reach its maximum of one.
		# Since the ace function reports an error of "'x' must be an array of at
        # least two dimensions" in this situation, I create this logical branch
        # to deal with this scenario.
        p <- 1
    }

    return(p)
}

.extractLMMoutputs <- function(lmms) {
    # This a subordinate function of findPhysLink and it works when evalPL is assessing
    # evidence of physical linkage without considering allelic physical distances.
    elements <- names(lmms)  # c("dif", "idd") or "dif"
    out <- vector(mode = "list", length = length(elements))
    names(out) <- elements
    for (i in elements) {
        out[[i]] <- lmms[[i]][["h1"]]
    }

    return(out)
}
