#' @title Summarising physical distances involving clusters of alleles
#'
#' @description This function is used for providing evidence for the presence of allele clusters as physical
#' clusters in bacterial strains where all of these alleles are co-occurring. It takes as input a data frame
#' of physical distance measurements and pulls out the distances corresponding to edges within each allele
#' cluster. It then summarises these distances for every strain where these distances are acquired as well as
#' for every cluster.
#'
#' @note All summaries and expected numbers are based on alleles that are actually co-occurring in a set of
#' strains. For instance, outputs only refer to five alleles if five out of six alleles in a cluster are consistently
#' co-occurring in corresponding strains. As such, explanations to the results should be related to alleles in the
#' 'alleles_max_co' element of the data frame 'cls.distr'.
#'
#' @param cls.distr A data frame produced by the function getClusterMemberCooccurrence for cluster distributions.
#' It shows strains where the most number of alleles are consistently co-occurring.
#' @param cls.col The name or index for the column of cluster IDs in cls.distr. Default: 1 (the first column).
#' @param allele.col The name or index for the column of allele IDs in cls.distr. Default: 2 (the second column).
#' @param cls A GraphSet-class object produced by the function extractSubgraphs, which lists edges per cluster.
#' @param ds A data frame of physical distances measured in all strains. It can be obtained from the data
#' frame "ds" in the output list of the function findPhysLink. The distances may be pre-filterred for a
#' maximal number of nodes or a maximal distance.
#' @param bidirectional A logical value specifying whether there are always two edges of opposing directions
#' connecting a pair of vertices. Default: TRUE.
#' @param sample.dists (optional) A square matrix for distances between samples. It can be acquired through the function projectSamples.
#' @param clade.pam (optional) a matrix for the presence/absence of samples in each clade of an input tree. It can be obtained using the function tree2Clades of phylix.
#' @param clade.sizes (optional) A named vector of integers for the number of samples in each clade. It can be obtained from the element "sizes"
#' in the outputs of the function tree2Clades. Optional.
#'
#' @examples
#' assoc <- findPhysLink(...)
#' clusters <- ...  # from a network package of your preference
#' com <- getClusterMemberCooccurrence(com = clusters, pam = assoc[["alleles"]][["A"]], cluster.colname = "community",
#' clade.pam = assoc[["struc"]][["clades"]][["pam"]], clade.sizes = assoc[["struc"]][["clades"]][["sizes"]],
#' sample.dists = assoc[["struc"]][["C"]][["d"]], n.cores = 4)
#' g_com <- extractSubgraphs(V = network[["V"]], E = network[["E"]], clusters = edges)  # cf. the documentation of this function
#' ds <- subset(assoc[["ds"]], node_number <= 3 & distance <= 2.8e6)
#' ds.com <- summarisePhysDistForClusters(cls.distr = com, cls.col = "community", cls.edges = g_com, ds = ds)
#'
#' Dependency: data.table
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export summarisePhysDistForClusters
#'
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 2 Oct 2017, the latest edition: 29 Oct 2018

summarisePhysDistForClusters <- function(cls.distr, cls.col = 1, allele.col = 2,
                                         cls, ds, bidirectional = TRUE,
                                         sample.dists = NULL, clade.pam = NULL,
                                         clade.sizes = NULL) {
    require(data.table)

    if (is.numeric(cls.col)) {  # is.integer does not work here.
        cls.col <- names(cls.distr)[cls.col]  # It is easier to have the column name rather than the index for programming.
    }
    if (is.numeric(allele.col)) {
        allele.col <- names(cls.distr)[allele.col]
    }

    # split the input data frame for empty and non-empty rows
    not_coocur <- is.na(cls.distr$clade)  # Some rows may not have constant co-occurring alleles recorded.
    na_rows <- cls.distr[not_coocur, ]
    na_n <- nrow(na_rows)
    cls.distr <- cls.distr[!not_coocur, ]
    cluster.ids <- cls.distr[, cls.col]

    # summarise physical distances for every cluster at the strain level
    phys <- mapply(.summarisePhysDistPerCluster, cl = cluster.ids,
                   alleles = cls.distr$alleles_max_co, strains = cls.distr$strains,
                   MoreArgs = list(ds, cls, bidirectional),
                   SIMPLIFY = FALSE)
    phys <- as.data.frame(rbindlist(phys), stringsAsFactors = FALSE)
    names(phys)[1] <- cls.col

    # append empty rows to phys for na_rows
    phys_add <- data.frame(cl = na_rows[, cls.col], strain = rep(NA, times = na_n), connectivity = rep(NA, times = na_n),
                           d_max = rep(NA, times = na_n), n_obs = rep(NA, times = na_n), n_exp = rep(NA, times = na_n),
                           missing_edges = rep(NA, times = na_n), stringsAsFactors = FALSE)  # n_exp is not calculated for simplicity although it can be determined.
    names(phys_add)[1] <- cls.col
    phys <- rbind.data.frame(phys, phys_add, stringsAsFactors = FALSE)
    phys <- phys[order(phys[, cls.col], decreasing = FALSE), ]

    # summarise the distances at the cluster level based on the data frame phys
    if (!is.null(sample.dists) & !is.null(clade.pam)) {
        if (is.null(clade.sizes)) {
            clade.sizes <- colSums(x = clade.pam)
        }
        cluster_summary <- lapply(unique(phys[, cls.col]), .clusterLevelSummary, phys, sample.dists, clade.pam, clade.sizes)
        cluster_summary <- as.data.frame(rbindlist(cluster_summary), stringsAsFactors = FALSE)
        names(cluster_summary)[1] <- cls.col

        # merge cluster_summary and cls.distr
        cluster_summary <- merge(x = cls.distr, y = cluster_summary, by = cls.col, all = TRUE, sort = TRUE)

        # rearrange columns to highlight key information
        new_order <- c(cls.col, "size", "co_max", "co_perc",
                       "strains_num", "strains_max_div", "clade_co_freq", "clade",
                       "clade_size", "strains100_num", "min_connectivity",
                       "strains100_max_div", "dmax_top", "dmax_median", "dmax_min",
                       "clade100_freq", "clade100", "clade100_size", allele.col,
                       "alleles_max_co", "strains", "strains100")
        colname_mismatch <- setdiff(new_order, names(cluster_summary))
        if (length(colname_mismatch) == 0) {  # ideal situation
            cluster_summary <- cluster_summary[, new_order]
        } else {
            stop(paste("Error: columns", colname_mismatch,
                       "are absent in the data frame new_order.", sep = " "))
        }

        out <- list(strain = phys, cluster = cluster_summary)  # strain-level and cluster-level summaries of physical distances
    } else {  # not run the summary process
        out <- phys
    }

    return(out)
}

# a subordinate function of summarisePhysDistForClusters
# Summarise distances for each cluster.
# cl: a single cluster ID, which is either a character or an integer value
.summarisePhysDistPerCluster <- function(cl, alleles, strains, ds, cls, bidirectional) {
    alleles_parsed <- strsplit(x = alleles, split = ",", fixed = TRUE)[[1]]  # including merged allele groups
    strains <- strsplit(x = strains, split = ",", fixed = TRUE)[[1]]

    # Sometimes in a linkage network we have a few idd. alleles that have been merged into a single group in the argument 'alleles'.
    # They do not change n_obs as they do not contribute to additional edges in networks where idd. alleles are merged into their
    # respective vertices. However, they may change the connectivity and maximal distance in their allele clusters. The following
    # section of codes lists identically distributed (idd.) alleles.
    if (grepl(pattern = "&", x = alleles, fixed = TRUE)) {
        alleles_idd <- .listIddAlleles(alleles = alleles_parsed, delim = "&")  # returns a named list of idd alleles
    } else {
        alleles_idd <- NULL
    }

    # expected edges of the current cluster
    edges_exp <- cls@E[[cl]]
    n <- nV(object = cls, id = cl)  # number of alleles (including merged groups) of the current subgraph

    # remove edges involving alleles that are not consistantly occurring in the current set of strains
    m <- length(alleles_parsed)  # number of singleton alleles + number of allele groups
    if (m < n) {
        edges_exp <- edges_exp[(edges_exp[, 1] %in% alleles_parsed) & (edges_exp[, 2] %in% alleles_parsed), ]
        n <- m
    }
    edges_num <- nrow(edges_exp)  # expected number of all edges

    # get the expected number of distance measurements (the same as the number of undirected edges connecting pairs of vertices) in every strain
    if (bidirectional) {
        n_exp <- edges_num / 2
    } else {  # undirected edges or a mixture of edge types
        n_exp <- length(unique(assignPairID(lmms = edges_exp, from = 1, paired.rows = TRUE)$pair))
    }

    # prepare the output data frame
    c <- data.frame(cl = vector(mode = ifelse(class(cl) == "character", "character", "integer"), length = 0),
                    strain = character(0), connectivity = numeric(0), d_max = integer(0),
                    n_obs = integer(0), n_exp = integer(0), missing_edges = character(0),
                    stringsAsFactors = FALSE)

    # extract distance measurements of edges in edges_exp for the current set of strains to reduce the number of comparisons
    ds <- subset(ds, sample %in% strains)

    # connectivity between these nodes (alleles) in each strain
    for (s in strains) {
        ds1 <- subset(ds, sample == s)
        selection <- rep(FALSE, times = nrow(ds1))  # selecting rows in ds
        edges_observed <- rep(FALSE, times = edges_num)  # selecting rows in edges_exp
        for (i in 1 : edges_num) {
            rw <- edges_exp[i, ]
            a1 <- rw[[1]]  # a1 and/or a2 may be an/both allele group(s), such as "TetRG_605.1265&TetG_632".
            a2 <- rw[[2]]

            # convert a merged item into a vector of individual alleles
            # Groups a1 and a2 cannot share any alleles, otherwise, they will be merged into a single vertex because of transitivity of distribution identity.
            if (grepl(pattern = "&", x = a1, fixed = TRUE)) {  # when this is a merged allele group
                a1 <- alleles_idd[[a1]]
            }
            if (grepl(pattern = "&", x = a2, fixed = TRUE)) {  # when this is a merged allele group
                a2 <- alleles_idd[[a2]]
            }

            # find rows in ds and flag rows in edges_exp
            observed <- FALSE
            for (q1 in a1) {
                for (q2 in a2) {
                    hits <- (ds1$query1 == q1 & ds1$query2 == q2) | (ds1$query1 == q2 & ds1$query2 == q1)  # 'hits' is a logical vector.
                    observed <- observed | any(hits)  # TRUE when this edge in edges_exp has any distance measured.
                    selection <- selection | hits
                }
            }
            edges_observed[i] <- observed
        }
        ds1 <- ds1[selection, ]  # May be empty when sum(selection) = 0 and equivalently, sum(edges_observed) = 0.

        # determine the number of observed edges
        if (bidirectional) {  # assuming edges are bidirectional between every pair of vertices
            n_obs <- sum(edges_observed) / 2  # Bandage returns one row of distance measurements per pair of alleles.
        } else {  # Edges are counted individually when they are admixed or undirected.
            n_obs <- sum(edges_observed)  # including scenarios where every subgraph is complete but its edges are undirected
        }

        # store summary statistics in data frames
        if (n_obs == n_exp) {  # The complete graph is present in the current strain.
            c <- rbind.data.frame(c, data.frame(cl = cl, strain = s, connectivity = round(n_obs / n_exp * 100, digits = 2),
                                                d_max = max(ds1$distance), n_obs = n_obs, n_exp = n_exp,
                                                missing_edges = "", stringsAsFactors = FALSE), stringsAsFactors = FALSE)
        } else if (n_obs > 0) {  # and nrow(ds1 < n_exp), then find out missing edges
            missing_pairs <- edges_exp[!edges_observed, ]
            c <- rbind.data.frame(c, data.frame(cl = cl, strain = s, connectivity = round(n_obs / n_exp * 100, digits = 2),
                                                d_max = max(ds1$distance), n_obs = n_obs, n_exp = n_exp,
                                                missing_edges = paste(paste(missing_pairs[, 1], missing_pairs[, 2], sep = "~"), collapse = ","),
                                                stringsAsFactors = FALSE), stringsAsFactors = FALSE)
        } else {  # all distances are not measurable in the current strain
            missing_pairs <- edges_exp
            c <- rbind.data.frame(c, data.frame(cl = cl, strain = s, connectivity = 0,
                                                d_max = NA, n_obs = n_obs, n_exp = n_exp,
                                                missing_edges = paste(paste(missing_pairs[, 1], missing_pairs[, 2], sep = "~"), collapse = ","),
                                                stringsAsFactors = FALSE), stringsAsFactors = FALSE)
        }
    }

    return(c)
}

.listIddAlleles <- function(alleles, delim = "&") {
    # This is a subordinate function for .summarisePhysDistPerCluster, which finds out identically distributed alleles
    # from a vector of allele/group names.

    alleles <- alleles[grepl(pattern = delim, x = alleles, fixed = TRUE)]  # keep items where the delimiter is found
    dict <- lapply(alleles, function(a) strsplit(x = a, split = delim, fixed = TRUE)[[1]])
    names(dict) <- alleles

    return(dict)
}

.clusterLevelSummary <- function(i, phys, sample.dists = NULL, clade.pam = NULL, clade.sizes = NULL) {
    # This is a subordinate function of summarisePhysDistForClusters. It works for a given cluster.
    # phys: a data frame summarising physical distances per cluster at the strain level
    # cls.distr: a data frame for cluster distributions in strains

    if (is.null(sample.dists) | is.null(clade.pam) | is.null(clade.sizes)) {
        stop("None of the arguments sample.dists, clade.pam and clade.sizes can be NULL.")
    }

    phys <- phys[phys[, 1] == i, ]  # only keep strain-level summaries of the current cluster
    if (!is.na(phys$strain[1])) {  # The cluster is physically present in at least one strain.
        min_connectivity <- min(phys$connectivity)
        if (nrow(phys) == 1) {  # when the cluster is only present in a single strain
            if (phys$connectivity[1] == 100) {
                strains100 <- phys$strain[1]
                strains100_num <- 1
                minc <- findMinIncCladeOfStrains(strains = strains100, clade.pam = clade.pam, clade.sizes = clade.sizes)
                clade100 <- minc[["clade"]]
                clade100_size <- minc[["size"]]
                clade100_freq <- minc[["freq"]]
                strains100_max_div <- 0  # maximal sample distance (divergence)
                dmax_top <- phys$d_max[1]
                dmax_median <- phys$d_max[1]
                dmax_min <- phys$d_max[1]
            } else {
                strains100 <- NA
                strains100_num <- 0
                clade100 <- NA
                clade100_size <- NA
                clade100_freq <- NA
                strains100_max_div <- NA
                dmax_top <- NA
                dmax_median <- NA
                dmax_min <- NA
            }
        } else {  # nrow(phys) > 1
            flag100 <- phys$connectivity == 100  # a logical vector for strains where connectivity equals 100%.
            if (any(flag100)) {  # There is at least one strain having physical distances of all edges measured.
                phys100 <- phys[flag100, ]
                strains100 <- phys100$strain
                strains100_num <- length(strains100)
                minc <- findMinIncCladeOfStrains(strains = strains100, clade.pam = clade.pam, clade.sizes = clade.sizes)
                clade100 <- minc[["clade"]]
                clade100_size <- minc[["size"]]
                clade100_freq <- minc[["freq"]]
                strains100_max_div <- max(sample.dists[strains100, strains100])  # maximal sample distance
                dmax_top <- max(phys100$d_max)
                dmax_median <- median(phys100$d_max)
                dmax_min <- min(phys100$d_max)
                strains100 <- paste(strains100, collapse = ",")
            } else {
                strains100 <- NA
                strains100_num <- 0
                clade100 <- NA
                clade100_size <- NA
                clade100_freq <- NA
                strains100_max_div <- NA
                dmax_top <- NA
                dmax_median <- NA
                dmax_min <- NA
            }
        }
    } else {  # an empty row
        min_connectivity <- NA
        strains100 <- NA
        strains100_num <- 0
        clade100 <- NA
        clade100_size <- NA
        clade100_freq <- NA
        strains100_max_div <- NA
        dmax_top <- NA
        dmax_median <- NA
        dmax_min <- NA
    }

    d <- data.frame(cl = i, strains100_num = strains100_num,
                    min_connectivity = min_connectivity, dmax_top = dmax_top, dmax_median = dmax_median,
                    dmax_min = dmax_min, strains100_max_div = strains100_max_div, clade100 = clade100,
                    clade100_size = clade100_size, clade100_freq = clade100_freq, strains100 = strains100,
                    stringsAsFactors = FALSE)

    return(d)
}
