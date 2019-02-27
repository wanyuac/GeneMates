#' @title Make a ring plot to show presence/absence of genotypes and allelic co-occurrence
#'
#' @description This function requires the ggtree package to run. RGBA (RGB + alpha channel) colours are used by default.
#' It is imortant to note that random effects are always unobservable in LMMs, however, we can deduce the posterior
#' distribution of the coefficents of random effects. Therefore, the significance calls from our Bayesian chi-square test
#' does not refer to the direction of the correlation between a PC and the response variable.
#' As a result, some clades that are void of the response allele (Y) may be shaded in the
#' picture when corresponding principal components (PCs) are negatively contributed to the observed responses. This is
#' not an error because the analysis of structural random effects only concerns whether a projection
#' vector (along a PC) contributes to presence-absence of the response allele.
#'
#' The rings are filled outwards from the tree. For example, assuming the genotype
#' of interest is defined as c("a1", "a2", "a3"), then the inner-most ring will
#' represent the presence/absence of the allele a1, and the out-most ring will
#' represent the third allele a3.
#'
#' @param pam A presence/absence matrix (PAM) of genotypes (either alleles or genes)
#' @param genotypes A named or un-named list of character vectors. Each vector
#' specifies genotypes (e.g. alleles) whose co-occurrence will be plotted.
#' @param tree A phylo object (cf. the ape package) for a bifurcating tree. It can be a phylogenetic tree for samples
#' or a neighbour-joining tree based on projections of samples.
#' @param struc.eff (optional) The data frame "eff" in the output list of the function testForStruEff.
#' @param clade.cor (optional) A data frame for the correlation between clades and sample projections. The output data frame of the function
#' corCladeProj is an expected input.
#' @param clade.sizes (optional) A named integer vector storing the size of each clade.
#' This function transfers the sizes into the output when this parameter is set, making
#' the result more informative (A user may find it easier to find out which shaded clade
#' is most correlated with which projection vector along a principal component).
#' @param struc.pmax (optional) Maximum of p-values for structural effects to call significant.
#' @param struc.nmax (optional) Maximal number of significant structural effects to be plotted.
#' Default: 10. It cannot exceed 25.
#' @param y (optional) Name for a single allele or gene that is considered as the response variable. Branches in the tree will be
#' coloured for it when there are sample projections significantly correlated with it.
#' @param y.pat (optional) Pattern ID of the y allele or gene. It can be retrieved from the data frame "alle.pat" in the element "alleles" in the
#' output of findPhysLink. This argument must be provided when a user wants to colour branches by significant structural effects.
#' @param genotype.cluster (optional) A logical value determining whether to perform hierarchical clustering of genotypes (columns of the PAM).
#' Not applicable when PAM has less than three columns. Notice the order of inner rings may not follow the original one specified in the genotypes list.
#' @param genotype.dist (optional) A string specifying which distance metric is used for the clustering. Default: binary
#' @param cluster.method (optional) A string specifying which clustering method is used. Default: single
#' @param y.colour (optional) A single colour for the y genotype variable.
#' @param x.colours (optional) A single colour or a named (by allele names) colour vector for all explanatory alleles. It must not be white.
#' @param co.colours (optional) A single or multiple colours for co-occurrence data. They must not contain white.
#' @param null.colour (optional) A single baseline colour for absence of every allele. Default: grey90. Users may choose "white" as an alternative
#' when co-occurrence events are relatively common.
#' @param highlight.tips (optional) A character vector of names of tips to be highlighted with coloured circles.
#' @param highlight.tip.colour (optional) One (an unnamed character vector) or more colours (a vector of colours named by tip labels) for highlighted tips.
#' @param highlight.tip.shape (optional) An integer specifying the shape of highlighted tips, which follows the standard pch argument for R plots. Default: 16.
#' @param highlight.tip.size (optional) An integer specifying the size of highlighted tips. Default: 1.
#' @param highlight.tip.alpha (optional) A numeric specifying the alpha of the tip symbol. Default: 0.75.
#' @param clade.colours (optional) A vector of colours for clades (at most 10) that are most correlated with projections that significantly
#' contribute to the response variable y.
#' @param output (optional) Path and name for the output image.
#' @param is.pdf (optional) A logical value specifying if a PDF will be generated instead of a PNG file. Default: FALSE.
#' @param res (optional) Resolution of the output figure. Default: 72 ppi.
#' @param width (optional) Width of the output image.
#' @param height (optional) Height of the output image.
#' @param unit (optional) The unit of the width and height of the output image. Valid values: "mm" and "px" (default).
#' @param htmap.width (optional) Directly passed to the width parameter of the pheatmap function.
#' @param offset (optional) A parameter directly passed to the pheatmap function in ggtree for the offset argument.
#' @param branch.width (optional) Width of branches in the tree.
#' @param font.size (optional) Size of column names printed on the heat map.
#' @param print.colnames (optional) A logical parameter specifying whether to print column names on the heat map.
#' @param show.ledgend (optional) A logical parameter specifying whether to display the legend for components that are significantly correlated with
#' the response allele y. This is a legend for the heat map. Default: FALSE.
#'
#' @examples
#' # Example 1
#' ringPlotPAM(pam = assoc[["alleles"]][["A"]], genotypes = list(c1 = c("SulI_1616", "DfrA12_1089"),
#' c2 = c("SulI_1616", "DfrA12_1089", "AadA2_1605.1158")), tree = tr, x.colours = "grey50",
#' co.colours = c("blue", "red"))
#'
#' # Example 2
#' y <- "CmlA5_1538"
#' y.pat = assoc[["alleles"]][["alle.pat"]]$pattern[which(assoc[["alleles"]][["alle.pat"]]$allele == y)]
#' View(subset(assoc[["struc"]][["eff"]], y_pat == y.pat & p_adj <= 0.05))
#'
#' rp <- ringPlotPAM(pam = assoc[["alleles"]][["A"]],
#'                  genotypes = list(c1 = c("Arr2_274", "CmlA5_1538"),
#'                                   c2 = c("CTX-M-15_150", "VEB-1_1435")),
#'                  y = y, y.pat = y.pat, tree = assoc[["struc"]][["C"]][["tr"]],
#'                  struc.eff = assoc[["struc"]][["eff"]],
#'                  clade.cor = assoc[["struc"]][["cor"]],
#'                  clade.sizes = assoc[["struc"]][["clades"]][["size"]],
#'                  struc.pmax = 0.05, struc.nmax = 20,
#'                  genotype.dist = "binary", cluster.method = "single",
#'                  x.colours = "grey50", null.colour = "gray90",
#'                  co.colours = c("#d73027", "#fc8d59"), y.colour = "grey10",
#'                  output = "ringPlot_arr2cmlA5Cluster_2017080203.png",
#'                  branch.width = 0.25, htmap.width = 0.5, width = 190, height = 190, unit = "mm",
#'                  res = 300, print.colnames = FALSE)
#'
#' View(rp[["top"]])
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Dependency: ape, ggplot2, ggtree
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 23 June 2017, lastest edition: 28 February 2018

ringPlotPAM <- function(pam, genotypes, tree, y = NULL, y.pat = NULL,
                        struc.eff = NULL, clade.cor = NULL, clade.sizes = NULL,
                        struc.pmax = 0.05, struc.nmax = 10,
                        genotype.cluster = TRUE, genotype.dist = "binary",
                        cluster.method = "single", x.colours = "grey50",
                        co.colours = "red", null.colour = "grey90", y.colour = "grey10",
                        highlight.tips = NULL, highlight.tip.colour = "red",
                        highlight.tip.shape = 16, highlight.tip.size = 1,
                        highlight.tip.alpha = 0.75, clade.colours = rainbow(10),
                        output = "ringPlot.png", is.pdf = FALSE,
                        res = 72, width = 1600, height = 1600, unit = "px",
                        htmap.width = 0.5, offset = -0.001, branch.width = 0.25,
                        font.size = 2, print.colnames = TRUE, show.legend = FALSE) {

    if (is.character(genotypes)) {  # when there is a single group of alleles whose co-occurrence is to be visualised
        genotypes <- list(genotypes)
    }
    gs <- .getAllGenotypes(genotypes)  # get an overall set of genotypes (Each element is a character vector)
    pam <- pam[tree$tip.label, gs]  # extract presence/absence status of the targets and matches rows of the PAM to the tree
    x.colour.num <- length(x.colours)  # >= 1

    # determine whether to highlight tips
    if (is.character(highlight.tips)) {
        if (any(!(highlight.tips %in% tree$tip.label))) {
            # Do not highlight tips in this circumstance, otherwise, users may not be aware of their mistakes.
            print("Warning: no tips will be highlighted because some of them are not present in tip labels.")
            paint.tips <- FALSE  # No tips will be highlighted.
        } else {
            paint.tips <- TRUE
            to.highlight <- tree$tip.label %in% highlight.tips
            tree.annot <- data.frame(strain = tree$tip.label,
                                     shape = as.integer(sapply(to.highlight, function(x) ifelse(x, highlight.tip.shape, NA))),
                                     size = as.integer(sapply(to.highlight, function(x) ifelse(x, highlight.tip.size, NA))),
                                     stringsAsFactors = FALSE)  # annotations of the input tree
            if (is.null(names(highlight.tip.colour))) {  # names = NULL when highlight.tip.colour is not a named vector.
                tree.annot$colour <- as.character(sapply(to.highlight, function(x) ifelse(x, highlight.tip.colour, NA)))
            } else {  # multiple colours for highlighted tips
                names(to.highlight) <- tree$tip.label
                tree.annot$colour <- as.character(sapply(tree.annot$strain, function(x) ifelse(to.highlight[[x]], highlight.tip.colour[[x]], NA)))
            }
        }
    } else {
        paint.tips <- FALSE
    }

    # perform hierarchical clustering to columns of pam
    if (genotype.cluster & ncol(pam) >= 3) {
        hc <- hclust(d = dist(t(pam), method = genotype.dist), method = cluster.method)
        pam <- pam[, hc$order]
    }

    co <- .calcCooccurrenceMatrix(pam = pam, genotypes = genotypes)  # calculate a co-occurrence matrix

    # configuration of colours for co-occurrence data
    co.colours.n <- length(co.colours)  # Notice co.colours may contain the same colour of x.colours.
    co.n <- ncol(co)  # number of co-occurrence groups
    if (co.colours.n < co.n) {
        print("Warning: there are less colours in co.colours than the actual number of co-occurrence groups.")
        print("Additional columns of co-occurrence data will be filled with the first colour in x.colours.")
        co.colours <- c(co.colours, rep(x.colours[1], times = co.n - co.colours.n))
    }

    # encode 0's and 1's in the resulting matrix with colour codes so that ggtree can make a heat map for the matrix
    sample.n <- nrow(co)  # number of samples in the co-occurrence matrix
    y.provided <- (class(y) == "character") && (length(y) > 0)
    if (y.provided) {
        y <- y[1]  # In case a user provides multiple names for y; other y's will be ignored.
        y.provided <- y.provided && (y %in% colnames(pam))  # Otherwise, pam[, y] fails.
    }

    # assign colours to explanatory and response alleles
    if (y.provided) {  # assign a colour code for y when it is specified
        if (x.colour.num > 1) {
            if (is.null(names(x.colours))) {
                stop("Argument error: x.colours must be a named vector when its length is greater than one.")
            }
            pam.colours <- append(x.colours, y.colour)  # assuming x.colours is a named vector
            pam.colour.n <- length(pam.colours)
            names(pam.colours)[pam.colour.n] <- y  # Y allele's name
            for (colour.code in 1 : pam.colour.n) {
                a <- names(pam.colours[colour.code])  # allele name
                pam[, a] <- pam[, a] * colour.code
            }
        } else if (x.colour.num == 1) {
            pam.colours <- unique(x.colours, y.colour)  # A single colour is returned if x.colours = y.colour. as.character(x.colours): remove names from the vector.
            colour.code <- length(pam.colours)  # codes: 0: absence; 1: explanatory alleles; 2: response allele.
            pam[, y] <- pam[, y] * colour.code  # Do not need to set entry values for x as they equal one where the Y allele is present.
        } else {
            stop("Argument error: x.colours must have at least one colour code.")
        }
        colour.code.max <- colour.code

        # assign a colour code each column of the co-occurrence matrix
        for (i in 1 : co.n) {
            co.colour <- co.colours[i]
            if (co.colour %in% pam.colours) {  # when a new colour is to be assigned
                colour.code <- which(pam.colours == co.colour)[1]
            } else {
                colour.code.max <- colour.code.max + 1
                colour.code <- colour.code.max
            }
            co[, i] <- co[, i] * colour.code  # eg. c(1, 0, 1) * 3 = c(3, 0, 3)
        }
        # colour levels for the heat map
        #   0: null; 1 - n: x's in PAM; n + 1: y in PAM when y.colour is not in x.colours; others: co-occurrence groups
        htmap.colours <- c(null.colour, pam.colours, co.colours[which(!(co.colours %in% pam.colours))])
    } else {  # All alleles are coloured the same when the y allele is not specified.
        print("The response genotype y is not properly specified.")
        colour.code <- 1  # for x.colours; colour.code = 0 denotes the colour of absence
        x.colour <- x.colours[1]
        for (i in 1 : co.n) {
            if (co.colours[i] != x.colours) {  # then assign a different colour
                colour.code <- colour.code + 1
                co[, i] <- co[, i] * colour.code
            }  # skip the column if the corresponding colour equals x.colours
        }

        # colour levels for the heat map
        # 0: null, 1: PAM, others: co-occurrence groups
        htmap.colours <- c(null.colour, x.colour, co.colours[which(co.colours != x.colours)])
    }

    lv <- 0 : colour.code  # levels of values
    lv.labels <- as.character(lv)
    names(htmap.colours) <- lv.labels  # convert this colours into a named vector

    # merge matrices to make a single heat map
    # This is a critical step for colouring the heat map.
    htmap <- as.data.frame(cbind(pam, co))  # Row names are retained during this conversion.
    for (i in names(htmap)) {
        htmap[, i] <- factor(htmap[, i], levels = lv, labels = lv.labels)
    }

    # assign colours to clades following significant structural effects
    if (y.provided & !is.null(y.pat) & !is.null(struc.eff) & !is.null(clade.cor)) {
        print("Colouring clades in the tree by significant structural effects.")

        # sanity check of struc.nmax
        if (!is.integer(struc.nmax)) {
            struc.nmax <- as.integer(struc.nmax)
        }
        if (struc.nmax < 1 | struc.nmax > 25) {
            print("Warning: struc.max is reset to 10 as the current specification is invalid.")
            struc.nmax <- 10
        }

        # find out clades to be highlighted
        sig.eff <- subset(struc.eff, y_pat == y.pat & p_adj <= struc.pmax)  # significant structural effects as a part of random effects in an LMM
        proj.n <- nrow(sig.eff)  # number of significant axes
        if (proj.n > 0) {
            # sort projection vectors by their adjusted p-values in an ascending order
            # Although the data frame is sorted in the same way in the function lmm, it is safe
            # to sort it again when users provide their own data frame for structural effects.
            sig.eff <- sig.eff[order(sig.eff$p_adj, decreasing = FALSE), ]
            sig.proj <- paste0("c", sig.eff$c)  # append a "c" character to axis IDs. For instance, "1" becomes "c1".
            if (proj.n > struc.nmax) {  # only retain the top-20 significant effects
                print(paste0("Only keep the top-", struc.nmax, " most significant projection vectors."))
                proj.n <- struc.nmax
                sig.eff <- sig.eff[1 : proj.n, ]
                sig.proj <- sig.proj[1 : proj.n]  # a character vector
            }

            # check validity of a user's colour specification
            if (length(clade.colours) < proj.n) {
                print("Warning: there are less user-specified colours than the number of clades to be coloured.")
                print("Replacing clade colours with rainbow colours.")
                clade.colours <- rainbow(n = struc.nmax, s = 1, v = 1)  # RGBA (RGB + alpha) colour code
            }

            # find out the most correlated clade for each significant projection axis (in parallel with a PC)
            # Notice sig.proj is sorted by significance of projections on each axis.
            top.clades.perAxis <- lapply(sig.proj, .findTopClade, clade.cor)  # lapply returns an unnamed list of single-row data frames
            top.clades <- do.call(rbind.data.frame, top.clades.perAxis)
            names(top.clades)[1 : 2] <- c("top_axis", "top_clade")

            # append clade sizes to the result data frame
            if (!is.null(clade.sizes)) {  # must be a named integer
                top.clades$size <- as.integer(clade.sizes[top.clades$top_clade])
            }

            # finally, assign colour codes to clades
            top.clades$colour <- clade.colours[1 : proj.n]  # append the colour column to top.clades
        } else {
            print("No clade will be coloured as none of structural effects is significant for explaining the response variable.")
            top.clades <- NULL
        }
    } else {
        print("Skip colouring clades by significant structural effects as relevant information is not provided.")
        top.clades <- NULL
    }

    # initialise the tree in the ggtree system
    require(ape)
    require(ggtree)
    p <- ggtree(tree, layout = "circular", ladderize = TRUE, size = branch.width)

    # highlight tips where specified
    if (paint.tips) {
        p <- p %<+% tree.annot  # update the tree view with annotations
        # scale_shape_identity enables mapping a continuous variable to shape.
        p <- p + geom_tippoint(aes(shape = shape, size = size, color = colour),
                               alpha = highlight.tip.alpha) +
            scale_shape_identity() + scale_size_identity() + scale_color_identity()
    }

    # highlight top clades where present
    if (!is.null(top.clades)) {
        rownames(top.clades) <- NULL  # Row name information is useless here.

        # Sort clades in the data frame so that colour assignments are consistent across data sets
        top.clades <- top.clades[order(top.clades$size, abs(top.clades$cor), decreasing = TRUE), ]

        # Generate figure panels
        for (i in 1 : proj.n) {
            r <- top.clades[i, ]  # take a row from the data frame
            p <- p + geom_hilight(node = as.integer(gsub(pattern = "N", replacement = "",
                                                         x = r$top_clade, fixed = TRUE)),
                                  fill = r$colour, alpha = 0.6)
        }
    }

    # add a heat map to the tree
    p <- gheatmap(p = p, data = htmap, color = "white", colnames = print.colnames,
                  width = htmap.width, offset = offset, font.size = font.size) +
        coord_polar(theta = "y") +
        scale_fill_manual(values = htmap.colours, breaks = lv.labels)

    # To-do: the legend display is not ideal.
    if (show.legend) {
        p <- p + theme(legend.position = "bottom")  # a legend of the heat map
    } else {
        p <- p + theme(legend.position = "none")
    }

    # draw the plot
    if (is.pdf) {
        pdf(filename = output, paper = "a4", width = width, height = height)
    } else {
        png(filename = output, width = width, height = height, res = res, units = unit)
    }
    par(oma = rep(0.1, times = 4), mar = rep(0.1, times = 4))
    print(p)
    dev.off()

    return(list(pam = pam, co = co, htmap = htmap, y = y, top = top.clades))
}

# A subordinate function of ringPlotPAM. It pools alleles specified in vectors of the list
# genotypes into a single vector so that we can extract presence/absence information from
# the genotype PAM.
.getAllGenotypes <- function(genotypes) {
    targets <- NULL
    if (is.list(genotypes)) {
        n <- length(genotypes)  # number of co-occurrence groups
        for (i in 1 : n) {
            gs <- genotypes[[i]]  # get alleles of a co-occurrence group
            if (is.character(gs)) {
                if (length(gs) > 1) {
                    targets <- union(targets, gs)
                } else {
                    stop("Error: every elements of the \"genotype\" must not be a single genotype.")
                }
            } else {
                stop("Error: genotype names must be characters.")
            }
        }
    } else if (is.character(genotypes)) {  # a single character vector of genotype names
        if (length(genotypes) > 1) {
            targets <- genotypes
        } else {
            stop("Error: \"genotype\" must not contain a single genotype.")
        }
    } else {
        stop("Error: \"genotype\" must be either a list or a character vector.")
    }

    return(targets)  # returns a single vector of genotype names
}

# A subordinate function of ringPlotPAM. It returns a co-occurrence matrix (binary columns).
# The parameter "genotypes" must be a list.
.calcCooccurrenceMatrix <- function(pam, genotypes) {
    n <- length(genotypes)  # number of elements in the list
    groups <- names(genotypes)  # for calculating co-occurrence events between each group of genotypes
    if (is.null(groups)) {  # no group names available: use integers instead
        groups <- as.character(1 : n)
        names(genotypes) <- groups
    }
    co <- matrix(0, nrow = nrow(pam), ncol = n)  # initialise the binary co-occurrence matrix
    rownames(co) <- rownames(pam)
    colnames(co) <- groups
    for (g in groups) {
        m <- pam[, genotypes[[g]]]  # a sub-matrix only contains presence/absence status of current genotypes
        co[, g] <- as.integer(as.integer(rowSums(m)) == ncol(m))  # 1 (TRUE): genotypes are cooccurring in a sample; only works for binary matrices.
    }

    return(co)
}

# This is a subordinate function of ringPlotPAM. It returns a single-row data frame consisting of a clade name (the most correlated clade with
# projections on each axis) and its correlation coefficient to the projections.
.findTopClade <- function(proj, clade.cor) {
    clades <- subset(clade.cor, axis == proj)
    i <- which.max(abs(clades$cor))  # choose the first one if there are multiple hits (unlikely to happen)
    top.clade <- clades[i, c("axis", "clade", "cor")]

    return(top.clade)
}
