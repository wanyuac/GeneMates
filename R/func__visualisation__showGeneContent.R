#' @title Draw a bubble plot and two bar plots to summarise gene and allele frequencies.
#'
#' @description Use the function countAlleles and calcGeneFreq to know allele and
#' gene content before using this function, in order to set appropriate graphic
#' arguments. The diameter of each bubble is porportional to the allele number.
#'
#' Dependency: ggplot2.
#'
#' @param af A table of allele frequencies, produced by the function countAlleles.
#' @param gf A table of gene frequencies, produced by the function calcGeneFreq.
#' @param d.min Diameter for genes of the minimum allele count.
#' @param d.max Diameter for genes of the largest allele count.
#' @param fill.colour A named vector of colours filling bubbles.
#' @param border.colour A single colour for borders of circles.
#' @param core.genes Names of core genes
#' @param bar.acc Whether to draw a bar plot of gene frequencies for accessory genes
#' (TRUE) or all genes (FALSE, default).
#' @param panel.names A character vector of three elements for panel names.
#' Default: a, b, c.
#' @param panelA.xmax Maximum on the X axis of panel a for gene frequencies.
#' Default: 120.
#' @param panelB.xinterv Distance between major ticks on the X axis of the first
#' bar plot (panel b). Default: 5.
#' @param panelB.yinterv Distance between major ticks on the Y axis of the first
#' bar plot (panel b). Default: 20.
#' @param panelC.xinterv Distance between major ticks on the X axis of the second
#' bar plot (panel b). Default: 50.
#' @param prev.out Output from a previous run of this function.
#' @param panelB.col Colour for a bar plot in the panel b. Default: grey50.
#' @param panelC.lwd An integer specifying the width of bars for allele frequencies
#' in the panel c. Default: 1.
#' @param f Name for the output PNG file. Default: gene_content.png.
#' @param w Image width.
#' @param h Image height.
#' @param r Image resolution.
#' @param u Unit for the image width and height.
#' @param img.oma The argument oma for the function par during plotting.
#' @param img.mar The argument mar for the function par.
#' @param img.mgp The argument mgp for the function par.
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 4 Sep 2018, the lastest edition: 7 Oct 2018

showGeneContent <- function(af, gf, pam.g, d.min = 2, d.max = 20,
                            fill.colour, border.colour = "grey90",
                            core.genes = NULL, bar.acc = FALSE,
                            panel.names = c("a", "b", "c"), panelA.xmax = 120,
                            panelB.xinterv = 5, panelB.yinterv = 20,
                            panelB.col = "grey50", panelC.lwd = 1, panelC.xinterv = 20,
                            prev.out = NULL, f = "gene_content.png",
                            w = 170, h = 200, r = 300, u = "mm",
                            img.oma = c(0.1, 0.1, 0.1, 0.1),
                            img.mar = c(3, 2.6, 1, 0.5),
                            img.mgp = c(1.8, 0.6, -0.2)) {
    if (is.null(prev.out)) {
        # Panel a: Gene frequencies and allele count per gene ===============
        # mapping allele counts to diameters via a linear model
        gf$diam <- .calcDiameters(gf$n_a, d.min, d.max)  # lib__visualisation.R
        gf <- gf[order(gf$class, decreasing = FALSE), ]
        gf$fill <- as.character(fill.colour[gf$class])
        gf$border <- border.colour

        # Try to draw hollow circules for core genes
        with_core <- !is.null(core.genes)
        if (with_core) {
            n_c <- length(core.genes)
            is_core <- gf$gene %in% core.genes
            gf$border[is_core] <- gf$fill[is_core]
            gf$fill[is_core] <- NA  # make circles unfilled
        }

        # Panel b: gene count per strain ===============
        n_all <- apply(pam.g, 1, function(r) sum(r != "-"))  # number of all genes per sample
        if (with_core) {
            # Core genes
            if (n_c > 1) {
                n_core <- apply(pam.g[, core.genes], 1, function(r) sum(r != "-"))  # number of intrinsic genes per sample
            } else {  # a single core gene
                n_core <- sapply(pam.g[, core.genes], function(r) as.integer(r != "-"))
            }

            # Accessory genes
            if (ncol(pam.g) - n_c > 1) {  # multiple accessory genes
                n_acc <- apply(pam.g[, !(colnames(pam.g) %in% core.genes)], 1, function(r) sum(r != "-"))
            } else {
                n_acc <- sapply(pam.g[, !(colnames(pam.g) %in% core.genes)], function(r) as.integer(r != "-"))
            }

            gene_num <- data.frame(strain = names(n_all), n = n_all, n_core = n_core,
                                   n_acc = n_acc, stringsAsFactors = FALSE)
        } else {
            gene_num <- data.frame(strain = names(n_all), n = n_all, stringsAsFactors = FALSE)
        }
        gene_num <- gene_num[order(gene_num$n, decreasing = TRUE), ]

        # Panel c: allele frequency ===============
        af$colour <- as.character(fill.colour[af$class])
        af <- af[order(af$class, af$freq, decreasing = TRUE), ]
    } else {
        gene_num <- prev.out[["gene_num"]]
        af <- prev.out[["af"]]
        gf <- prev.out[["gf"]]
        strain_num <- prev.out[["strain_num"]]
        with_core <- prev.out[["with_core"]]
        bar.acc <- prev.out[["bar_acc"]]
    }

    # Gene count per strain (continue), and also prepare allele counts for panel c
    if (with_core && bar.acc) {
        print("Draw bar plots of gene counts and allele frequencies of accessory genes.")
        strain_num <- table(gene_num$n_acc)  # to make a bar plot of the number of accessory ARGs per strain
        af <- subset(af, !gene %in% core.genes)
    } else {
        print("Draw bar plots of gene counts and allele frequencies of all genes.")
        strain_num <- table(gene_num$n)  # number of strains per gene count
    }
    xs <- as.integer(names(strain_num))  # gene count per strain
    ys <- as.integer(strain_num)
    y_max <- ceiling(max(ys) / panelB.yinterv) * panelB.yinterv
    x_max <- ceiling(max(xs) / panelB.xinterv) * panelB.xinterv
    n_g <- nrow(gf)  # number of genes
    n_a <- nrow(af)
    y_max_c <- ceiling(max(af$freq) / 10) * 10  # for panel c
    x_max_c <- ceiling(n_a / panelC.xinterv) * panelC.xinterv

    # Make the figure covering alleles ===============
    png(filename = f, width = w, height = h, res = r, units = u)
    layout(mat = matrix(c(1, 1, 2, 3), byrow = FALSE, ncol = 2))
    par(oma = img.oma, mar = img.mar, mgp = img.mgp)

    # Panel a: a bubble chart showing gene frequency, AMR classes and allele diversity
    plot(1, axes = FALSE, xlim = c(0, panelA.xmax), ylim = c(1, n_g), type = "n",
         ylab = "Gene", xlab = "Gene frequency (%)", cex.lab = 1)  # xlim > 100 to leave some space for circle
    symbols(x = gf$freq, y = 1 : n_g, circles = gf$diam, fg = gf$border,
            bg = gf$fill, add = TRUE, inches = 0.4)
    text(x = gf$freq, y = 1 : n_g, labels = gf$n_a, offset = 0, cex = 0.8, col = "black")
    axis(side = 1, at = seq(0, 100, by = 20), las = 1, cex.axis = 0.8)
    rug(side = 1, x = seq(10, 90, by = 20), ticksize = -0.025)
    title(main = panel.names[1], adj = 0, cex.main = 1)

    # Panel b at the bottom: a bar plot of gene count per strain
    # https://stackoverflow.com/questions/38002360/increasing-the-width-of-type-h-r-plot
    #plot(x = xs, y = ys, col = "black", type = "h", lwd = panelB.lwd, lend = 1,
    #     xlab = "Gene count", ylab = "Strain count", xlim = c(0, x_max),
    #     ylim = c(0, y_max), axes = FALSE, cex.lab = 1)
    # Determine bar heights at continuous integers
    ys_b <- rep(0, times = x_max + 1)  # default heights in the panel b
    for (i in 1 : length(xs)) {
        x <- xs[i]
        ys_b[x + 1] <- ys[i]  # xs and ys are matched.
    }

    # Produce labels on the X axis
    x_labs <- rep(NA, times = length(ys_b))
    prev_val <- 0  # the previous value written into x_labs
    x_labs[1] <- as.character(prev_val)  # the first label
    lab_counter <- 0
    for (i in 2 : length(x_labs)) {
        lab_counter <- lab_counter + 1
        # Since x_max is determined by panelB.xinterv, the last element written
        # into the x_labs must equal x_max.
        if (lab_counter == panelB.xinterv) {
            prev_val <- prev_val + panelB.xinterv
            x_labs[i] <- as.character(prev_val)
            lab_counter <- 0
        }
    }

    # Make the bar plot for panel b
    barplot(height = ys_b, names.arg = x_labs, col = panelB.col, axes = FALSE,
            xlab = "Gene count", ylab = "Strain count", cex.lab = 1,
            ylim = c(0, y_max))  # bars is a matrix of a single column, which stores the midpoint of each bar.
    #axis(side = 1, at = bars[, 1], labels = x_labs, cex.axis = 0.8, vadj = -0.5)
    axis(side = 2, at = seq(0, y_max, by = panelB.yinterv), cex.axis = 0.8)
    minor_y_interv_width <- panelB.yinterv / 2
    rug(side = 2, x = seq(minor_y_interv_width, y_max - minor_y_interv_width,
                          by = panelB.yinterv), ticksize = -0.025)
    title(main = panel.names[2], adj = 0, cex.main = 1)

    # Panel c: a bar plot showing allele frequency coloured by AMR classes
    plot(x = 1 : n_a, y = af$freq, type = "h", lwd = panelC.lwd, lend = 1,
         col = af$colour, xlab = "Allele index", ylab = "Allele frequency (%)",
         cex.lab = 1, xlim = c(0, x_max_c), ylim = c(0, y_max_c), axes = FALSE)
    axis(side = 1, at = seq(0, x_max_c, by = panelC.xinterv), las = 1, cex.axis = 0.8)
    axis(side = 2, at = seq(0, y_max_c, by = 10), las = 2, cex.axis = 0.8)
    rug(side = 2, x = seq(5, y_max_c - 5, by = 10), ticksize = -0.025)
    title(main = panel.names[3], adj = 0, cex.main = 1)

    dev.off()

    out <- list(gene_num = gene_num, af = af, gf = gf, strain_num = strain_num,
                with_core = with_core, bar_acc = bar.acc)

    return(out)
}
