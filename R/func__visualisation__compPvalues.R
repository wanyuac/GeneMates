#' @title Draw a scatter plot and histograms to compare p-values from linear mixed models
#' (LMMs) and penalised logistic models (PLMs).
#'
#' @description This function only works at the pattern level. P-values are grouped
#' by log10 lambda0 of X and Y.
#' Dependencies: ggplot2, grid and gridExtra.
#'
#' @param p.lmm A data frame from the function lmm for association status between
#' patterns. Assuming assoc = lmm(...), then p.lmm = assoc$lmms.pat$dif$h1.
#' @param p.plm A data frame from the function plr for association status between
#' patterns. Assuming assoc = plr(...), then p.plm = assoc$pat. Pattern orders for
#' association tests must match to those in p.lmm. Otherwise, both data frames
#' cannot merge correctly.
#' @param lmm.h0 A data frame from the function lmm for null models. lmm.h0 =
#' assoc$lmms.pat$dif$h0.
#' @param p.min Minimum of raw p-values that can be represented in computers
#' precisely. Default: 2.2e-16.
#' @param p.max Maximum of Bonferroni-corrected p-values for concluding
#' significance. Default: 0.05.
#' @param L.weak Maximum log10 lambda for defining weak structural random effects.
#' Default: 1.
#' @param bks Breaks for X and Y axes in the plot. Default: 0, 2, 4, 6, ..., 16.
#' @param show.p.max A logical argument determine whether to draw a grey dashed
#' line on each axis to indicate the cutoff for significance. Default: TRUE.
#' @param cols A character vector with names "0", "1", "2" and "3" for definition
#' of colours for lambda groups.
#' @param plot.title.size Size of the title of each panel.
#' @param axis.text.size Size of axis labels.
#' @param axis.title.size Size of axis titles.
#' @param axis.title.size.panelB Title sizes for panel b, which should be much
#' smaller than axis.title.size and controls the alignment of panels.
#' @param img Name and path for the output PNG file. Default: comp_pvals.png.
#' @param img.w Image width.
#' @param img.h Image height.
#' @param img.r Image resolution.
#' @param img.u Unit for the image width and height.
#' @param img.oma The argument oma for the function par during plotting.
#' @param img.mar The argument mar for the function par.
#' @param img.mgp The argument mgp for the function par.
#' @param p.prev Output of this function from a previous run. It overrides p.lmm
#' and p.plm.
#'
#' @return A data frame underlying the output plot.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 13 Sep 2018, the lastest edition: 13 Sep 2018

compPvalues <- function(p.lmm, p.plm, lmm.h0, p.min = 2.2e-16, p.max = 0.05,
                        L.weak = 1, bks = seq(0, 16, by = 2), show.p.max = TRUE,
                        cols = c("0" = "red", "1" = "orange", "2" = "blue",
                                 "3" = "grey50"), plot.title.size = 10,
                        axis.text.size = 8, axis.title.size = 10,
                        axis.title.size.panelB = 3, img = "comp_pvals.png",
                        img.w = 150, img.h = 150, img.r = 300, img.u = "mm",
                        img.oma = c(0.1, 0.1, 0.1, 0.1), img.mar = c(3, 2.6, 1, 0.5),
                        img.mgp = c(1.8, 0.6, -0.2),
                        p.prev = NULL) {
    require(ggplot2)
    require(grid)
    require(gridExtra)

    # Merge data frames
    if (is.null(p.prev)) {
        p.lmm <- p.lmm[, c("y_pat", "x_pat", "p_wald")]
        p.plm <- p.plm[, c("y_pat", "x_pat", "p_chisq")]
        p <- merge(x = p.lmm, y = p.plm, by = c("y_pat", "x_pat"), all = TRUE, sort = FALSE)
        if (any(is.na(p$p_wald) | is.na(p$p_chisq))) {
            stop("Error: there are unmatched rows in p.lmm and p.plm.")
        }

        print(paste(nrow(p), "tests for each kind of models.", sep = " "))

        # Calculate negative log10 of p-values
        p$p_wald_comp <- p$p_wald >= p.min
        p$p_chisq_comp <- p$p_chisq >= p.min
        p <- subset(p, p_wald_comp & p_chisq_comp)
        p$p_wald_nevLog10 <- -log10(p$p_wald)
        p$p_chisq_nevLog10 <- -log10(p$p_chisq)

        # Retrieve lambda0 and perform Log10 transformation to lambdas
        p$L0_y <- lmm.h0$lambda0_remle[match(p$y_pat, lmm.h0$y_pat)]
        p$L0_x <- lmm.h0$lambda0_remle[match(p$x_pat, lmm.h0$y_pat)]

        # Group allele pairs based on their lambda0 values
        p$L0_y_group <- as.integer(p$L0_y <= L.weak)  # zero: strong structural random effect
        p$L0_x_group <- as.integer(p$L0_x <= L.weak) * 2
        p$L0_group <- p$L0_y_group + p$L0_x_group
        p$L0_group <- as.factor(p$L0_group)  # levels: 0, 1, 2, 3

        # Get the final data frame
        p <- p[, c("y_pat", "x_pat", "p_wald", "p_chisq", "p_wald_nevLog10", "p_chisq_nevLog10",
                   "L0_y", "L0_x", "L0_group")]
    } else {
        p <- p.prev
    }
    n <- nrow(p)
    print(paste(n, "pairs of p-values are comparable."))


    # Calculate an upper bound of p-values for significance based on remaining
    # pattern pairs.
    bks.max <- max(bks)
    p.max <- -log10(p.max / n)
    print(paste("Minimum -log10 p-value for significance in comparable tests:",
                round(p.max, digits = 4), sep = " "))

    # Prepare three ggplot objects
    a <- ggplot(data = p, mapping = aes(x = L0_group, y = p_wald_nevLog10,
                                        fill = L0_group))
    if (show.p.max) {
        a <- a + geom_hline(yintercept = p.max, linetype = "dashed", size = 0.5,
                            colour = "grey25")
    }
    a <- a + geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black") +
        labs(title = "a", x = "Group", y = expression(-log[10](p[2]))) +
        scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                           labels = as.character(bks), expand = c(0, 0)) +
        scale_fill_manual(breaks = as.character(0 : 3), values = cols) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold", size = plot.title.size),
              panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_line(colour = "grey95"),
              axis.text = element_text(size = axis.text.size, colour = "black"),
              axis.title = element_text(size = axis.title.size, colour = "black"),
              legend.position = "none")

    b <- ggplot(data = p,
                mapping = aes(x = p_chisq_nevLog10, y = p_wald_nevLog10,
                              colour = L0_group)) +
        geom_point(shape = 16, size = 1)
    if (show.p.max) {
        b <- b + geom_hline(yintercept = p.max, colour = "grey25", size = 0.5,
                            linetype = "dashed") +
            geom_vline(xintercept = p.max, colour = "grey25", size = 0.5,
                       linetype = "dashed")
    }
    b <- b + geom_abline(intercept = 0, slope = 1, size = 0.5, colour = "black") +
        labs(title = "b", x = expression(-log[10](p[1])),
             y = expression(-log[10](p[2]))) +
        scale_colour_manual(breaks = c("0", "1", "2", "3"), values = cols) +
        scale_x_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                           labels = as.character(bks), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                           labels = as.character(bks), expand = c(0, 0)) +
        theme_bw() +
        theme(plot.title = element_text(face = "bold", size = plot.title.size),
              panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_line(colour = "grey95"),
              axis.text = element_text(size = axis.text.size, colour = "black"),
              axis.title = element_text(size = axis.title.size.panelB,
                                        face = "bold", colour = "white"),
              legend.position = "none")

    c <- ggplot(data = p, mapping = aes(x = L0_group, y = p_chisq_nevLog10,
                                        fill = L0_group))
    if (show.p.max) {
        c <- c + geom_hline(yintercept = p.max, linetype = "dashed", size = 0.5,
                            colour = "grey25")
    }
    c <- c + geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black") +
        labs(title = "c", x = "Group", y = expression(-log[10](p[1]))) +
        scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                           labels = as.character(bks), expand = c(0, 0)) +
        scale_fill_manual(breaks = as.character(0 : 3), values = cols) +
        coord_flip() + theme_bw() +
        theme(plot.title = element_text(face = "bold", size = plot.title.size),
              panel.grid.major = element_line(colour = "grey90"),
              panel.grid.minor = element_line(colour = "grey95"),
              axis.text = element_text(size = axis.text.size, colour = "black"),
              axis.title = element_text(size = axis.title.size, colour = "black"),
              legend.position = "none")

    # Draw a figure
    png(filename = img, width = img.w, height = img.h, units = img.u, res = img.r)
    par(oma = img.oma, mar = img.mar, mgp = img.mgp)
    grid.arrange(a, b, textGrob(""), c, ncol = 2, widths = c(2, 5),
                 heights = c(5, 2))
    dev.off()

    return(p)
}
