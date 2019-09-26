#' @title Draw a scatter plot and histograms to compare p-values from linear mixed models
#' (LMMs) and penalised logistic models (PLMs).
#'
#' @description This function only works at the pattern level. There are two ways to group
#' p-values: (1) grouping by p-values when L.weak <= 0 (default); and (2) grouping by log10
#' lambda0 of X and Y when parameter L.weak > 0. Each way results in four groups, which by
#' default are coloured coded as red (group 0), purple (group 1), blue (group 2) and grey
#' (group 3). Note that this function assumes that the number of hypothesis tests equals the
#' number of rows in the input data frames p.lmm and p.plm, which must have the same number
#' of rows as well.
#'
#' Dependencies: ggplot2, grid, gridExtra, and reshape2.
#'
#' @param p.lmm A data frame from the function lmm for association status between
#' patterns. Assuming assoc = lmm(...), then p.lmm = assoc$lmms.pat$dif$h1.
#' @param p.plm A data frame from the function plr for association status between
#' patterns. Assuming assoc = plr(...), then p.plm = assoc$pat. Pattern orders for
#' association tests must match to those in p.lmm. Otherwise, both data frames
#' cannot merge correctly.
#' @param lmm.h0 A data frame from the function lmm for null models. lmm.h0 =
#' assoc$lmms.pat$dif$h0.
#' @param p.min Minimum of raw p-values that can be represented in users' computers
#' precisely. Default: 2.2e-16. Any p-value that is smaller than this cut-off will
#' be substituted by it in this function.
#' @param p.adj.max Maximum of Bonferroni-corrected p-values for concluding
#' significance. Default: 0.05.
#' @param L.weak Maximum of estimate lambda0 for defining weak structural random
#' effects (default: 1, whose log10 equals zero). This parameter is used for grouping
#' data points when L.weak > 0.
#' @param bks Breaks for X and Y axes in the plot. They correspond to -log10 transformed
#' raw p-values. Default value: c(0, 2, 4, 6, ..., 16).
#' @param show.p.adj.max A logical argument determine whether to draw a grey dashed
#' line on each axis to indicate the cutoff for significance. Default: TRUE.
#' @param cols A character vector with names "0", "1", "2" and "3" for definition
#' of colours for all four groups. When data points are grouped by p-values, 0: not
#' significant in either LMMs or PLMs; 1: only significant in PLMs; 2: only significant
#' in LMMs; and 3: significant in both LMMs and PLMs. When data points are grouped by
#' lambda0, 0: weak structural random effects in both X and Y; 1: strong structural
#' random effects only in X; 2: strong structural random effects only in Y; 3: strong
#' structural random effects in both X and Y.
#' @param panel.num An integer argument specifying how many panels are to be drawn in the
#' figure. By default, three panels (panel.num = 3: two box plots and one scatter plot)
#' are created. Otherwise, two panels (one scatter plot and one box plot with groups names
#' LMM and PLM) are drawn.
#' @param panel.titles A vector of two (panel.num = 2) or three (panel.num = 3) characters
#' for titles of figure panels. Default: c("a", "b", "c").
#' @param boxplot.show.group Whether to draw box plots showing every group of data points.
#' Default: TRUE. Two box plots each has a single box will be drawn if this parameter
#' is FALSE.
#' @param boxplot.default.fill Fill colour when boxplot.group = FALSE. Default: grey80.
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
#' @param p.prev Output of this function from a previous run. Other inputs will be
#' ignored if p.prev is specified.
#'
#' @return A data frame underlying the output plot.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export
#
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First version: 13 Sep 2018, the lastest edition: 25 Sep 2019

comparePvalues <- function(p.lmm, p.plm, lmm.h0 = NULL, p.min = 2.2e-16, p.adj.max = 0.05,
                           L.weak = 0, bks = seq(0, 16, by = 2), show.p.adj.max = TRUE,
                           cols = c("0" = "grey50", "1" = "blue", "2" = "red", "3" = "purple"),
                           panel.num = 3, panel.titles = c("a", "b", "c"),
                           boxplot.show.group = TRUE, boxplot.default.fill = "grey80",
                           plot.title.size = 10, axis.title.size = 10,
                           axis.title.size.panelB = 3, axis.text.size = 8,
                           img = "pvalue_comparison.png",
                           img.w = 150, img.h = 150, img.r = 300, img.u = "mm",
                           img.oma = c(0.1, 0.1, 0.1, 0.1), img.mar = c(3, 2.6, 1, 0.5),
                           img.mgp = c(1.8, 0.6, -0.2),
                           p.prev = NULL) {
    require(ggplot2)
    require(grid)
    require(gridExtra)

    if (length(panel.titles) < panel.num) {
        print("Warning: panel titles are fewer than the panel count. Replaced the titles with default values.")
        panel.titles <- c("a", "b", "c")
    }

    # Merge data frames ===============
    if (is.null(p.prev)) {
        p.lmm <- p.lmm[, c("y_pat", "x_pat", "p_wald")]  # raw p-values from LMMs
        p.plm <- p.plm[, c("y_pat", "x_pat", "p_chisq")]  # raw p-values from penalised logistic regression
        p <- merge(x = p.lmm, y = p.plm, by = c("y_pat", "x_pat"), all = TRUE, sort = FALSE)  # Assumption: both p.lmm and p.plm consist of the same pattern pairs.
        if (any(is.na(p$p_wald) | is.na(p$p_chisq))) {
            stop("Error: there are unmatched rows in p.lmm and p.plm.")
        }

        # Get the number of hypothesis tests
        n <- nrow(p)  # Variable p is the data frame of all p-values.
        p.adj.max <- p.adj.max / n  # take Bonferroni correction into account
        print(paste(n, "hypothesis tests have been performed for each kind of models.", sep = " "))

        # Replace p-values that are too small to be accurate with the minimum acceptable p-value p.min
        p$p_wald_acc <- p$p_wald  # make a column for "accurate" p-values
        p$p_wald_acc[p$p_wald_acc < p.min] <- p.min
        p$p_chisq_acc <- p$p_chisq
        p$p_chisq_acc[p$p_chisq_acc < p.min] <- p.min

        # Negative-log10 transformation of remaining raw p-values
        p$p_wald_nevLog10 <- -log10(p$p_wald_acc)  # for accurate p-values from LMMs
        p$p_chisq_nevLog10 <- -log10(p$p_chisq_acc)  # for accurate p-values from PLMs

        # Grouping data points
        group.by.lambda0 <- (L.weak > 0) && (!is.null(lmm.h0))
        if (group.by.lambda0) {
            # Retrieve lambda0 and perform Log10 transformation to lambdas
            p$L0_x <- lmm.h0$lambda0_remle[match(p$x_pat, lmm.h0$y_pat)]
            p$L0_y <- lmm.h0$lambda0_remle[match(p$y_pat, lmm.h0$y_pat)]

            # Group data points (pattern pairs) based on their lambda0 values
            p$L0_x_group <- as.integer(p$L0_x > L.weak)  # 1: strong structural random effect in X; otherwise, group code = 0.
            p$L0_y_group <- as.integer(p$L0_y > L.weak) * 2  # 0: weak structural random effect; otherwise, group code = 2.
            p$group <- p$L0_x_group + p$L0_y_group  # 0, 1, 2, 3 (strong structural random effect in both X and Y)
            p$group <- as.factor(p$group)  # convert group codes to levels: 0, 1, 2, 3

            # Get the final data frame, which contains the group codes
            p <- p[, c("y_pat", "x_pat", "p_wald", "p_chisq", "p_wald_acc", "p_chisq_acc",
                       "p_wald_nevLog10", "p_chisq_nevLog10", "L0_y", "L0_x", "group")]
        } else {  # grouping by raw p-values
            p$group <- as.integer(p$p_chisq_acc <= p.adj.max) + as.integer(p$p_wald_acc <= p.adj.max) * 2  # possible values: 0, 1, 2, 3
            p$group <- as.factor(p$group)
        }
    } else {
        p <- p.prev
    }

    # Calculate an upper bound of p-values for significance ===============
    bks.max <- max(bks)  # the maximum of break points
    p.adj.max <- -log10(p.adj.max)  # -log10 transformed cut-off of raw p-values
    print(paste("Minimum -log10 p-value for significance in comparable tests:",
                round(p.adj.max, digits = 4), sep = " "))

    # Prepare three ggplot objects that constitute the resulting figure ===============
    group.names <- names(cols)
    three.panels <- panel.num == 3

    if (three.panels) {
        # Panel a: a box plot for coordinates on the Y axis (p-values from LMMs)
        if (boxplot.show.group) {
            a <- ggplot(data = p, mapping = aes(x = group, y = p_wald_nevLog10, fill = group)) +
                geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black") +
                scale_fill_manual(breaks = group.names, values = cols)
            lab.x <- "Group"
        } else {
            a <- ggplot(data = p, mapping = aes(x = "", y = p_wald_nevLog10)) +
                geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black", fill = boxplot.default.fill)
            lab.x <- ""
        }

        if (show.p.adj.max) {
            a <- a + geom_hline(yintercept = p.adj.max, linetype = "dashed", size = 0.5, colour = "grey25")
        }

        a <- a + labs(title = panel.titles[1], x = lab.x, y = expression(-log[10](p[LMM]))) +
            scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            theme_bw() +
            theme(plot.title = element_text(face = "bold", size = plot.title.size),
                  panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey95"),
                  axis.text = element_text(size = axis.text.size, colour = "black"),
                  axis.title = element_text(size = axis.title.size, colour = "black"),
                  legend.position = "none")

        # Panel b: a scatter plot showing how raw p-values differ between the two kinds of models
        b <- ggplot(data = p, mapping = aes(x = p_chisq_nevLog10, y = p_wald_nevLog10, colour = group)) +
            geom_point(shape = 16, size = 1)

        if (show.p.adj.max) {
            b <- b + geom_hline(yintercept = p.adj.max, colour = "grey25", size = 0.5, linetype = "dashed") +
                geom_vline(xintercept = p.adj.max, colour = "grey25", size = 0.5, linetype = "dashed")
        }

        b <- b + geom_abline(intercept = 0, slope = 1, size = 0.5, colour = "black") +
            labs(title = panel.titles[2], x = expression(-log[10](p[PLM])), y = expression(-log[10](p[LMM]))) +
            scale_colour_manual(breaks = group.names, values = cols) +
            scale_x_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            theme_bw() +
            theme(plot.title = element_text(face = "bold", size = plot.title.size),
                  panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey95"),
                  axis.text = element_text(size = axis.text.size, colour = "black"),
                  axis.title = element_text(size = axis.title.size.panelB, face = "bold", colour = "white"),
                  legend.position = "none")  # Set the text colour of axis titles to white to hide them

        # Panel c: a box plot for coordinates on the X axis (p-values from PLMs)
        if (boxplot.show.group) {
            c <- ggplot(data = p, mapping = aes(x = group, y = p_chisq_nevLog10, fill = group)) +
                geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black") +
                scale_fill_manual(breaks = group.names, values = cols)
        } else {
            c <- ggplot(data = p, mapping = aes(x = "", y = p_chisq_nevLog10)) +
                geom_boxplot(lwd = 0.5, outlier.size = 0.8, colour = "black", fill = boxplot.default.fill)
        }

        if (show.p.adj.max) {
            c <- c + geom_hline(yintercept = p.adj.max, linetype = "dashed", size = 0.5, colour = "grey25")
        }

        c <- c + labs(title = panel.titles[3], x = lab.x, y = expression(-log[10](p[PLM]))) +
            scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            coord_flip() + theme_bw() +
            theme(plot.title = element_text(face = "bold", size = plot.title.size),
                  panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey95"),
                  axis.text = element_text(size = axis.text.size, colour = "black"),
                  axis.title = element_text(size = axis.title.size, colour = "black"),
                  legend.position = "none")
    } else {
        # Panel a: a scatter plot showing how raw p-values differ between the two kinds of models
        b <- ggplot(data = p, mapping = aes(x = p_chisq_nevLog10, y = p_wald_nevLog10, colour = group)) +
            geom_point(shape = 16, size = 1)

        if (show.p.adj.max) {
            b <- b + geom_hline(yintercept = p.adj.max, colour = "grey25", size = 0.5, linetype = "dashed") +
                geom_vline(xintercept = p.adj.max, colour = "grey25", size = 0.5, linetype = "dashed")
        }

        b <- b + geom_abline(intercept = 0, slope = 1, size = 0.5, colour = "black") +
            labs(title = panel.titles[1], x = expression(-log[10](p[PLM])), y = expression(-log[10](p[LMM]))) +
            scale_colour_manual(breaks = group.names, values = cols) +
            scale_x_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            theme_bw() +
            theme(plot.title = element_text(face = "bold", size = plot.title.size),
                  panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey95"),
                  axis.text = element_text(size = axis.text.size, colour = "black"),
                  axis.title = element_text(size = axis.title.size.panelB, face = "bold", colour = "black"),
                  legend.position = "none")

        # Panel b: a box plot with two boxes comparing p-values from PLMs and LMMs
        require(reshape2)

        D <- p[, c("p_chisq_nevLog10", "p_wald_nevLog10")]
        names(D) <- c("PLM", "LMM")
        D <- melt(data = D, measure.vars = c("PLM", "LMM"), variable.name = "Group", value.name = "p",
                  factorsAsStrings = FALSE)
        d <- ggplot(data = D) + geom_boxplot(mapping = aes(x = Group, y = p), fill = boxplot.default.fill,
                                             lwd = 0.5, outlier.size = 0.8, colour = "black")
        if (show.p.adj.max) {
            d <- d + geom_hline(yintercept = p.adj.max, linetype = "dashed", size = 0.5, colour = "grey25")
        }

        d <- d + labs(title = panel.titles[2], x = "Model", y = expression(-log[10](p))) +
            scale_y_continuous(limits = c(0, bks.max), trans = "sqrt", breaks = bks,
                               labels = as.character(bks), expand = c(0, 0)) +
            coord_flip() + theme_bw() +
            theme(plot.title = element_text(face = "bold", size = plot.title.size),
                  panel.grid.major = element_line(colour = "grey90"),
                  panel.grid.minor = element_line(colour = "grey95"),
                  axis.text = element_text(size = axis.text.size, colour = "black"),
                  axis.title = element_text(size = axis.title.size, colour = "black"),
                  legend.position = "none")
    }

    # Draw a figure ===============
    png(filename = img, width = img.w, height = img.h, units = img.u, res = img.r)
    par(oma = img.oma, mar = img.mar, mgp = img.mgp)
    if (three.panels) {
        if (boxplot.show.group) {
            grid.arrange(a, b, textGrob(""), c, nrow = 2, ncol = 2, widths = c(2, 5), heights = c(5, 2))
        } else {
            grid.arrange(a, b, textGrob(""), c, nrow = 2, ncol = 2, widths = c(1, 4), heights = c(4, 1))
        }
    } else {
        grid.arrange(b, d, nrow = 2, ncol = 1, heights = c(5, 2))
    }
    dev.off()

    return(p)
}
