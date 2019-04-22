#' @title An expansion of ggtree's gheatmap function for displaying an allelic
#' presence-absence matrix.
#'
#' @description To append a heatmap of a presence-absence matrix (PAM) to the
#' right side of a phylogenetic tree. This function is modified from the function
#' gheatmap in the R package ggtree (github.com/GuangchuangYu/ggtree) through
#' replacing scale_fill_gradient with scale_fill_discrete. As a result, columns
#' in the heat map can be coloured by group information. This function is useful
#' when users need draw a heat map of gene or allele contents across strains.
#'
#' In the PAM, absence status is coded as zero and presence status is coded as one.
#'
#' @param p A tree view from the function ggtree in the package ggtree
#' @param data A binary presence-absence matrix (PAM)
#' @param col_colours Either a vector of colour codes named by column names of
#' the PAM, or a single colour code for all columns.
#' @param null_colour The colour for absence status.
#' @param border_colour Colour of heatmap cell border
#' @param cluster_cols A logical argument determining if columns of PAM will be
#' clustered or not. Default: FALSE. Clustering will not run when there are less
#' than three columns in the PAM.
#' @param cluster_method Method for clustering. Default: binary.
#' @param cluster_distance The distance used for clustering columns. Default: binary.
#' @param rev_cols A logical argument specifying whether to reverse the column order
#' when it is necessary for a better visualisation.
#' @param colnames Logical, add matrix colnames or not
#' @param colnames_position One of 'bottom' or 'top'
#' @param colnames_angle Angle of column names
#' @param colnames_level Levels of colnames
#' @param set_label_colours Whether colour column names in the same way as painting columns.
#' By default, this option is turned off, which result in all column names to be printed
#' in black.
#' @param colnames_offset_x x offset for column names
#' @param colnames_offset_y y offset for column names
#' @param font.size Font size of matrix colnames
#' @param hjust hjust for column names (0: align left, 0.5: align center, 1: align righ)
#' @param offset Offset of heatmap to tree
#' @param width Total width of heatmap, compare to width of tree
#' @param show_legend A logical argument determining whether to show the legend.
#'
#' @return A list of three elements: p, a tree view that can be plotted directly;
#' mapping: a data frame mapping allele names to colours; data, a data frame of
#' weighted presence-absence status.
#'
#' @examples htmap <- heatMapPAM(...)
#' plot(htmap$p)
#'
#' @export
#' @author Guangchuang Yu, Yu Wan (\email{wanyuac@@gmail.com})
#
# First edition of this function: 5 Sep 2018; the latest edition: 6 Feb 2019
# Licence: Artistic License 2.0 (follow the licence of the package ggtree)

heatMapPAM <- function(p, data, col_colours = "black", null_colour = "grey90",
                       border_colour = "white",
                       cluster_cols = FALSE, cluster_method = "single",
                       cluster_distance = "binary", rev_cols = FALSE,
                       colnames = TRUE, colnames_position = "bottom", colnames_angle = 0,
                       colnames_level = NULL, set_label_colours = FALSE,
                       colnames_offset_x = 0, colnames_offset_y = 0,
                       font.size = 4, hjust = 0.5, offset = 0, width = 1,
                       show_legend = FALSE) {
    # The first two packages are dependencies of the package ggtree.
    require(ggplot2)
    require(tidyr)  # for the function gather
    require(tibble)
    require(magrittr)  # for operators "%<>%" and "%>%"(github.com/GuangchuangYu/ggtree/blob/master/R/operator.R)
    require(ggtree)

    colnames_position %<>% match.arg(c("bottom", "top"))
    variable <- value <- lab <- y <- NULL

    ## convert relative width of the PAM to width of each cell
    width <- width * (p$data$x %>% range(na.rm = TRUE) %>% diff) / ncol(data)
    isTip <- x <- y <- variable <- value <- from <- to <- NULL

    df <- p$data
    df <- df[df$isTip, ]
    start <- max(df$x, na.rm = TRUE) + offset

    # Column-wise clustering
    if (cluster_cols && ncol(data) > 2) {
        hc <- hclust(d = dist(t(data), method = cluster_distance),
                     method = cluster_method)
        data <- data[, hc$order]  # reorder columns
        clustered <- TRUE
    } else {
        clustered <- FALSE
    }

    # Inverse columns for improving visualisation
    if (rev_cols) {
        data <- data[, ncol(data) : 1]
    }

    # weight cells before converting the PAM (dd) into a data frame
    n_colours <- length(col_colours)
    is_multi_colour <- n_colours > 1 && !is.null(names(col_colours))
    if (is_multi_colour) {
        colours_uniq <- sort(unique(as.character(col_colours)), decreasing = FALSE)  # colours for positive values in PAM
        colour_codes <- 1 : length(colours_uniq)  # codes for colours
        names(colour_codes) <- colours_uniq  # e.g., c("red" = 1, "blue" = 2, ...)
        column_names <- colnames(data)  # The matrix product loses column names, that is, the allele names.
        if (clustered) {
            col_colours <- col_colours[column_names]
        } else if (n_colours > length(column_names))  {
            col_colours <- col_colours[column_names]
        }
        data <- data %*% diag(as.integer(colour_codes[as.character(col_colours)]))  # convert colour characters into integer codes
        colnames(data) <- column_names
        names(colours_uniq) <- as.character(colour_codes)
        colours_uniq <- append(c("0" = null_colour), colours_uniq)

        # Convert values in the matrix into factors so as to map colours to the levels
        dd <- as.data.frame(data)
        for (i in names(dd)) {
            dd[, i] <- factor(dd[, i], levels = c(0, as.integer(colour_codes)),
                              labels = c("0", as.character(colour_codes)))
        }
    } else {
        dd <- as.data.frame(data)
    }

    i <- order(df$y)

    ## handle collapsed tree
    ## https://github.com/GuangchuangYu/ggtree/issues/137
    i <- i[!is.na(df$y[i])]

    lab <- df$label[i]
    ## https://github.com/GuangchuangYu/ggtree/issues/182
    dd <- dd[match(lab, rownames(dd)), , drop = FALSE]
    dd$y <- sort(df$y)
    dd$lab <- lab

    # tibble::set_tidy_names(dd) solves the problem of the error: "Can't bind data because some arguments have the same name"
    # see https://github.com/tidyverse/tidyr/issues/472
    dd <- gather(data = tibble::set_tidy_names(dd), key = variable, value = value, -c(lab, y))

    i <- which(dd$value == "")
    if (length(i) > 0) {
        dd$value[i] <- NA
    }
    if (is.null(colnames_level)) {
        dd$variable <- factor(dd$variable, levels = colnames(data))
    } else {
        dd$variable <- factor(dd$variable, levels = colnames_level)
    }
    V2 <- start + as.numeric(dd$variable) * width

    # Create a data frame for label attributes
    mapping <- data.frame(from = as.character(dd$variable), to = V2, stringsAsFactors = FALSE)  # from: label texts
    mapping <- unique(mapping)
    paint_labels <- set_label_colours && is_multi_colour
    if (paint_labels) {
        mapping$col <- as.character(col_colours[mapping$from])  # variable label colour
        mapping$col <- factor(mapping$col, levels = as.character(colours_uniq),
                              labels = names(colours_uniq))
    }
    mapping$from = as.factor(mapping$from)  # change back to factors

    # Create the coloured tile matrix
    dd$x <- V2
    dd$width <- width
    dd[[".panel"]] <- factor("Tree")
    if (is.null(border_colour)) {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                            inherit.aes = FALSE)
    } else {
        p2 <- p + geom_tile(data = dd, aes(x, y, fill = value), width = width,
                            colour = border_colour, inherit.aes = FALSE)
    }

    # Print column names
    if (colnames) {
        if (colnames_position == "bottom") {
            y <- 0
        } else {
            y <- max(p$data$y) + 1
        }
        mapping$y <- y
        mapping[[".panel"]] <- factor("Tree")
        if (paint_labels) {
            p2 <- p2 + geom_text(data = mapping,
                                 aes(x = to, y = y, label = from, colour = col),
                                 size = font.size, inherit.aes = FALSE,
                                 angle = colnames_angle, nudge_x = colnames_offset_x,
                                 nudge_y = colnames_offset_y, hjust = hjust)
        } else {  # Labels are printed in black.
            p2 <- p2 + geom_text(data = mapping,
                                 aes(x = to, y = y, label = from),
                                 size = font.size, inherit.aes = FALSE,
                                 angle = colnames_angle, nudge_x = colnames_offset_x,
                                 nudge_y = colnames_offset_y, hjust = hjust)
        }
    }

    # Scale colours for both tiles and column labels
    if (is_multi_colour) {
        # This command does not scale label colours when paint_labels = FALSE.
        # Class: title of the legend
        # scale_fill_manual: for tile colours; scale_colour_manual: for label colours
        p2 <- p2 + scale_fill_manual(name = "Class",
                                     breaks = c("0", as.character(colour_codes)),
                                     values = colours_uniq, na.value = NA) +
            scale_colour_manual(name = "Class",
                                breaks = c("0", as.character(colour_codes)),
                                values = colours_uniq, na.value = NA)
    } else {
        p2 <- p2 + scale_fill_gradient(low = null_colour, high = col_colours,
                                       na.value = NA)
    }

    # Print legend
    if (show_legend) {
        p2 <- p2 + theme(legend.position = "right", legend.title = element_blank())
    }

    attr(p2, "mapping") <- mapping

    return(list(p = p2, mapping = mapping, data = dd))
}
