#' Draw a heat map for a given variable
#'
#' This is a generic function for plotting a heat map given two ID variables and a value variable.
#'
#' @param data A data frame of at least three columns, such as x, y and z, or a matrix to be plotted directly as a heat map.
#' @param x Name (character) of the first variable of IDs.
#' @param y Name (character) of the second variable of IDs.
#' @param val Name (character) or index (integer) of the variable for values in the heat matrix.
#' @param diag The value of diagonal cells in the unclustered heat map.
#' @param replace.na The value used to replace NA's.
#' @param cluster.row Whether cluster rows hierarchically.
#' @param cluster.col Whether cluster columns hierarchically.
#' @param dist.row The distance metric for hierarchical clustering of rows.
#' @param dist.col The distance metric for hierarchical clustering of columns.
#' @param cluster.method The method for hierarchical clustering.
#' @param colour.low The colour for the lowest value in the heat matrix.
#' @param colour.mid The colour for the middle (that is, the "zero") value. Set to NULL if only a single colour transition are to be included.
#' @param colour.high The colour for the highest value.
#' @param colour.grad The number of colour grades.
#' @param colour.breaks Break points (colour.grad + 1) of values for assigning colours.
#' @param display.val A boolean value determining whether to overlay values on the heat map or not.
#' @param font.size Font size of labels in the heat map.
#' @param filename Path and name for the output image. The filename extension (in lower case) determines the output format.
#' @param res Resolution of the output figure. Default: 72 ppi.
#' @param width Width of the output image.
#' @param height Height of the output image.
#' @param unit The unit of the width and height of the output image. Valid values: "mm" and "px" (default).
#'
#' @author Yu Wan (\email{wanyuac@gmail.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  Development history: 12 - 13 April 2017

drawHeatMap <- function(data = NULL, x = "x", y = "y", val = "val", diag = 1, replace.na = 0,
                        cluster.row = TRUE, cluster.col = TRUE, cluster.method = "complete",
                        dist.row = "euclidean", dist.col = "euclidean",
                        colour.low = "white", colour.mid = NULL, colour.high = "red", colour.grad = 50,
                        colour.breaks, display.val = FALSE, font.size = 10,
                        filename = "heatmap.png", res = 72, width = 800, height = 800, unit = "px") {
    # data preparation
    if (is.data.frame(data)) {
        # remove unnecessary columns from the data frame
        if (is.character(x)) {
            x <- .retriveIndex(data, x)
        }
        if (is.character(y)) {
            y <- .retriveIndex(data, y)
        }
        if (is.character(val)) {
            val <- .retriveIndex(data, val)
        }

        # convert the data frame into a matrix
        data <- .df2matrix(df = data[, c(x, y, val)], diag = diag, replace.na = replace.na)
    } else if (!is.matrix(data)) {
        stop("Argument error: the data must be either a data frame or a matrix.")
    }

    # initialise a heat map
    require(pheatmap)
    fn.len <- nchar(filename)
    fn.ext <- tolower(substr(filename, start = fn.len - 2, stop = fn.len))
    if (fn.ext == "png") {
        png(filename = filename, width = width, height = height, res = res, units = unit)
    } else {
        pdf(filename = filename, paper = "a4", width = width, height = height)
    }

    # set colours
    if (is.null(colour.mid)) {  # a single colour transition, such as white <-> red
        colGenerator <- colorRampPalette(c(colour.low, colour.high), space = "rgb")
    } else {  # two transitions, such as blue <-> white <-> red
        colGenerator <- colorRampPalette(c(colour.low, colour.mid, colour.high), space = "rgb")
    }
    colours <- colGenerator(n = colour.grad)  # evenly divide the whole range of colours

    # draw the heat map
    htmap <- pheatmap(mat = data, color = colours, breaks = colour.breaks, display_numbers = display.val, fontsize = font.size,
                      cluster_rows = cluster.row, cluster_cols = cluster.col, clustering_distance_rows = dist.row,
                      clustering_distance_cols = dist.col, clustering_method = cluster.method)

    dev.off()

    return(list(data = data, htmap = htmap))
}
