#' @title Creating a temporal network from graphs
#'
#' @description This function creates a temporal network for visualisation from
#' a list of Graph objects. Graphs in this list are named by years. Every member
#' graph must consist of the same vertices and edges. In particular, when the
#' network is undirectional, node order of each edge must be the same in the data
#' frame defining edges.
#'
#' Dependency: package network and networkDynamic. The package ndtv is required
#' to display the output dynamic network as an animation.
#'
#' @param gs A GraphSet object with years as graph IDs.
#' @param v.label A single name or index for the column from which vertex labels
#' are drawn.
#' @param e.tail Column name for tail vertices.
#' @param e.head Column name for head vertices.
#' @param directed A logical parameter specifying whether the network is directed
#' or not. Default: FALSE (undirected).
#' @param t.gap An integer or numeric argument determining the non-inclusive
#' minimum difference (namely, t > t.gap rather than t >= t.gap) between two
#' consecutive time points to be considered as not continuous. This arguments
#' determines each onset and terminus range for edges.
#' @param v.value A character or integer specifying the column in V(gs) used as
#' the variable defining vertex sizes.
#' @param v.value.base An integer specifying the minimum value for the variable
#' (column) defining vertex sizes. For instance, you may set v.size.base for a
#' dynamic co-occurrence count where the expected minimum count is zero.
#' @param v.size.min An integer for the minimum vertex size (when v.size.base is
#' reached).
#' @param v.size.max An integer for the maximum vertex size.
#' @param v.colour A character vector of colours named by vertex labels.
#' @param e.weight Column name or index in E(gs[[t]]) for edge widths.
#' @param e.weight.base The minimum valid value of edge weights. It works in the
#' same way as v.value.base.
#' @param e.weight.cutoff Non-inclusive minimum edge weight that determines whether an edge
#' is present at a time point.
#' @param e.width.min Minimum edge width.
#' @param e.width.max Maximum edge width.
#' @param e.colour A column name or index for the attributed used for colour assignment.
#' @param e.colour.low A character specifying the colour for edges of the least weight.
#' @param e.colour.high A character specifying the colour for edges of the highest weight.
#' @param e.colour.num An integer specifying the number of colours to be assigned.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export tempNet
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 14 August 2017; the lastest edition: 18 August 2018

tempNet <- function(gs, v.label = "allele", e.tail = "node1", e.head = "node2",
                    directed = FALSE, t.gap = 1,
                    v.value = "count", v.value.base = 0, v.size.min = 1,
                    v.size.max = 25, v.colour = "red",
                    e.weight = "n_xy", e.weight.base = 1, e.weight.cutoff = 0,
                    e.width.min = 1, e.width.max = 10, e.colour = "m",
                    e.colour.low = "#FFA0A0", e.colour.high = "#FF0000",
                    e.colour.num = 9) {
    require(network)
    require(networkDynamic)

    # Validity assessment ###############
    if (class(gs) != "list" || class(gs[[1]]) != "Graph") {
        stop("Argument error: gs must be a list of Graphs.")
    }

    # Get time points from the names ###############
    ts <- sort(as.integer(names(gs), decreasing = FALSE))
    t.start <- min(ts)
    t.end <- max(ts)
    print(paste(length(ts), "time points to display.", sep = " "))
    print(paste0("Start and end time: ", t.start, " and ", t.end, "."))

    # Establish a HASH table for vertex IDs ###############
    # Because the networkDynamic only accepts numeric vertex IDs.
    # Here, I assume every Graph object has the same vertices.
    v.names <- sort(Vn(gs[[1]]), decreasing = FALSE)  # vertex names
    nv <- length(v.names)  # node number per graph
    vid.mapping <- data.frame(id = 1 : nv, label = v.names, stringsAsFactors = FALSE)
    if (!is.null(v.colour)) {
        if (length(v.colour) > 1) {
            vid.mapping$colour <- as.character(v.colour[vid.mapping$label])
        } else {
            print("A uniform colour is applied to vertices.")
            vid.mapping$colour <- rep(v.colour, times = nv)
        }
    } else {
        print("No vertex colour is specified. Use the default colour grey50.")
        vid.mapping$colour <- rep("grey50", times = nv)
    }


    # Pull edge and vertex lists together ###############
    # Assume all edges are placed in the same order in each list.
    E.attr <- data.frame(time = integer(0), tail = integer(0), head = integer(0),
                         weight = numeric(0), colour_attr = numeric(0),
                         edge.id = integer(0),
                         stringsAsFactors = FALSE)
    V.attr <- data.frame(time = integer(0), label = character(0), value = numeric(0),
                         stringsAsFactors = FALSE)
    ne <- nE(gs[[1]])  # edge number per graph and network slice
    print("Pulling edge and vertex information together.")
    for (t in ts) {
        # edge attributes | year = t
        g <- gs[[as.character(t)]]  # g is a Graph object.
        e <- E(g)
        e <- e[, c(e.tail, e.head, e.weight, e.colour)]  # edge attributes
        names(e) <- c("tail", "head", "weight", "colour_attr")
        e.t <- cbind.data.frame(time = rep(t, times = ne), e, edge.id = 1 : ne)
        E.attr <- rbind.data.frame(E.attr, e.t, stringsAsFactors = FALSE)  # edge.id is a reserved column name for edge IDs in the networkDynamic package.

        # vertex attributes | year = t
        v <- V(g)  # temporal vertex attributes
        v <- v[, c(v.label, v.value)]
        names(v) <- c("label", "value")
        v.t <- cbind.data.frame(time = rep(t, times = nrow(v)), v)
        V.attr <- rbind.data.frame(V.attr, v.t, stringsAsFactors = FALSE)
    }
    E.attr <- subset(E.attr, weight > e.weight.cutoff)

    # Assign edge colours
    colourGenerator <- colorRampPalette(colors = c(e.colour.low, e.colour.high),
                                        space = "rgb")
    attr.levels <- cut(E.attr$colour_attr,
                       breaks = seq(from = min(E.attr$colour_attr),
                                    to = max(E.attr$colour_attr),
                                    length.out = e.colour.num + 1),
                       include.lowest = TRUE)
    E.attr$colour <- colourGenerator(e.colour.num)[attr.levels]

    # Replace vertex names with integer IDs
    print("Converting vertex names to integer IDs.")
    E.attr$tail <- vid.mapping$id[match(E.attr$tail, vid.mapping$label)]
    E.attr$head <- vid.mapping$id[match(E.attr$head, vid.mapping$label)]
    V.attr$vertex.id <- vid.mapping$id[match(V.attr$label, vid.mapping$label)]  # five columns: time, label, value, colour, vertex.id

    # Linearly convert vertex values and edge weights into diameters and widths, respectively
    print("Performing linear transformation to vertex values and edge weights.")
    V.attr$size <- .linearTrans(x = V.attr$value, x0 = v.value.base,
                                y.min = v.size.min, y.max = v.size.max,
                                x.outlier = 0, y.outlier = 0)
    E.attr$width <- .linearTrans(x = E.attr$weight, x0 = e.weight.base,
                                 y.min = e.width.min, y.max = e.width.max,
                                 x.outlier = 0, y.outlier = 0)  # In fact, weights in E.attr should not contain any outliers at this stage.

    # Establish edge spells (activity script) from the pooled edge list
    act <- .mkEdgeSpells(E.attr, t.gap)

    # Create the temporal network, assuming all Graph objects contain the same
    # number of vertices.
    print("Initialising the temporal network.")
    net <- network.initialize(n = nv, directed = directed)
    net <- networkDynamic(base.net = net, edge.spells = act, start = t.start,
                          end = t.end)  # The edge spells define the range of time for visualisation in the dynamic network.

    # Append vertex attributes to the dynamic network ###############
    # Static attributes
    print("Appending vertex attributes.")
    for (i in 1 : nrow(vid.mapping)) {
        r <- vid.mapping[i, ]
        v <- r$id
        # Since the terminus is not included for appending the vertex attribute,
        # the argument terminus should be either t.end + 1 or Inf.
        activate.vertex.attribute(x = net, prefix = "label", value = r$label,
                                  v = v, onset = t.start, terminus = Inf)
        activate.vertex.attribute(x = net, prefix = "colour", value = r$colour,
                                  v = v, onset = t.start, terminus = Inf)
    }

    # Dynamic attribute: size
    for (t in ts) {  # Notice ts may have gaps.
        Vs.t <- subset(V.attr, time == t)
        for (i in 1 : nrow(Vs.t)) {
            r <- Vs.t[i, ]
            activate.vertex.attribute(x = net, prefix = "size", value = r$size,
                                      at = t, v = r$vertex.id)
        }
    }
    for (t in setdiff(x = t.start : t.end, y = ts)) {  # fill the gaps
        activate.vertex.attribute(x = net, prefix = "size", value = 0, at = t)  # for all vertices in the current year
    }

    # Add edge attributes to the network ###############
    print("Appending edge attributes. It may take a while to finish.")
    for (t in ts) {
        Es.t <- subset(E.attr, time = t)
        for (i in 1 : nrow(Es.t)) {
            r <- Es.t[i, ]
            activate.edge.attribute(x = net, prefix = "width", value = r$width,
                                    at = t, e = r$edge.id)
            activate.edge.attribute(x = net, prefix = "colour", value = r$colour,
                                    at = t, e = r$edge.id)
        }
    }

    # Return a temporal graph and the action table
    return(list(G = net, A = act, E = E.attr, V = V.attr, M = vid.mapping))
}

.mkEdgeSpells <- function(Es, t.gap) {
    print("Generating an activity list for edges.")
    act <- data.frame(onset = integer(0), terminus = integer(0),
                      tail = integer(0), head = integer(0),
                      edge.id = integer(0), stringsAsFactors = FALSE)

    # Go through all edges
    # Assumes that at each time point, there is one and only one edge of the
    # same alleles present.
    eids <- sort(unique(Es$edge.id), decreasing = FALSE)
    for (i in eids) {
       Ei <- subset(Es, edge.id == i)  # all time points of the current edge
       if (nrow(Ei) > 0) {
           v.tail <- Ei$tail[1]  # tail vertex of the current edge
           v.head <- Ei$head[1]  # head vertex
           ts <- sort(unique(Ei$time), decreasing = FALSE)  # 1, 2, 3, 5, 6, 9, ...
           if (length(ts) > 1) {  # This edge appeared multiple times.
               t.on <- ts[1]  # the lower bound of an interval
               t.max <- max(ts)
               t0 <- t.on  # lagging pointer
               for (t in ts[-1]) {
                   if (t == t.max) {
                       if (t - t0 > t.gap) {  # In the meantime, the last time point separates from the others.
                           spell <- data.frame(onset = c(t.on, t), terminus = c(t0, t),
                                               tail = rep(v.tail, times = 2),
                                               head = rep(v.head, times = 2),
                                               edge.id = rep(i, times = 2), stringsAsFactors = FALSE)
                       } else {  # force the current interval to be saved
                           spell <- data.frame(onset = t.on, terminus = t,
                                               tail = v.tail, head = v.head,
                                               edge.id = i, stringsAsFactors = FALSE)
                       }
                       act <- rbind.data.frame(act, spell, stringsAsFactors = FALSE)  # Then the "for" loop terminates.
                   } else if (t - t0 > t.gap) {
                       # a gap is found and an interval is hence determined
                       # Now t is in the next interval.
                       spell <- data.frame(onset = t.on, terminus = t0,
                                           tail = v.tail, head = v.head,
                                           edge.id = i, stringsAsFactors = FALSE)
                       act <- rbind.data.frame(act, spell, stringsAsFactors = FALSE)
                       t.on <- t  # move the onset to the current point
                       t0 <- t
                   } else {  # The increment of time points remains stable.
                       t0 <- t  # move to the next time point
                   }
               }
           } else {  # This edge only appear once.
               spell <- data.frame(onset = ts, terminus = ts,
                                   tail = v.tail, head = v.head, edge.id = i,
                                   stringsAsFactors = FALSE)
           }
       } else {
           print(paste("Warning: no edge is present for the edge ID", i,
                       "in Es, which is abnormal. Skip this edge.", sep = " "))
       }
    }

    return(act)
}

.linearTrans <- function(x, x0, y.min, y.max, x.outlier = 0, y.outlier = 0.5) {
    # Apply to a single value x
    y.range <- y.max - y.min
    x.range <- max(x) - x0
    b <- y.range / x.range
    y <- sapply(x, function(z) ifelse(z > x.outlier, round((z - x0) * b + y.min, digits = 4), y.outlier))

    return(y)
}
