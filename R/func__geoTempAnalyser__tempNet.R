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
#' @param directed A logical parameter specifying whether the network is directed
#' or not. Default: FALSE (undirected).
#' @param tail_col Column name for tail vertices.
#' @param head_col Column name for head vertices.
#' @param t_gap An integer or numeric argument determining the non-inclusive
#' minimum difference (namely, t > t_gap rather than t >= t_gap) between two
#' consecutive time points to be considered as not continuous.
#' @param e_col Column name for edge weights.
#' @param e_width Non-inclusive minimum edge weight that determines whether an edge
#' is present at a time point.
#' @param v_col A single or multiple names for the columns from which temporal
#' vertex attributes are drawn and incorporated into the network. However, since
#' this function reserves "label" as an attribute, this name must not present in
#' this argument. Furthermore, va_col must not point to the first column, which
#' is used for vertex labels.
#' @param default_va A named list of default attribute values for va_col.
#'
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#' @export tempNet
#
#  Copyright 2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 14 August 2017; the lastest edition: 15 August 2018

tempNet <- function(gs, tail_col = "node1", head_col = "node2", directed = FALSE,
                    t_gap = 1, t_min = 1, t_max = 10,
                    v_size = NULL, v_col = NULL,
                    e_width = NULL, e_col = NULL,
                    default_va = list()) {
    require(network)
    require(networkDynamic)

    # Validity assessment
    if (class(gs) != "list" || class(gs[[1]]) != "Graph") {
        stop("Argument error: gs must be a list of Graphs.")
    }
    if (length(va_col) != length(default_va) || !is.list(default_va) || any(! names(default_va) %in% va_col)) {
        stop("Argument error: default_va does not match to va_col.")
    }

    # Get time points from the names
    ts <- sort(as.integer(names(gs), decreasing = FALSE))
    print(paste(length(ts), "time points to display.", sep = " "))

    # Establish a HASH table for vertex IDs as networkDynamic only accepts numeric
    # vertex IDs.
    v_names <- sort(Vn(gs[[1]]), decreasing = FALSE)  # vertex names
    nv <- length(v_names)  # node number per graph
    Vids <- data.frame(id = 1 : nv, name = v_names, stringsAsFactors = FALSE)

    # Pull edge and vertex lists together, assuming all edges are placed in the
    # same order in each list.
    Es <- NULL
    Vs <- NULL
    ne <- nE(gs[[1]])  # edge number per graph
    eids <- 1 : ne  # edge indices
    print("Pulling edge and vertex information together.")
    for (t in ts) {
        g <- gs[[as.character(t)]]  # g is a Graph object.
        e <- E(g)
        e <- e[, c(tail_col, head_col, ew_col)]  # edge list
        names(e) <- c("tail", "head", "weight")
        e_t <- cbind.data.frame(time = rep(t, times = ne), e, edge.id = eids)
        Es <- rbind.data.frame(Es, e_t, stringsAsFactors = FALSE)  # edge.id is a reserved column name for edge IDs in the networkDynamic package.
        v <- V(g)  # temporal vertex attributes
        v_attr <- !is.null(va_col)  # whether there are temporal vertex attributes or not
        if (v_attr) {
            v <- v[, c(1, which(names(v) %in% va_col))]
            names(v)[1] <- "label"
        } else {  # no temporal attribute to be drawn
            v <- data.frame(label = v[, 1], stringsAsFactors = FALSE)
        }
        v_t <- cbind.data.frame(time = rep(t, times = nrow(v)), v)
        Vs <- rbind.data.frame(Vs, v_t, stringsAsFactors = FALSE)
    }

    # Replace vertex names with integer IDs
    print("Converting vertex names to integer IDs.")
    Es$tail <- Vids$id[match(Es$tail, Vids$name)]
    Es$head <- Vids$id[match(Es$head, Vids$name)]
    Vs$vertex <- Vids$id[match(Vs$label, Vids$name)]
    if (v_attr) {
        Vs <- Vs[, c("time", "vertex", "label", va_col)]
    } else {
        Vs <- Vs[, c("time", "vertex", "label")]
    }

    # Establish edge spells (activity script) from the pooled edge list
    Es <- subset(Es, weight > ew_min)
    ts <- sort(unique(Es$time), decreasing = FALSE)  # remaining time points
    print("Generating an activity list for edges.")
    act <- .mkEdgeSpells(Es, t_gap)

    # Create the temporal network, assuming all Graph objects contain the same
    # number of vertices.
    print("Initialising the temporal network.")
    net <- network.initialize(n = nv, directed = directed)
    net <- networkDynamic(base.net = net, edge.spells = act)

    # Append labels to vertices
    print("Appending vertex labels.")
    activate.vertex.attribute(x = net, prefix = "label", value = "",
                              onset = -Inf, terminus = Inf)  # null vertex label
    for (v in unique(Vs$vertex)) {
        activate.vertex.attribute(x = net, prefix = "label",
                                  value = Vids$name[Vids$id == v],
                                  v = v, onset = -Inf, terminus = Inf)
    }

    # Append temporal vertex attributes to the vertices
    if (v_attr) {
        print("Appending temporal vertex attributes.")

        # for other attributes
        for (va in va_col) {  # for each vertex attribute
            # set default values
            def_va <- default_va[[va]]
            for (t in t_min : t_max) {
                activate.vertex.attribute(x = net, prefix = va, value = def_va, at = t)
                activate.vertex.attribute(x = net, prefix = va, value = def_va, at = t)
            }

            # set individual values
            for (t in ts) {
                Vs_t <- subset(Vs, time == t)
                n <- nrow(Vs_t)
                if (n > 0) {
                    for (i in 1 : n) {
                        attr_t_i <- Vs_t[i, ]
                        activate.vertex.attribute(x = net, prefix = va,
                                                  value = attr_t_i[[va]],
                                                  at = t, v = attr_t_i[["vertex"]])
                    }
                } else {
                    print(paste0("There is no value for the temporal attribute ",
                                va, " at the time ", t, ". Skip this year."))
                }
            }
        }
    } else {
        print("No temporal vertex attributes to be attached.")
    }

    # Return a temporal graph and the action table
    # M: mapping
    return(list(G = net, A = act, E = Es, V = Vs, M = Vids))
}

.mkEdgeSpells <- function(Es, t_gap) {
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
           v_tail <- Ei$tail[1]  # tail vertex of the current edge
           v_head <- Ei$head[1]  # head vertex
           ts <- sort(unique(Ei$time), decreasing = FALSE)  # 1, 2, 3, 5, 6, 9, ...
           if (length(ts) > 1) {  # This edge appeared multiple times.
               t_on <- ts[1]  # the lower bound of an interval
               t_max <- max(ts)
               t0 <- t_on  # lagging pointer
               for (t in ts[-1]) {
                   if (t == t_max) {
                       if (t - t0 > t_gap) {  # In the meantime, the last time point separates from the others.
                           spell <- data.frame(onset = c(t_on, t), terminus = c(t0, t),
                                               tail = rep(v_tail, times = 2),
                                               head = rep(v_head, times = 2),
                                               edge.id = rep(i, times = 2), stringsAsFactors = FALSE)
                       } else {  # force the current interval to be saved
                           spell <- data.frame(onset = t_on, terminus = t,
                                               tail = v_tail, head = v_head,
                                               edge.id = i, stringsAsFactors = FALSE)
                       }
                       act <- rbind.data.frame(act, spell, stringsAsFactors = FALSE)  # Then the "for" loop terminates.
                   } else if (t - t0 > t_gap) {
                       # a gap is found and an interval is hence determined
                       # Now t is in the next interval.
                       spell <- data.frame(onset = t_on, terminus = t0,
                                           tail = v_tail, head = v_head,
                                           edge.id = i, stringsAsFactors = FALSE)
                       act <- rbind.data.frame(act, spell, stringsAsFactors = FALSE)
                       t_on <- t  # move the onset to the current point
                       t0 <- t
                   } else {  # The increment of time points remains stable.
                       t0 <- t  # move to the next time point
                   }
               }
           } else {  # This edge only appear once.
               spell <- data.frame(onset = ts, terminus = ts,
                                   tail = v_tail, head = v_head, edge.id = i,
                                   stringsAsFactors = FALSE)
           }
       } else {
           print(paste("Warning: no edge is present for the edge ID", i,
                       "in Es, which is abnormal. Skip this edge.", sep = " "))
       }
    }

    return(act)
}
