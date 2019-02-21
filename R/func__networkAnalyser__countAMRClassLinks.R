#' @title Count the number of links in a network at the level of antimicrobial classes
#'
#' @description This function is useful in summarising connections between
#' antimicrobial classes based on a network.
#'
#' @param net A Graph object, or a data frame whose first two columns are names
#' of allele Y (response variable) and X (explanatory variable).
#' @param mapping A data frame mapping allele names to antimicrobial classes. It
#' can be the data frame "mapping" in the output list of findPhysLink.
#' @param out.matrix A logical argument determining whether a square matrix or a
#' data frame is returned.
#'
#' @return A data frame of three columns: Class_y and Class_x, the target and
#' source antimicrobial class, respectively; Links, link counts per kind of edge
#' (Class_y, Class_x).
#'
#' @export countAMRClassLinks
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#
# First and the latest edition: 21 Feb 2019

countAMRClassLinks <- function(net, mapping, out.matrix = FALSE) {
    c <- class(net)
    if (c == "Graph") {
        cl <- E(net)
        cl <- cl[, c(1, 2)]  # target and source alleles
    } else if (c == "data.frame") {
        cl <- net[, c(1, 2)]
    } else {
        stop("Argument error: net must be a GeneMate Graph object or a generic data frame.")
    }
    names(cl) <- c("Class_y", "Class_x")
    cl$Class_y <- mapping$class[match(cl$Class_y, mapping$allele)]
    cl$Class_x <- mapping$class[match(cl$Class_x, mapping$allele)]
    pat <- unique(cl)  # link patterns: unique rows of cl
    link_n <- mapply(.countLinkPatterns, pat$Class_x, pat$Class_y,
                     MoreArgs = list(links = cl), SIMPLIFY = FALSE)
    link_n <- do.call("rbind.data.frame", args = link_n)

    # Convert a data frame into a matrix
    if (out.matrix) {
        classes <- sort(union(pat$Class_y, pat$Class_x), decreasing = FALSE)
        n <- length(classes)
        M <- matrix(data = 0, nrow = n, ncol = n, dimnames = list(classes, classes))  # row: source; column: target
        for (i in 1 : nrow(link_n)) {
            r <- link_n[i, ]
            M[r[["Class_x"]], r[["Class_y"]]] <- r[["Links"]]
        }
        link_n <- M
    }

    return(link_n)
}

.countLinkPatterns <- function(x, y, links) {
    # a subordinate function of countAMRClassLinks
    n <- sum(links$Class_y == y & links$Class_x == x)

    return(data.frame(Class_y = y, Class_x = x, Links = n, stringsAsFactors = FALSE))
}
