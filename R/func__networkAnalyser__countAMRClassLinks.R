#' @title Count the number of links in a network at the level of antimicrobial classes
#'
#' @description This function is useful in summarising connections between
#' antimicrobial classes based on a network.
#'
#' @param net A Graph object.
#' @param mapping A data frame mapping allele names to antimicrobial classes. It
#' can be the data frame "mapping" in the output list of findPhysLink.
#'
#' @return A data frame of three columns: Class_y and Class_x, the target and
#' source antimicrobial class, respectively; Links, link counts per kind of edge
#' (Class_y, Class_x).
#'
#' @export countAMRClassLinks
#' @author Yu Wan (\email{wanyuac@@gmail.com})
#
# First and the latest edition: 21 Feb 2019

countAMRClassLinks <- function(net, mapping) {
    cl <- E(net)
    cl <- cl[, c(1, 2)]  # target and source alleles
    names(cl) <- c("Class_y", "Class_x")
    cl$Class_y <- mapping$class[match(cl$Class_y, mapping$allele)]
    cl$Class_x <- mapping$class[match(cl$Class_x, mapping$allele)]
    pat <- unique(cl)  # link patterns: unique rows of cl
    link_n <- mapply(.countLinkPatterns, pat$Class_x, pat$Class_y,
                     MoreArgs = list(links = cl), SIMPLIFY = FALSE)
    link_n <- do.call("rbind.data.frame", args = link_n)

    return(link_n)
}

.countLinkPatterns <- function(x, y, links) {
    # a subordinate function of countAMRClassLinks
    n <- sum(links$Class_y == y & links$Class_x == x)

    return(data.frame(Class_y = y, Class_x = x, Links = n, stringsAsFactors = FALSE))
}
