#' @title Importing a table of physical distances as a data frame
#' @description This function reads a tabulated Bandage output into R for the
#' evaluation of distance consistency. The distance table can be created using
#' a helper pipeline physDist (github.com/wanyuac/physDist).
#' @param dists The distance table to be imported.
#' @param delim A character specifying the delimiter in the distance table. Default:
#' a tab character.
#' @param ingroup A vector of strain names to be included in the result.
#' @param outgroup A vector of strains whose physical distances will be excluded
#' from the result.
#' @return A data frame of physical distances.
#' @author Yu Wan, \email{wanyuac@@126.com}
#' @export
#  Copyright 2019 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 23 Janury - 27 March 2017; the lastest edition: 2 June 2019

# Read distances between alleles ###############
importPhysicalDists <- function(dists = NULL, delim = "\t", ingroup = NULL, outgroup = NULL) {
    # This function filters the distances for strain names, but it does not perform
    # further filter for the maximum distance or node number, which is done in
    # the function summariseDist.
    ds.class <- class(dists)
    if (ds.class == "character") {  # a path to the distance file
        print(paste0(Sys.time(), ": reading physical distances from ", dists))
        ds <- read.delim(dists, sep = delim, stringsAsFactors = FALSE)
        ds <- ds[, c("query1", "query2", "sample", "distance", "node_number",
                     "source", "orientation", "distance_path")]
    } else if (ds.class == "data.frame") {
        print("Skip reading distances as they have been imported.")
    } else {
        stop("Input error: the dists argument must be a file path or a data frame.")
    }

    # Filter distances to keep those in the SNP matrix
    # This step is pretty important. Otherwise, the function .estimateIBD will
    # return an error of "length of phenotypic and of phylogenetic data do not
    # match" because extra samples get selected for distance analysis.
    # A user may specify ingroup and/or outgroup strains for this filter.
    if (!is.null(ingroup)) {
        print("Filtering allelic physical distances for ingroup strains.")
        ds <- subset(ds, sample %in% ingroup)
    }
    if (!is.null(outgroup)) {
        print("Filtering allelic physical distances to remove outgroup strains.")
        ds <- subset(ds, !(sample %in% outgroup))
    }

    return(ds)
}
