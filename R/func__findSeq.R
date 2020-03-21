#' @title Search for a query sequence against a list of assemblies
#'
#' @description This function determines presence of a query sequence in a list of assembly graphs or FASTA files (readable to Bandage). It aims to
#' answer the question that how frequent a query is found in a collection of assemblies. This function does not work on Window OS
#' unless the Linux commandline "cut" is enabled.
#'
#' @param query Path to a FASTA file, which may contain multiple query sequences.
#' @param assemblies A data frame, a character matrix or a CSV file whose first two columns provide strain names and paths to assembly files. These two columns may be
#' named Strain and Assembly for instance. This argument can also be a path to a CSV file (with a header line for column names) for this data frame. For
#' Bandage, a valid assembly file can be either a SPAdes FASTG file or a FASTA file. This function searches the query in every assembly file. Users may use
#' a spreadsheet to create a CSV file for this data frame and import it into R.
#' @param bandage.path Path to Bandage, without any backslash or forward slash terminating this parameter.
#' @param blast.params Parameters passed directly to BLAST through the option "--blastp" of Bandage. Run "bandage --helpall" for details. Default: megablast.
#' @param bandage.params Parameters passed directly to Bandage. Run "bandage --helpall" as well to see all valid parameters. These parameters controls
#' how Bandage identifies a query.
#' @param n.cores Number of computational cores that will be used in parallel for this function. It follows the same convention defined in the function
#' findPhysLink. For simplicity, set it to zero to automatically detect and use all available cores; set it to -1 to leave one core out (recommended unless
#' this function is executed through an SLURM job system).
#' @param del.temp A logical parameter determing whether to keep all temporal files under the current working directory. Default: removing all of these files.
#'
#' @return A single data frame of identified query paths, one (the top hit) for each assembly. NA values are present if no query path is found at all in an
#' assembly.
#'
#' @examples
#' paths <- findSeq(query = "integrons.fna", assemblies = a, bandage.path = "apps/Bandage",
#' bandage.params = "--ifilter 95 --evfilter 1e-3 --pathnodes 6 --minhitcov 0.98 --minpatlen 0.98 --maxpatlen 1.02",
#' n.cores = 4, del.temp = FALSE)
#'
#' @author Yu Wan (\email{wanyuac@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  The first and the lastest edition: 22 August 2017

findSeq <- function(query = NULL, assemblies = NULL, bandage.path = "./bandage",
                    blast.params = "-task megablast",
                    bandage.params = "--ifilter 95 --evfilter 1e-3 --pathnodes 6 --minhitcov 0.98 --minpatlen 0.98 --maxpatlen 1.02",
                    n.cores = -1, del.temp = TRUE) {

    ##### Check if all paths are reachable ##########
    if (is.character(assemblies)) {  # when a path is supplied instead of a data frame or a matrix
        assemblies <- read.csv(file = assemblies, stringsAsFactors = FALSE)  # An error arises when this file is not accessible.
    }

    if (.checkFiles(inputs = c(query, assemblies[, 2], bandage.path))) {
        print("Pass: all input files are accessible.")
    } else {
        stop("Input error: some files are inaccessible.")
    }

    ##### Search for the query in every assembly ##########
    require(parallel)
    require(data.table)

    # convert the data frame assemblies into a character matrix
    if (class(assemblies) == "data.frame") {
        assemblies <- as.matrix(assemblies)
        rownames(assemblies) <- NULL
    }

    # add quotes to blast.params, otherwise, Bandage reports an error: "Invalid option: megablast" (or blastn)
    blast.params <- paste0("\"", blast.params, "\"")

    print(paste0(Sys.time(), ": searching for the query sequence in ", nrow(assemblies), " assemblies."))
    cl <- makeCluster(.setCoreNum(n.cores = n.cores, cores.avai = detectCores()))
    clusterExport(cl = cl, varlist = list("query", "bandage.path", "blast.params", "bandage.params",
                                          "del.temp"), envir = environment())
    paths <- parApply(cl, assemblies, 1, .searchForQuery, query, bandage.path, blast.params, bandage.params, del.temp)
    stopCluster(cl)
    print(paste0(Sys.time(), ": finished all searches successfully."))

    return(rbindlist(paths))
}

# This is a subordinate function of findSeq. It searches paths for every query in a FASTA file in an
# assembly file.
# r: a row in the data frame assemblies
.searchForQuery <- function(r, query, bandage.path, blast.params = "'-task megablast'", bandage.params, del.temp = TRUE) {
    strain <- r[[1]]  # Strain names are stored in the first column of the input data frame assemblies.
    assembly <- r[[2]]  # path to the assembly file (in FASTG, GFA or FASTA format)
    cmd <- paste(bandage.path, "querypaths", assembly, query, strain, "--blastp", blast.params,
                 bandage.params, sep = " ")
    system(cmd, wait = TRUE)  # run this command

    # chop off the last column (sequence) from the output file of Bandage
    out1 <- paste0(strain, ".tsv")  # output file of Bandage
    out2 <- paste0(strain, "__trimmed.tsv")  # the processed output file
    if (file.exists(out1)) {
        system(paste("cut -f1-11", out1, ">", out2, sep = " "), wait = TRUE)  # not supported on Windows OS unless the cut command is installed
    } else {
        stop(paste0("Executive error: output file ", out1, " is not produced successfully."))
    }

    # process Bandgae output
    paths <- read.delim(file = out2, stringsAsFactors = FALSE)  # 11 columns
    n <- nrow(paths)
    if (n == 0) {  # An empty data frame indicates absence of any query paths.
        paths <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(paths)))
        n <- 1
    }
    paths <- cbind.data.frame(Strain = rep(strain, times = n), paths, stringsAsFactors = FALSE)
    names(paths) <- c("Strain", "Query", "Path", "Length", "Query_covered_by_path", "Query_covered_by_hits",
                      "Mean_hit_identity", "Total_hit_mismatches", "Total_hit_gap_opens", "Relative_length",
                      "Length_discrepancy", "E_value_product")  # unify column names

    # clearance
    if (del.temp) {
        system(paste("rm", out1, out2, sep = " "), wait = TRUE)
    }

    return(paths)
}
