# Functions for environment configurations.
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 11 May 2017; the latest edition: 3 Apr 2018

##### File and directory management ###############
# Check if every file on a list exists or not
# inputs: a vector of characters (absolute file paths and names)
# file.n: expected number of inputs
.checkFiles <- function(inputs = NULL, file.n = length(inputs)) {
    print("Checking existance of input files.")
    pass <- TRUE
    if (length(inputs) == file.n) {  # All input files must be available as expected.
        for (f in inputs) {
            if (!file.exists(f)) {
                print(paste("Path error: the file", f, "does not exist.", sep = " "))
                pass <- FALSE
            }
        }
    } else {
        print(paste("Error: not all of", file.n, "input files are provided!"), sep = " ")
        pass <- FALSE
    }

    return(pass)
}

.checkAssocOut <- function(assoc.out = NULL) {
    # A subordinate function of findPhysLink, which check if the argument for
    # assoc.out is valid.
    if (is.list(assoc.out)) {
        # check if assoc.out contains all obligate elements
        if (any(!(c("outputs", "alleles", "genes", "mapping", "tests", "lmms.pat",
                    "lmms", "struc") %in% names(assoc.out)))) {
            pass <- FALSE  # miss one or more elements
        } else {
            pass <- TRUE
        }
    } else {
        pass <- FALSE
    }

    return(pass)
}

# a generic function for writting the data frame or matrix x into storage
.writeData <- function(x, output, sep = "\t", message = NULL, skip = TRUE,
                       row.names = FALSE, col.names = FALSE, quote = FALSE) {
    if (file.exists(output) & skip) {
        print(paste0(message, " Skip writing the output file ", output, "."))
    } else {
        print(paste0(message, " Writing the file ", output, "."))
        write.table(x, file = output, sep = sep,
                    row.names = row.names, col.names = col.names, quote = quote)
    }
}

##### Initialise stages ###############
.initialiseStageRecords <- function(func = "lmm", save.stages = TRUE,
                                    prefix = NULL) {
    # func: which function does this function produces stage filenames for.
    # func = "findPhysLink" or "lmm"
    if (save.stages) {
        print("Stages will be saved under the /temp directory.")
        if (!file.exists("temp")) {  # check the presence of \temp
            dir.create("temp")
        }
        if (func == "lmm") {
            temp.files <- c("snps" = paste0("temp/", prefix, "__snps.rds"),
                            "genes" = paste0("temp/", prefix, "__genes.rds"),
                            "alleles" = paste0("temp/", prefix, "__alleles.rds"),
                            "gene.alleles" = paste0("temp/", prefix, "__gene_alleles.rds"),
                            "allele.pairs" = paste0("temp/", prefix, "__allele_pairs.rds"),
                            "lmms.pat.dif" = paste0("temp/", prefix, "__lmms_pat_dif.rds"),
                            "lmms.pat.idd" = paste0("temp/", prefix, "__lmms_pat_idd.rds"),
                            "lmms" = paste0("temp/", prefix, "__lmms.rds"),
                            "struc" = paste0("temp/", prefix, "__population_structure.rds"))
        } else {
            temp.files <- c("ds" = paste0("temp/", prefix, "__ds.rds"),
                            "ds.summary" = paste0("temp/", prefix, "__ds_summary.rds"),
                            "assoc" = paste0("temp/", prefix, "__assoc.rds"))
        }
    } else {
        print("Warning: stage outputs will not be saved.")
        temp.files <- NULL
    }

    return(temp.files)
}

# Read an RDS file into a variable
.recoverHistory <- function(f) {
    print(paste0("Loading previous record ", f))
    record <- readRDS(f)  # import the list snps

    return(record)
}

# Generate the name of an output file
.makeFileName <- function(dir.path, prefix, basename, delim = "__") {
    return(paste(dir.path, paste0(prefix, delim, basename), sep = "/"))
}

# Create a new directory if "d" is not found
# d: either an absolute path or a relative path of the target directory
.checkDir <- function(dirs) {
    for (d in dirs) {
        if (!dir.exists(d)) {
            dir.create(d)
        }
    }
}

# Determine the number of cores for parallel computing
# n.cores: a user-specified number of cores (cf. the parameter n.cores of the findPhysLink function)
# n.avai: number of cores detected using the function detectCores() of the parallel function
.setCoreNum <- function(n.cores, cores.avai) {
    if (n.cores < 1) {  # expect n.cores = -1 or 0 in this case
        if (n.cores < -1) {  # In case a user gives a wrong argument.
            n.cores <- -1
        }
        n.cores <- cores.avai + n.cores  # n.cores = -1: leave one core free for system management (recommended)
        if (n.cores == 0) {  # when there is only a single core available
            n.cores <- 1
        }
        print(paste(cores.avai, "cores have been detected and", n.cores, "are used to launch parallel tasks.", sep = " "))
    } else {  # when n.cores >= 1
        if (n.cores > cores.avai) {  # In case a user specifies an excessive number of cores
            print(paste("Warning:", n.cores, "cores are specified, but only", cores.avai, "cores are available.", sep = " "))
            n.cores <- cores.avai - 1  # Do not use all cores for safty.
            if (n.cores == 0) {
                n.cores <- 1  # when cores.avai = 1
            }
            print(paste("Adjusted to use", n.cores, "to launch parallel tasks.", sep = " "))
        } else {
            print(paste("Use", n.cores, "cores as specified to launch parallel tasks.", sep = " "))
        }
    }

    return(n.cores)
}
