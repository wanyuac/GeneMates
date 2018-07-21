# A module (private functions) of the phylix package for association analysis
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 17 March - 10 April 2017; the latest edition: 21 July 2018

##### Processing core-genome SNP data ###############
# Construct a matrix of major and minor alleles of every biallelic SNPs
.getAllelePairs <- function(bi) {
    snps <- colnames(bi)  # copy SNP names
    a <- matrix(nrow = 2, ncol = ncol(bi))
    colnames(a) <- snps
    rownames(a) <- c("major", "minor")

    for (snp in snps) {
        tab <- table(bi[, snp])
        alleles <- names(tab)
        count.1 <- as.integer(tab[alleles[1]])  # number of the first allele
        count.2 <- as.integer(tab[alleles[2]])  # number of the second allele

        if (count.1 > count.2) {  # when allele 1 is the major allele
            a["major", snp] <- alleles[1]
            a["minor", snp] <- alleles[2]
        } else {
            a["major", snp] <- alleles[2]
            a["minor", snp] <- alleles[1]
        }
    }

    return(a)
}

# Substitute every major allele with naught and every minor allele with one for every SNP ###############
# Arguments:
#   bi: a matrix of biallelic SNPs. SNP names define columns.
#   a: a matrix generated using .getAllelePairs for minor and major alleles per SNP.
.encodeAlleles <- function(bi, a) {
    m <- matrix(0, nrow = nrow(bi), ncol = ncol(bi))  # of the same dimension as the matrix bi
    snps <- colnames(bi)
    rownames(m) <- rownames(bi)
    colnames(m) <- snps
    for (snp in snps) {
        m[, snp] <- as.integer(bi[, snp] == a["minor", snp])  # Minor alleles are coded as one.
    }

    return(m)
}

###### Processing genotype data of bacterial genes ###############
.getAlleleNames <- function(pam) {
    # works for the function importGeneticPAM
    # extracts allele names from every column
    # pam: a matrix whose row names are strain names and columns are gene names
    a <- NULL
    if (is.matrix(pam)) {
        for (i in 1 : ncol(pam)) {
            col <- pam[, i]  # take a column of PAM
            a.new <- unique(col[col != "-"])
            if (length(a.new) > 0) {
                a <- c(a, a.new)
            }
        }
    } else {
        stop("Argument error: the input must be a matrix.")
    }

    return(a)
}

# Compress a presence/absence matrix into patterns of allelic presence/absence
# pam: a presence/absence matrix, where columns are alleles and rows are strains
# pat.id keeps the original number and order of alleles.
# In the pattern matrix B, patterns are arranged from 1 to the maximum of pattern IDs.
.getPatterns <- function(pam) {
    # compress the allelic matrix into a pattern matrix by retriving the first row index from pam that matches to each pattern
    patterns <- factor(apply(pam, 2, paste, collapse = ""))  # assign a pattern to each column, eg., c(1,0,1) => "101"
    profiles.unique <- match(levels(patterns), patterns)
    pat.m <- pam[, profiles.unique]  # a pattern matrix with columns ordered by levels from 1 to the maximum
    n.pat <- ncol(pat.m)  # number of patterns
    colnames(pat.m) <- paste0("pat_", 1 : n.pat)

    # store the mapping from allele names to pattern IDs
    alleles.mapping <- data.frame(allele = colnames(pam),
                                  pattern = as.integer(patterns),
                                  stringsAsFactors = FALSE)  # get the level ID that every allele maps to

    # identify patterns that refer to identically distributed (idd) alleles
    pat.sizes <- data.frame(pattern = 1 : n.pat, size = integer(n.pat))
    level.ids <- as.integer(patterns)  # 1, 2, ..., n.pat
    pat.sizes$size <- as.integer(sapply(pat.sizes$pattern, function(i) sum(level.ids == i)))  # the number of alleles each pattern has

    return(list(B = pat.m, alle.pat = alleles.mapping, pat.sizes = pat.sizes))  # remove redundant alleles
}

# Preparing allele pairs for LMMs ###############
# Calculates how many times does an explantory allele co-occur with the response allele
.countCoocurrence <- function(alleles.x, allele.y, pam) {
    y <- pam[, allele.y]
    pam.x <- pam[, alleles.x]
    n.x <- length(alleles.x)
    if (n.x > 1) {
        counts <- as.integer(apply(pam.x, 2, function(x) sum(as.logical(x) & as.logical(y))))
    } else if (n.x == 1) {  # pam.x becomes a vector of integers
        counts <- sum(as.logical(pam.x) & as.logical(y))
    } else {
        stop("Argument error: alleles.x must not be of zero length.")
    }

    return(counts)
}

# Pairs each response pattern with its possible explanatory patterns for assessing between-pattern associations
# x: an p-by-n matrix about m resistance alleles, where p equals the number of patterns (p <= m)
# y.pat: ID of the pattern that is used as the response variable
# ag: a data frame containing columns "allele", "gene" and "pattern"
# prt: set TRUE to actually print results into text files; set FALSE when debugging
# compress: set TRUE to compress genotype files using gzip
# skip: whether to avoid overwriting existing output files
.printBtwPatternGenFile <- function(x, y.pat, ag, pam, min.co = 2, prefix, output.dir,
                                    prt = TRUE, skip = TRUE) {
    alleles.y <- subset(ag, pattern == y.pat)  # response alleles under the current pattern ID
    alleles.x <- subset(ag, pattern != y.pat)  # explanatory alleles that are not identically distributed with the response allele.
    ax.df <- data.frame(y = character(0), x = character(0), x_pat = integer(0), stringsAsFactors = FALSE)  # initialisation

    # list possible explanatory alleles of each response allele
    for (i in 1 : nrow(alleles.y)) {  # for every response allele under this pattern, find out all candidate explanatory alleles
        ln <- alleles.y[i, ]  # pull out one line from the data frame
        ay <- ln$allele  # choose an allele under the current y pattern; ln$allele = alleles.y$allele[i]
        gy <- ln$gene  # find out the gene corresponding to this allele
        ax <- subset(alleles.x, gene != gy)  # valid explanatory alleles to test for = those in different patterns and genes & each co-occur with the response allele at least min.co times
        if (min.co > 0) {
            ax <- ax[.countCoocurrence(alleles.x = ax$allele, allele.y = ay, pam = pam) >= min.co, ]  # ax = NULL when the original ax = NULL.
        }
        ax.n <- nrow(ax)
        if (ax.n > 0) {  # If there are any potential explanatory alleles left, add them to the data frame.
            ax <- ax[, c("allele", "pattern")]  # drop the column of genes
            names(ax) <- c("x", "x_pat")
            ax.df <- rbind.data.frame(ax.df, cbind.data.frame(y = rep(ay, times = ax.n), ax, stringsAsFactors = FALSE),
                                      stringsAsFactors = FALSE)
        }  # Otherwise, omit these explanatory alleles.
    }

    # make a BIMBAM-format genotype file for GEMMA
    if (nrow(ax.df) > 0) {
        pat.ids <- unique(ax.df$x_pat)  # all patterns of X to be included under the current Y pattern
        pat.ids <- sort(pat.ids, decreasing = FALSE)  # patterns to be used as predictors for independent tests; sorted by levels
        n <- length(pat.ids)  # number of remaining patterns in X
        if (n == 1) {  # x[pat.ids,] becomes a named numeric vector when n = 1
            X.bimbam <- data.frame(id = paste0("pat_", pat.ids), presence = 1, absence = 0,
                                   matrix(x[pat.ids, ], nrow = 1), stringsAsFactors = FALSE)
            names(X.bimbam)[4 : ncol(X.bimbam)] <- colnames(x)
        } else {
            X.bimbam <- data.frame(id = paste0("pat_", pat.ids),
                                   presence = rep(1, times = n),
                                   absence = rep(0, times = n),
                                   x[pat.ids, ], stringsAsFactors = FALSE)
        }

        # make a pattern annotation file in the BIMBAM format
        # Herein I deliberately make the "SNP positions" (namely, the index column) to equal pattern IDs in order to simplify the processing of GEMMA outputs.
        X.annots <- data.frame(id = X.bimbam$id, index = pat.ids, chr = rep(24, times = n),
                               stringsAsFactors = FALSE)
        if (prt) {
            .writeData(x = X.bimbam,
                       output = paste0(output.dir, "/", prefix, "y__", y.pat, "__xs.txt"),
                       message = "[Genotype matrix X for inter-pattern associations]",
                       skip = skip)  # already under "output.dir/gene/"
            .writeData(x = X.annots,
                       output = paste0(output.dir, "/", prefix, "y__", y.pat, "__xs_annot.txt"),
                       message = "[Genotype annotations for inter-pattern associations]",
                       skip = skip)
        }
    }
    else {  # if nothing of X patterns is found at all for the current Y pattern
        print(paste0("Warning: no X patterns can be tested with the Y pattern of ", y.pat))
		ax.df <- NULL
    }

    return(ax.df)
}

.printIddPatternGenFile <- function(x, y.pat, ag, pam, min.co = 2, prefix, output.dir, prt = TRUE) {
    ag <- subset(ag, pattern == y.pat)  # for both X and Y patterns
    ax.df <- data.frame(y = character(0), x = character(0), x_pat = integer(0), stringsAsFactors = FALSE)  # initialisation

    # list possible explanatory alleles of each response allele
    n <- nrow(ag)  # obviously, n >= 2 for patterns of idd alleles. This is guaranteed.
    for (i in 1 : (n - 1)) {
        ay <- ag$allele[i]
        gy <- ag$gene[i]
        for (j in (i + 1) : n) {
            ax <- ag$allele[j]
            gx <- ag$gene[j]
            if (min.co > 0) {
                # Push a valid pair of alleles into this stack; When the PAM is produced by SRST2, gy != gx actually always true for idd. alleles.
                if ((gy != gx) & (.countCoocurrence(alleles.x = ax, allele.y = ay, pam = pam) >= min.co)) {
                    ax.df <- rbind(ax.df, data.frame(y = c(ay, ax), x = c(ax, ay), x_pat = c(y.pat, y.pat), stringsAsFactors = FALSE))  # symmetric pairwise tests
                }
            } else {
                if (gy != gx) {  # The only requirement here is that gene X != gene Y.
                    ax.df <- rbind(ax.df, data.frame(y = c(ay, ax), x = c(ax, ay), x_pat = c(y.pat, y.pat), stringsAsFactors = FALSE))
                }
            }
        }
    }

    # make a BIMBAM-format genotype file for GEMMA
    # In fact, we only need the log file to assess the population structure for idd. alleles.
    # Therefore, no genotype file will be created here if the pattern is already tested for between-pattern association;
    # otherwise, a new genotype file of a randomised (of the order) x will be saved to get the log file.
    # The .assoc.txt file produced for arbitary X's will not be used.
    # Randomisation of x (y) is taken because GEMMA cannot fit the LMM when x = y.
    if (nrow(ax.df) > 0) {  # if the Y allele has any explanatory X alleles (of the same pattern)
        genotype.file <- paste0(output.dir, "/", prefix, "y__", y.pat, "__xs.txt")  # already under "output.dir/gene/"
        if (!file.exists(genotype.file)) {  # If this pattern will be tested for between-pattern associations, then just read the output log file.
            x.ori <- x[y.pat, ]  # get the original Y pattern as X
            n <- length(x.ori)
            set.seed(100)
            x.random <- sample(x.ori, size = n)  # then randomise it so that estimates for LMM parameters will converge
            while (sum(x.random == x.ori) == n) {  # resample until x.random != x.ori
                x.random <- sample(x.ori, size = n)
            }
            X.bimbam <- data.frame(id = paste0("pat_", y.pat), presence = 1, absence = 0,
                                   matrix(x.random, nrow = 1), stringsAsFactors = FALSE)
            names(X.bimbam)[4 : ncol(X.bimbam)] <- colnames(x)  # copy isolate names into X.bimbam

            # make a pattern annotation file in the BIMBAM format
            # Herein I deliberately make the "SNP positions" (namely, the index column) to equal pattern IDs in order to simplify the processing of GEMMA outputs.
            X.annots <- data.frame(id = X.bimbam$id, index = y.pat, chr = 24, stringsAsFactors = FALSE)
            if (prt) {  # should not use .writeData here as the skip option violates the expected behaviour of this function
                write.table(X.bimbam, file = genotype.file,
                            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
                write.table(X.annots, file = paste0(output.dir, "/", prefix, "y__", y.pat, "__xs_annot.txt"),
                            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
            }
        }
    }
    else {  # if nothing of X patterns is found at all for the current Y pattern
        print(paste0("Warning: no X patterns are found under the Y pattern of ", y.pat))
		ax.df <- NULL
    }

    return(ax.df)  # nrow(ax.df >= 0)
}

# Iteratively make BIMBAM-format genotype and annotation files of a specific pattern for GEMMA
# Alleles of the same genes will not show co-occurrence, which is a limitation of SRST2 (Alleles of the same gene never co-occur in allele profiles).
# Identically distributed alleles (sharing the same pattern) will not subject to LMM-based analysis.
# .pairTestAlleles calls .printPatternGenFile to pair each response pattern with its possible explanatory patterns.
.pairTestAlleles <- function(X, gene.alleles, pam, min.co = 2, prefix, prt = FALSE,
                             output.dir, skip = TRUE) {
    # initialisation
    print(paste0("Making genotype and annotation files of valid explanatory alleles under ", output.dir, "."))
    n.pat <- nrow(X)  # number of all patterns, assuming rows are pattern values and column names are isolate names
    y.patterns <- vector(mode = "list", length = 2)
    names(y.patterns) <- c("dif", "idd")
    y.patterns[["dif"]] <- NULL
    y.patterns[["idd"]] <- NULL
    alleles.test <- vector(mode = "list", length = 2)
    names(alleles.test) <- c("dif", "idd")  # differently and identically distributed alleles
    alleles.test[["dif"]] <- NULL
    alleles.test[["idd"]] <- NULL  # create a null data frame
    if (!is.null(prefix)) {
        prefix <- paste0(prefix, "__")
    }

    # tests for associations between alleles that are not identically distributed
    # go through every pattern, selecting one as the response pattern, even if the pattern size > 1
    # This is because we use patterns as representatives of alleles and test for between-pattern associations,
    # which is irrelevant to the size of each pattern.
    for (i in 1 : n.pat) {
        new.df <- .printBtwPatternGenFile(x = X, y.pat = i, ag = gene.alleles, pam = pam, min.co = min.co,
                                          prefix = prefix, output.dir = output.dir, prt = prt,
                                          skip = skip)  # set the i-th pattern as y and the rest of patterns as x's
        if (!is.null(new.df)) {  # Sometimes there is no X patterns to be tested for the current Y pattern due to various reasons.
            new.df <- cbind.data.frame(y_pat = rep(i, times = nrow(new.df)), new.df, stringsAsFactors = FALSE)
            alleles.test[["dif"]] <- rbind.data.frame(alleles.test[["dif"]], new.df[, c("y", "x", "y_pat", "x_pat")],
                                                      stringsAsFactors = FALSE)
            y.patterns[["dif"]] <- append(y.patterns[["dif"]], i)
        }
    }

    # for idd. alleles
    ag.idd <- subset(gene.alleles, pat.size > 1)
    if (nrow(ag.idd) > 0) {  # when idd. alleles are present
        pat.idd <- unique(ag.idd$pattern)  # patterns having idd alleles
        for (i in pat.idd) {
            new.df <- .printIddPatternGenFile(x = X, y.pat = i, ag = gene.alleles, pam = pam, min.co = min.co,
                                              prefix = prefix, output.dir = output.dir, prt = prt)
            if (!is.null(new.df)) {  # Should always be true for idd alleles.
                new.df <- cbind.data.frame(y_pat = rep(i, times = nrow(new.df)), new.df, stringsAsFactors = FALSE)
                alleles.test[["idd"]] <- rbind.data.frame(alleles.test[["idd"]], new.df[, c("y", "x", "y_pat", "x_pat")],
                                                          stringsAsFactors = FALSE)
                y.patterns[["idd"]] <- append(y.patterns[["idd"]], i)
            }
        }
    } else {
        print("There is no identically distributed alleles.")  # Both alleles.test[["idd"]] and y.patterns[["idd"]] equal NULL.
    }

    return(list(dif = list(tests = alleles.test[["dif"]], y.pats = y.patterns[["dif"]]),
                idd = list(tests = alleles.test[["idd"]], y.pats = y.patterns[["idd"]])))
}

##### Run GEMMA iteratively for Y patterns ###############
# idd: whether alleles are identically distributed or not
# This function launch an independent GEMMA task for every Y pattern where corresponding output files are not present.
# Input files of GEMMA for each Y pattern have been generated in the previous function .pairTestAlleles.
# This function only needs to recursively compose and submit GEMMA command lines using each Y pattern.
# Here, assuming there are patterns to test.
.runGEMMA <- function(y.patterns, y.file, k.file, genotype.dir, prefix, n.cores = -1,
                      gemma, skip = TRUE) {
    n.pat <- length(y.patterns)  # Actual number of GEMMA runs <= length(y.patterns) (Some results may had been generated before).

    if (n.pat > 0) {  # when there are some patterns to be tested
        print(paste0("Running GEMMA for ", n.pat, " Y patterns."))

        # initialisation
        if (!is.null(prefix)) {
            prefix <- paste0(prefix, "__")
        }
        fn <- data.frame(y_pat = y.patterns,
                         lmm_file = sapply(y.patterns, function(i) paste0("output/", prefix, "y__", i, ".assoc.txt")),
                         log_file = sapply(y.patterns, function(i) paste0("output/", prefix, "y__", i, ".log.txt")),
                         stringsAsFactors = FALSE)  # names of output files from GEMMA

        # check whether the Y pattern has been tested for associations or not
        # This function does not automatically overwrite existant outputs.
        # Delete corresponding output files if you want to regenerate some results.
        # This behaviour may be changed in the future, depending on users' demand.
        run <- rep(TRUE, times = n.pat)
        for (i in 1 : n.pat) {
            if (file.exists(fn$lmm_file[i]) & file.exists(fn$log_file[i]) & skip) {
                run[i] <- FALSE  # skip this pattern for LMMs (notice "i" is not the pattern ID but the index)
            }
        }
        if (sum(run) > 0) {  # There are still some patterns to run.
            y.patterns <- y.patterns[run]
            print(paste("Actually to run", length(y.patterns), "Y patterns.", sep = " "))

            # determine the number of cores
            require(parallel)

            # run GEMMA for the remaining Y patterns
            # Since we are only interested in log files for within-pattern associations,
            # GEMMA only runs for idd. alleles whose patterns have not been used to test for between-pattern associations.
            n.cores <- .setCoreNum(n.cores = n.cores, cores.avai = detectCores())  # cores.avai >= 1
            cl <- makeCluster(n.cores)  # n.cores >= 1
            clusterExport(cl = cl,
                          varlist = list("gemma", "genotype.dir", "prefix", "y.file", "k.file"),
                          envir = environment())  # make variables accessible to different cores
            parLapply(cl, y.patterns, function(i) system(paste(gemma,
                                                               "-g", paste0(genotype.dir, "/", prefix, "y__", i, "__xs.txt"),
                                                               "-p", y.file, "-n", i,
                                                               "-a", paste0(genotype.dir, "/", prefix, "y__", i, "__xs_annot.txt"),
                                                               "-k", k.file, "-lmm 1 -notsnp -o", paste0(prefix, "y__", i),
                                                               sep = " "), wait = TRUE))  # wait until all jobs finish
            stopCluster(cl)
        } else {
            print("Skip this stage as all Y patterns have been run.")
        }
    } else {
        fn <- NULL
        stop("Argument error: the variable y.patterns should contain at least one element.")
    }

    return(fn)
}

##### Summarise results ###############
.readLine <- function(line) {  # a subordinate function of .readGEMMALogs
    a <- unlist(strsplit(line, " "))
    a <- as.numeric(a[length(a)])  # take the last word, which is the value of interest

    return(a)
}

# returns a single-line data frame
# Expect log files from GEMMA 0.96.
.readGEMMALogs <- function(log.file, i) {
    f <- scan(log.file, what = character(0), sep = "\n")
    logl_H0 <- .readLine(f[15])  # REMLE log-likelihood in the null model
    vg <- .readLine(f[19])  # the 17-th line "log-ML in the null (linear mixed) model"
    ve <- .readLine(f[20])  # take the last word as well
    l.remle <- round(vg / ve, digits = 9)  # keep the same precision level as vg and ve in GEMMA's outputs

    # return a data frame of five variables
    return(data.frame(y_pat = i, lambda0_remle = l.remle, logl_H0 = logl_H0, vg = vg, ve = ve))
}

# lmm.outputs equals outputs[["lmms.pat"]]
.readFittedLMMs <- function(lmm.outputs) {
    print("Reading results from GEMMA.")
    allele.types <- names(lmm.outputs)  # c("dif", "idd") or "dif" (when there are no idd. alleles)
    lmms.h1 <- vector(mode = "list", length = length(allele.types))
    names(lmms.h1) <- allele.types
    lmms.h0 <- lmms.h1

    # reading LMMs fitted under the alternative hypothesis H1
    for (t in allele.types) {
        lmms.h1[[t]] <- NULL
        lmms.h0[[t]] <- NULL
        n <- nrow(lmm.outputs[[t]])  # How many pairs of output files are there.
        if (n > 0) {  # No action is taken when lmms.outputs[["idd"]] = NULL (n always >0 for lmms.outputs[["dif"]]).
            for (i in 1 : n) {  # read every file
                r <- lmm.outputs[[t]][i, ]  # pull out a single row out of the data frame
                y.pat <- r$y_pat

                # for *.assoc.txt (associations under alternative hypothesis)
                if (t == "dif") {
                    h1 <- read.delim(r$lmm_file, stringsAsFactors = FALSE)  # read a tab-delimited file including its header line (by default)
                    h1 <- h1[, c("ps", "beta", "se", "l_remle", "p_wald")]  # discard useless columns; Actually, the ps column consists of the pattern ID of each X.
                    names(h1)[1] <- "x_pat"  # relace "ps" with "x_pat"
                    h1 <- data.frame(y_pat = rep(y.pat, times = nrow(h1)), h1)  # attach y_pat to the beginning of this data frame
                } else {
                    # We need to correct parameters under H1 for idd. alleles because LMMs cannot be fitted in these alleles.
                    # For idd. allele pairs, results under H1 are not useful. We only consider the results under H0.
                    # Basically, we do not test for any association between idd. alleles through LMMs.
                    # The purpose of constructing a data frame under H1 is to simplify the organisation of physical gene distances.
                    h1 <- data.frame(y_pat = y.pat, x_pat = y.pat, beta = 1, se = 0, l_remle = NA, p_wald = NA)
                }
                lmms.h1[[t]] <- rbind.data.frame(lmms.h1[[t]], h1, stringsAsFactors = FALSE)

                # for *.log.txt (null hypothesis)
                # Obviously, lmms.h0[["dif"]] and lmms.h0[["idd"]] share the same log file for any pair of patterns.
                h0 <- .readGEMMALogs(r$log_file, y.pat)  # h0 contains only a single row.
                lmms.h0[[t]] <- rbind.data.frame(lmms.h0[[t]], h0, stringsAsFactors = FALSE)
            }
        }
    }

    # Bonferroni correction of P values under H1 only for between-pattern associations
    # There are not P values in LMMs under H0.
    # The findPhysLink function guarantees that lmms.h1[[i]] is not empty.
    for (i in names(lmms.h1)) {
        if (i == "dif") {
            lmms.h1[[i]]$p_adj <- p.adjust(p = lmms.h1[[i]]$p_wald, method = "bonferroni")
        } else  {  # within-pattern associations: arbitrarily assign zeros to all P values when there are patterns of idd. alleles
            lmms.h1[[i]]$p_adj <- rep(NA, times = nrow(lmms.h1[[i]]))
        }
        lmms.h1[[i]] <- lmms.h1[[i]][, c("y_pat", "x_pat", "beta", "se", "l_remle", "p_wald", "p_adj")]  # rearrange columns
    }

    # return values
    if (length(lmms.h1) == 2) {
        z <- list(dif = list(h1 = lmms.h1[["dif"]], h0 = lmms.h0[["dif"]]),
                  idd = list(h1 = lmms.h1[["idd"]], h0 = lmms.h0[["idd"]]))
    } else {
        z <- list(dif = list(h1 = lmms.h1[["dif"]], h0 = lmms.h0[["dif"]]))
    }

    return(z)
}

# Recover allele-level associations from pattern-level associations ==========
# retrieve allele names for LMM results under the alternative hypothesis H1
# This is a subordinate function of the .restoreAlleleNames function.
# lmms: a data frame; tests: a list of two elements (tests and y.pats)
# h1: whether LMM are fitted under the alternative hypothesis or not.
# mapping: the data frame for allele-gene mappings is not necessary when h1 = TRUE.
.patternToAlleles <- function(lmms, tests, h1 = TRUE, mapping = NULL) {
    lmms.a <- NULL  # allelic results from LMMs

    if (h1) {
        print("Processing LMM results obtained under the alternative hypothesis.")
        y.pats <- tests[["y.pats"]]  # an integer vector of all Y patterns in the LMM results

        for (p in y.pats) {  # for each Y pattern, retrieve information about its corresponding X patterns
            lmms.y <- subset(lmms, y_pat == p)  # extract results from LMMs corresponding to a specific Y pattern
            tests.y <- subset(tests[["tests"]], y_pat == p)  # extract information of x, y and x_pat
            alleles.y <- unique(tests.y$y)  # Y alleles under the current Y pattern. There may be one or more alleles under the same Y pattern.

            # recover X alleles that are paired with each Y allele of the current Y pattern p
            for (a in alleles.y) {
                tests.xy <- subset(tests.y, y == a)  # the allele-pattern mapping table for X under the current Y pattern

                # Then retrieve results from lmms.y according to this allele-pattern mapping table.
                # This step recovers results from LMMs following the order of original alleles in the data frame "tests".
                df <- lmms.y[match(tests.xy$x_pat, lmms.y$x_pat), ]  # Note that nrow(df) > nrow(lmms.y) when the cardinality of at least one X pattern > 1.
                names(df)[1 : 2] <- c("y", "x")  # rename "y_pat" and "x_pat" to "y" and "x", respectively
                df$y_pat <- df$y  # make a copy of Y patterns
                df$x_pat <- df$x

                # substitute pattern IDs with allele names
                df$y <- rep(a, times = nrow(df))
                df$x <- tests.xy$x  # replace pattern IDs of X with allele names in an order of that in a.p$x
                lmms.a <- rbind.data.frame(lmms.a, df, stringsAsFactors = FALSE)
            }
        }
    } else {  # results under H0: only involve with Y alleles and patterns
        if (is.null(mapping)) {
            print("Argument error: a data frame of allele-gene mappings must be provided when h1 is FALSE.")
        } else {
            print("Processing LMM results obtained under the null hypothesis.")
            mapping <- subset(mapping, pattern %in% tests[["y.pats"]])
            lmms.a <- lmms[match(mapping$pattern, lmms$y_pat), ]  # rearrange rows of lmms according to mapping$pattern
            lmms.a <- cbind.data.frame(y = lmms.a, lmms.a, stringsAsFactors = FALSE)  # keep the y_pat column
            lmms.a$y <- mapping$allele
        }
    }

    return(lmms.a)
}

# lmms.pat: results from .readFittedLMMs; tests: alleles tested with every allele (from different genes) of each Y pattern
# The first two column names of lmm.a must be c("y", "x"). lmm.a eventually becomes the result data frame.
.restoreAlleleNames <- function(lmms.pat, tests, mapping) {
    print("Recovering allele names from pattern IDs.")

    # initialisation
    if ("idd" %in% names(lmms.pat)) {
        lmms <- list(dif = list(h1 = NULL, h0 = NULL),
                     idd = list(h1 = NULL, h0 = NULL))
    } else {
        lmms <- list(dif = list(h1 = NULL, h0 = NULL))
    }

    for (i in names(lmms.pat)) {  # c("dif", "idd") or "dif"
        lmms[[i]][["h1"]] <- .patternToAlleles(lmms = lmms.pat[[i]][["h1"]], tests = tests[[i]],
                                               h1 = TRUE, mapping = NULL)
        lmms[[i]][["h0"]] <- .patternToAlleles(lmms = lmms.pat[[i]][["h0"]], tests = tests[[i]],
                                               h1 = FALSE, mapping = mapping)
    }

    return(lmms)
}

# Retrieve allele counts and calculate the number of allelic co-occurrence events ==========
# Arguments
#   lmms: results of inter-pattern and intra-pattern allelic LMMs, follows the structure [["dif/idd"]][["h1/h0"]]
#   A: a presence/absence matrix of patterns
#   num: a data frame showing alleles and their numbers
.appendAlleleCounts <- function(lmms, A, num) {
    for (i in names(lmms)) {  # i = "dif" or "idd"
        for (j in names(lmms[[i]])) {
            d <- lmms[[i]][[j]]
            n <- nrow(d)
            d$n_y <- num$count[match(d$y, num$allele)]  # retrieve allele numbers from the data frame "num"
            if (j == "h1") {
                d$n_x <- num$count[match(d$x, num$allele)]
                d$n_xy <- integer(n)
                for (k in 1 : n) {  # calculate the number of co-occurrence events
                    ln <- d[k, c("y", "x")]  # extract a single line from the data frame
                    d$n_xy[k] <- sum(as.logical(A[, ln$y]) & as.logical(A[, ln$x]))
                }
                d <- d[, c("y", "x", "y_pat", "x_pat", "pair", "n_y", "n_x", "n_xy", "beta", "se", "l_remle",
                           "p_wald", "p_adj")]
            } else if (j == "h0") {
                d <- d[, c("y", "y_pat", "n_y", "lambda0_remle", "logl_H0", "vg", "ve")]
            } else {
                stop("Key error: the inner element name must be either \"h1\" or \"h0\".")
            }
            lmms[[i]][[j]] <- d
        }
    }

    return(lmms)
}

# Import a phylogenetic tree. This is a subordinate function of lmm.
# To use this function when tree = NULL, tree.proj cannot be NULL.
.importTree <- function(tree = NULL, outliers = NULL, tree.proj = NULL,
                        ref.rename = NULL) {
    require(ape)
    require(phytools)

    # IMPORT THE RAW TREE ###############
    tree.class <- class(tree)
    if (tree.class == "NULL") {
        if (!is.null(tree.proj)) {
            tree <- tree.proj  # use the projection tree when the external tree is not specified
        } else {
            stop("[Import tree] Error: either an external tree or a projection tree must be provided.")
        }
    } else if (tree.class == "character") {
        tree <- read.tree(tree)  # import the tree from a file
    }  # Otherwise, GeneMates uses the user's tree (an phylo object).

    # MODIFY THE TREE FOR AN APPROPRIATE STRUCTURE ###############
    # 1. Drop outlier tips
    if (!is.null(outliers)) {
        # Notice the drop.tip function ignores outliers that are not present
        # in the original tip labels. This is a desirable behaviour as the tree
        # may be estimated from a subset of the SNP table where some samples are
        # removed due to contamination or so. Therefore, we do not check if there
        # are any outliers absent on the tree.
        print("Dropping tips corresponding to outlier samples from the tree.")
        tree <- drop.tip(phy = tree, tip = outliers)  # It keeps the rooting condition.
    }

    # 2. Rename "Ref" when it is present and the ref.rename argument is given.
    if (!is.null(ref.rename)) {
        i <- which(tree$tip.label == "Ref")
        n <- length(i)
        if (n == 1) {
            print(paste("Replacing the tip label Ref with", ref.rename,
                        "on the tree.", sep = " "))
            tree$tip.label[i] <- ref.rename
        } else if (n > 1) {
            print("Warning: no tip label on the tree is changed because there are more than one \'Ref\'.")
        } else {
            print("No tip label on the tree is changed because the label \'Ref\' is not found.")
        }
    }

    # 3. The tree must be rooted for analysing the structural random effects as
    # this package assumes there are N - 1 internal nodes for a tree of N tips,
    # which only holds for a rooted tree. By contrast, an unrooted tree only has
    # N - 2 internal nodes.
    if (!is.rooted(tree)) {
        print("Midpoint rerooting the tree as it is unrooted.")
        tree <- midpoint.root(tree)
    }

    # 4. Check the minimum branch length and ensure it is positive, which is a
    # prerequisite for the ape::ace function. Otherwise, an error of "some branches
    # have length zero or negative" arises.
    len.nonpve <- tree$edge.length <= 0
    if (any(len.nonpve)) {  # There are N + (N - 1) - 1 = 2N - 2 edges for a rooted tree of N tips.
        print("Warning: some branches in the tree have zero or negative lengths.")
        print("Solution: arbitrarily adjust these branch lengths to 1/10 of the minimum positive length.")
        min.pve <- min(tree$edge.length[which(!len.nonpve)])
        tree$edge.length[which(len.nonpve)] <- min.pve / 10
    }

    return(tree)
}

# Determine sample distances
.determineSampleDists <- function(sample.dists = NULL, proj.dists = NULL,
                                  external.tree = FALSE, tree = NULL,
                                  outliers = NULL) {
    require(ape)

    if (is.null(sample.dists)) {
        if (external.tree & (!is.null(tree))) {
            # when a tree is provided and this tree is prepared by the user
            sample.dists <- cophenetic.phylo(tree)  # tr may be filtered against outlier samples
        } else {
            # when the tree is either a projection tree or absent
            sample.dists <- proj.dists
        }
    }  # Else, do nothing.

    # Exclude rows and columns from the distance matrix against outliers
    if (!is.null(outliers)) {
        samples <- rownames(sample.dists)
        samples <- samples[!(samples %in% outliers)]
        sample.dists <- sample.dists[samples, samples]
    }  # Otherwise, GeneMates uses the user's distance matrix, which can be, but not necessarily be a distance from the phylogenetic tree.

    return(sample.dists)
}
