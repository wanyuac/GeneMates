#' @title Testing for associations between occurrence of bacterial genes/alleles with linear mixed model
#'
#' @description This is one of the two main functions of this package (the other main function is findPhysLink). It estimates parameters
#' of univariate linear mixed models (LMMs) to test for associations between occurrence of bacterial genes or alleles when control for
#' bacterial population structure. This function does not utilise allelic physical distances at all.
#'
#' Dependency: packages data.table and parallel
#'
#' @param snps Core-genome SNPs used for estimating the relatedness matrix.
#      Valid values: a complete path to a SNP table, or an un-centred, encoded, biallelic SNP matrix G
#' @param snps.delim (optional) Delimiters of fields in the SNP table (pam.delim and dist.delim are defined similary)
#' @param pos.col (optional) An integer (column index) or a string (column name) specifying which column contains SNP positions
#' @param ingroup (optional) A vector of characters for names of isolates to be analysed. Isolates may be sorted, such as according to the phylogeny. The function includes all isolates by default.
#' @param outliers (optional) A vector of characters for isolate/strain names to be excluded from snps, pam
#' @param ref (optional) name of the reference genome. The column name "Ref" in the SNP matrix will be replaced with this argument.
#' @param min.mac (optional) An integer specifying the minimal number of times required for the minor allele of every biallelic SNP to occur across all isolates. SNPs failed this criterion will be removed from this analysis.
#' @param genes.excl (optional) Genes to be excluded from PAMs. For example, genes.excl = c("AmpH_Bla", "OqxBgb_Flq", "OqxA_Flq", "SHV.OKP.LEN_Bla").
#' @param genetic.pam A presence/absence matrix of genes. It may be a compiled table from SRST2.
#' @param allelic.pam A presence/absence matrix of alleles. The matrix may be a compiled table of SRST2 results.
#' @param genetic.pam.delim (optional) A delimiter character in the genetic PAM. Default: tab.
#' @param allelic.pam.delim (optional) A delimiter character in the allelic PAM. Default: tab.
#' @param min.count (option) The minimum count of alleles/genes in the current data set to be included
#' for analysis.
#' @param min.co (optional) The minimum number of allelic co-occurrence events. Set it to zero to specify the tests even though
#' the corresponding alleles never co-occur. However, this will increase the number of tests tremendously.
#' @param mapping (optional) A data frame mapping alleles to genes and patterns, etc, which equals the "mapping" element within findPhysLink's output list.
#'      This argument is only used when a user reruns a previous analysis.
#' @param tree (optional) A path to a tree file or a phylo object for a tree of all samples.
#' The format of the tree file must be compartile to the read.tree function in the ape package.
#' @param sample.dists (optional) A numeric matrix of distances between samples. The distances can be Euclidean distances between projections, phylogenetic
#' tip distances or SNP distances (the number of SNPs between any two samples). A matrix of Euclidean distances will be computed for projections
#' of samples if this option is left NULL.
#' @param output.dir (optional) Path of the output directory.
#       A relative path "output" under the current working directory is recommended as GEMMA always create a directory named output to store its outputs.
#       Otherwise, you will end up with two output directories: one for yours, and the other for GEMMA.
#' @param prefix (optional) For names of all output files
#' @param gemma.path Path to GEMMA. No forward slash should be attached at the end of the path.
#' @param n.cores Number of cores used to run GEMMA in parallel where possible.
#'      -1: automatically detect the number of available cores N, but use N - 1 cores (recommended)
#'      0: automatically detect the number of available cores and use all of them. Be careful when the current R session is not running through SLURM.
#'      >= 1: use the number of cores as specified. n.cores is reset to the maximal number of available cores N when n.cores > N.
#' @param save.stages (optional) Whether to turn on stage control or not. Recommend to turn it on when you are not sure whether the pipeline will be finished smoothly.
#' @param skip (optional) Whether to avoid overwriting existing output files.
#'
#' @author Yu Wan, \email{wanyuac@@gmail.com}
#' @export
#
#  Copyright 2017-2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 15 Dec 2017; the latest edition: 20 June 2018

# LMM-based association tests
lmm <- function(snps = NULL, snps.delim = ",", pos.col = "Pos", min.mac = 1,
                genetic.pam = NULL, genetic.pam.delim = "\t", genes.excl = NULL,
                allelic.pam = NULL, allelic.pam.delim = "\t", mapping = NULL,
                min.count = 2, min.co = 2,
                ingroup = NULL, outliers = NULL, ref = NULL,
                tree = NULL, sample.dists = NULL,
                output.dir = "output", prefix = NULL,
                gemma.path = "gemma", n.cores = -1,
                save.stages = TRUE, skip = TRUE) {

    # 1. Process core-genome SNPs ###############
    print(paste0(Sys.time(), ": initialising the lmm function."))

    # check whether these input files exist
    inputs <- list(snps, allelic.pam, genetic.pam, sample.dists, gemma.path)  # load file paths into a list
    input.files <- NULL
    for (item in inputs) {
        if (class(item) == "character") {
            input.files <- append(input.files, item)
        }  # Non-character inputs (such as a list or the NULL value) will not be checked for their paths.
    }
    if (length(input.files) > 0) {
        if (.checkFiles(inputs = input.files)) {
            print("Pass: all input files are accessible.")
        } else {
            stop("Input error: some file names are missing or some files are inaccessible.")
        }
    }

    # set up output files
    output.dirs <- c(main = output.dir,
                     snp = paste(output.dir, "snp", sep = "/"),
                     gene = paste(output.dir, "gene", sep = "/"))
    outputs <- list(snps = .makeFileName(output.dirs[["snp"]], prefix, "G.txt"),
                    snp.annots = .makeFileName(output.dirs[["snp"]], prefix, "snps_annot.txt"),
                    pam = .makeFileName(output.dirs[["gene"]], prefix, "pam.txt"),
                    Y = .makeFileName(output.dirs[["gene"]], prefix, "Y.txt"),
                    K = .makeFileName("output", paste0(prefix, "_K"), "cXX.txt", delim = "."),
                    U = .makeFileName("output", paste0(prefix, "_SVD"), "eigenU.txt", delim = "."),
                    D = .makeFileName("output", paste0(prefix, "_SVD"), "eigenD.txt", delim = "."),
                    lmms.pat = NULL)
    .checkDir(output.dirs)  # set up the main output directory

    # check whether any stages were completed or not
    stage.outputs <- .initialiseStageRecords("lmm", save.stages, prefix)  # returns a named vector of characters for paths to stage records

    # extract non-isolate-specific biallelic core-genome SNPs
    stage.record <- stage.outputs[["snps"]]
    if (save.stages & file.exists(stage.record)) {
        snps <- .recoverHistory(stage.record)
    } else {
        # elements of snps: G, S, G.bimbam, annots, snp.alleles, mac, core, bi and var
        snps <- importCoreGenomeSNPs(snps = snps, snps.delim = snps.delim,
                                     pos.col = pos.col, replace.ref = ref,
                                     min.mac = min.mac, ingroup = ingroup,
                                     outliers = outliers, G.file = outputs[["snps"]],
                                     annots.file = outputs[["snp.annots"]],
                                     skip = skip)
        if (save.stages) {
            saveRDS(snps, file = stage.record)
        }
    }

    # 2. Import allelic and genetic presence/absence matrices of bacterial genes ###############
    stage.record <- stage.outputs[["genes"]]
    if (save.stages & file.exists(stage.record)) {
        genes <- .recoverHistory(stage.record)
    } else {
        genes <- importGeneticPAM(pam = genetic.pam, pam.delim = genetic.pam.delim,
                                  outliers = outliers, min.count = min.count,
                                  genes.rm = genes.excl,
                                  sample.order = rownames(snps[["S"]]))
        if (save.stages) {
            saveRDS(genes, file = stage.record)
        }
    }

    stage.record <- stage.outputs[["alleles"]]
    if (save.stages & file.exists(stage.record)) {
        alleles <- .recoverHistory(stage.record)
    } else {
        # elements of alleles: X, Y, A, B, alle.pat, pat.sizes and alleles.idd
        alleles <- importAllelicPAM(pam = allelic.pam, pam.delim = allelic.pam.delim,
                                    outliers = outliers, min.count = min.count,
                                    alleles.inc = genes[["alleles.inc"]],
                                    output.y = outputs[["Y"]],
                                    sample.order = rownames(snps[["S"]]),
                                    skip = skip)
        if (save.stages) {
            saveRDS(alleles, file = stage.record)
        }
    }

    # 3. Make a data frame for three variables: genes, alleles and allele frequencies ###############
    stage.record <- stage.outputs[["gene.alleles"]]
    if (save.stages & file.exists(stage.record)) {
        gene.alleles <- .recoverHistory(stage.record)
    } else {
        if (class(mapping) == "data.frame") {
            if (sum(names(mapping) %in% c("allele", "gene", "pattern", "count", "freq", "pat.size")) == 6) {
                print("Skip computing the mapping data frame as it has been provided.")
                gene.alleles <- mapping
            } else {
                print("Regenerating the mapping data frame as it is incomplete.")
                # match every allele to a gene and count its frequency
                gene.alleles <- countAlleles(pam.genes = genes[["pam"]],
                                             pam.alleles = alleles[["A"]],
                                             patterns = alleles[["alle.pat"]],
                                             pat.sizes = alleles[["pat.sizes"]])
            }
        } else {  # for each fresh run
            gene.alleles <- countAlleles(pam.genes = genes[["pam"]],
                                         pam.alleles = alleles[["A"]],
                                         patterns = alleles[["alle.pat"]],
                                         pat.sizes = alleles[["pat.sizes"]])
        }

        if (save.stages) {
            saveRDS(gene.alleles, file = stage.record)
        }
    }

    # 4. Estimate the relatedness matrix K from core-genome SNPs ###############
    stage.record <- outputs[["K"]]  # equals the file name of the relatedness matrix
    if (file.exists(stage.record) & (save.stages | skip)) {
        print("Skip generating the relatedness matrix as it has been produced.")
    } else {
        cmd <- paste(gemma.path, "-g", outputs[["snps"]], "-p", outputs[["Y"]], "-n 1 -a", outputs[["snp.annots"]],
                     "-gk 1 -maf 0 -o", paste0(prefix, "_K"), sep = " ")  # a centred but not standardised relatedness matrix; "-n i" is arbitary and the value of i does not make a difference.
        print(paste0("Estimating the relatedness matrix through ", cmd))
        system(cmd, wait = TRUE)  # wait until the job finishes

        # check if the relatedness matrix was generated successfully or not
        if (!file.exists(stage.record)) {
            stop("Command line error: the relatedness matrix is not produced successfully!")
        }
    }

    # 5. Iteratively make genotype and annotation files for LMMs ###############
    # From now on, all data processing are divided by whether alleles are identically distributed (idd) or not.
    stage.record <- stage.outputs[["allele.pairs"]]
    if (save.stages & file.exists(stage.record)) {
        allele.pairs <- .recoverHistory(stage.record)
    } else {
        # When there are not any idd alleles:
        # allele.pairs = list(dif = list(tests, y.pats), idd = list(tests = NULL, y.pats = NULL))
        # The filter of allele pairs for co-occurrence counts applies here. Set min.co = 0 to turn it off.
        allele.pairs <- .pairTestAlleles(X = alleles[["X"]], gene.alleles = gene.alleles,
                                         pam = alleles[["A"]], min.co = min.co,
                                         prefix = prefix, prt = TRUE,
                                         output.dir = paste0(output.dir, "/gene"))  # elements: dif, idd
        if (save.stages) {
            saveRDS(allele.pairs, file = stage.record)
        }
    }
    has.idd <- !is.null(allele.pairs[["idd"]][["y.pats"]])  # a logical indicator for whether there are idd. alleles

    # 6. Run GEMMA iteratively to fit LMMs in accordance with allele.pairs ###############
    if (has.idd) {
        outputs[["lmms.pat"]] <- list(dif = NULL, idd = NULL)
    } else {
        outputs[["lmms.pat"]] <- list(dif = NULL)
    }

    # association between patterns
    stage.record <- stage.outputs[["lmms.pat.dif"]]
    if (save.stages & file.exists(stage.record)) {
        outputs[["lmms.pat"]][["dif"]] <- .recoverHistory(stage.record)
    } else {
        outputs[["lmms.pat"]][["dif"]] <- .runGEMMA(y.patterns = allele.pairs[["dif"]][["y.pats"]],
                                                    y.file = outputs[["Y"]],
                                                    k.file = outputs[["K"]],
                                                    genotype.dir = paste0(output.dir, "/gene"),
                                                    prefix = prefix, n.cores = n.cores,
                                                    gemma = gemma.path, skip = skip)  # a data frame with variables y.pat, lmm.file and log.file
        if (save.stages) {
            saveRDS(outputs[["lmms.pat"]][["dif"]], file = stage.record)
        }
    }

    # association between identically distributed alleles (self-tests of patterns)
    if (has.idd) {
        stage.record <- stage.outputs[["lmms.pat.idd"]]
        if (save.stages & file.exists(stage.record)) {
            outputs[["lmms.pat"]][["idd"]] <- .recoverHistory(stage.record)
        } else {
            # No GEMMA sessions will be launched by the following command line but only returns a data frame of file names.
            # The returned value may be NULL when there are not idd. alleles.
            outputs[["lmms.pat"]][["idd"]] <- .runGEMMA(y.patterns = allele.pairs[["idd"]][["y.pats"]],
                                                        y.file = outputs[["Y"]],
                                                        k.file = outputs[["K"]],
                                                        genotype.dir = paste0(output.dir, "/gene"),
                                                        prefix = prefix, n.cores = n.cores,
                                                        gemma = gemma.path, skip = skip)
            if (save.stages) {
                saveRDS(outputs[["lmms.pat"]][["idd"]], file = stage.record)
            }
        }
    }

    # 7. Import results of GEMMA ###############
    # This stage is not saved as it is not time consuming.
    # import fitted pattern-based LMMs under both hypotheses (H1 and H0) for non-idd and idd. alleles
    # For associations between idd. alleles, only log files are actually useful.
    # Structure of the output list:
    # list(dif = list(h1, h0), idd = list(h1, h0)) when there are idd. alleles;
    # otherwise, list(dif = list(h1, h0)).
    lmms.pat <- .readFittedLMMs(outputs[["lmms.pat"]])

    # recover allele-level associations from pattern-level associations
    stage.record <- stage.outputs[["lmms"]]
    if (save.stages & file.exists(stage.record)) {
        lmms <- .recoverHistory(stage.record)
    } else {
        # return a list with elements "dif" and "idd" or  with "dif" only
        lmms <- .restoreAlleleNames(lmms.pat = lmms.pat, tests = allele.pairs, mapping = gene.alleles)

        # assign an ID to each pair of alleles
        lmms[["dif"]][["h1"]] <- assignPairID(lmms = lmms[["dif"]][["h1"]], from = 1)
        if (has.idd) {
            lmms[["idd"]][["h1"]] <- assignPairID(lmms = lmms[["idd"]][["h1"]], from = nrow(lmms[["dif"]][["h1"]]) + 1)
        }

        # add allele counts and numbers of co-occurrence events
        lmms <- .appendAlleleCounts(lmms = lmms, A = alleles[["A"]], num = gene.alleles)

        if (save.stages) {
            saveRDS(lmms, file = stage.record)
        }
    }

    # 8. Resolving bacterial population structure ###############
    stage.record <- stage.outputs[["struc"]]  # temp/[prefix]__population_structure.rds
    if (file.exists(stage.record) & file.exists(outputs[["U"]]) & file.exists(outputs[["D"]]) & (save.stages | skip)) {
        print("Skip projecting samples and calculating sample distances as relevant results have been produced.")
        struc <- .recoverHistory(stage.record)
    } else {
        # projection of samples with singular-value decomposition of the biallelic cgSNP matrix
        # append "_SVD" to GEMMA's outputs at this stage to avoid overwriting the log file produced when generating the relatedness matrix.
        C <- projectSamples(K = outputs[["K"]], G = outputs[["snps"]], Y = outputs[["Y"]],
                            L = ncol(snps[["G"]]), samples = rownames(snps[["G"]]),
                            prefix = paste0(prefix, "_SVD"), get.dists = TRUE,
                            dist.method = "euclidean", get.tree = TRUE,
                            gemma.path = gemma.path)  # returns a list of elements C, sv, ev, U, K, d and tr

        # import a tree and convert its topology into a PAM of samples in each clade
        # The returned tree equals the projection tree when the argument tree = NULL.
        external.tree <- !is.null(tree)  # if the tree is a user tree
        tree <- .importTree(tree = tree, outliers = outliers, tree.proj = C[["tr"]])
        clades <- tree2Clades(tr = tree, sample.order = rownames(C[["C"]]))

        # find out the minimal inclusive clade for each pair of co-occurring alleles
        # The distribution of statistics, such as f_xy (co-occurrence frequency) and d_max (maximal sample distance),
        # help users to decide appropriate thresholds for inferring identity-by-descent.
        sample.dists <- .determineSampleDists(sample.dists = sample.dists,
                                              proj.dists = C[["d"]],
                                              external.tree = external.tree,
                                              tree = tree, outliers = outliers)
        mc.xy <- findMinIncClade(lmms = lmms, allele.pam = alleles[["A"]],
                                 clade.pam = clades[["pam"]],
                                 clade.sizes = clades[["sizes"]],
                                 sample.dists = sample.dists, n.cores = n.cores)

        # test for random structural effects on every response pattern y with Bayesian posterior distributions
        if (has.idd) {
            pat.h0 <- list(lmms.pat[["dif"]][["h0"]], lmms.pat[["idd"]][["h0"]])
        } else {
            pat.h0 <- list(lmms.pat[["dif"]][["h0"]])
        }
        struc.eff <- testForStruEff(pat.h0 = pat.h0, Y = alleles[["Y"]],
                                    C = C[["C"]], K = C[["K"]],
                                    L = ncol(snps[["G"]]), n.cores = n.cores)

        # get correlations between clades and projections
        # Notice clades[["pam]] may be derived from a phylogenetic tree rather than a projection tree.
        r <- corCladeProj(clades = clades[["pam"]], projections = C[["C"]],
                          clade.sizes = clades[["sizes"]], n.cores = n.cores)

        # create the outcome of this stage
        struc <- list(C = C, clades = clades, mc = mc.xy, eff = struc.eff,
                      cor = r, tree = tree)  # save the tree actually used for structural analysis
        if (save.stages) {
            saveRDS(struc, file = stage.record)
        }
    }

    # No clearance step in this function as the function findPhysLink usually follows.

    # Finally, construct and return an interface for users or other functions
    return(list(outputs = outputs, stage.outputs = stage.outputs,
                snps = snps, alleles = alleles, genes = genes, mapping = gene.alleles,
                tests = allele.pairs, lmms.pat = lmms.pat, lmms = lmms, struc = struc))
}
