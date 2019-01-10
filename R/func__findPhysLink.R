#' @title Testing for associations and physical linkage between every pair of
#' bacterial genes at the allele level
#'
#' @description This function has two behaviours. It always tests for associations
#' between alleles of bacterial genes with univariate linear mixed models. In
#' addition, when allelic physical distances are provided, it determines evidence
#' of physical linkage between the alleles. This is the main function of the
#' package GeneMates. For this package, the term "sample" refers to either a
#' bacterial isolate or a strain.
#'
#' Dependent libraries: data.table, parallel, ape, phytools
#'
#' @param assoc.out A previous output of findPhysLink when it was used only for
#' association tests. Equivalently, this output is the same as the function lmm.
#' This list may not include the large element snps for convenience. This parameter
#' is useful when users need to incorporate distance information into the result
#' of association analysis. Other arguments, such as snps, snps.delim, ...,
#' allelic.pam and genetic.pam, etc., will not be used when this argument is valid.
#
#  =============== Parameters for the lmm function ===============
#' @param snps Core-genome SNPs used for estimating the relatedness matrix.
#' Valid values: a complete path to a SNP table, or an un-centred, encoded,
#' biallelic SNP matrix G
#' @param snps.delim (optional) Delimiters of fields in the SNP table (pam.delim
#' and dist.delim are defined similary)
#' @param pos.col (optional) An integer (column index) or a string (column name)
#' specifying which column contains SNP positions
#' @param ingroup (optional) A vector of characters for names of isolates to be
#' analysed. Isolates may be sorted, such as according to the phylogeny. The function
#' includes all isolates by default.
#' @param outliers (optional) A vector of characters for isolate/strain names to
#' be excluded from snps, pam and ds
#' @param ref (optional) name of the reference genome. The column name "Ref" in
#' the SNP matrix will be replaced with this argument.
#' @param min.mac (optional) An integer specifying the minimal number of times
#' required for the minor allele of every biallelic SNP to occur across all isolates.
#' SNPs failed this criterion will be removed from this analysis.
#' @param genes.excl (optional) Genes to be excluded from PAMs. For example,
#' genes.excl = c("AmpH_Bla", "OqxBgb_Flq", "OqxA_Flq", "SHV.OKP.LEN_Bla").
#' @param genetic.pam A presence/absence matrix of genes. It may be a compiled
#' table from SRST2.
#' @param allelic.pam A presence/absence matrix of alleles. The matrix may be a
#' compiled table of SRST2 results.
#' @param genetic.pam.delim (optional) A delimiter character in the genetic PAM.
#' Default: tab.
#' @param allelic.pam.delim (optional) A delimiter character in the allelic PAM.
#' Default: tab.
#' @param min.count (option) The minimum count of alleles/genes in the current
#' data set to be included for analysis.
#' @param min.co (optional) The minimal number of allelic co-occurrence events.
#' @param mapping (optional) A data frame mapping alleles to genes and patterns,
#' etc, which equals the "mapping" element within findPhysLink's output list.
#'      This argument is only used when a user reruns a previous analysis.
#' @param tree (optional) A path to a tree file or a phylo object for a tree of
#' all samples.
#' The format of the tree file must be compartile to the read.tree function in
#' the ape package.
#' @param sample.dists (optional) A numeric matrix of distances between samples.
#' The distances can be Euclidean distances between projections, phylogenetic
#' tip distances or SNP distances (the number of SNPs between any two samples).
#' A matrix of Euclidean distances will be computed for projections
#' of samples if this option is left NULL.
#' @param gemma.path Path to GEMMA. No forward slash should be attached at the
#' end of the path.
#
#  =============== Parameters for findPhysLink ===============
#' @param phys.dists A table of physical distances between targets. GeneMates
#' matches sample names between this table and those in the SNP matrix, hence it
#' is okay to have extra or less samples in phys.dists.
#' @param dist.delim The delimit character in the table of physical distances.
#' Default: tab.
#' @param max.node.num (optional) An integer specifying the maximal number of
#' nodes per path in which we can trust their distance measurements.
#'      Set max.node.num = NULL to turn off this filter of distance measurements.
#' @param max.dist (optional) An inclusive upper bound for filterring distance
#' measurements. Measurements above this threshold will be ignored.
#'      Set it to NULL to turn off this filter.
#' @param d.qs (optional) Quantile probabilities for allelic distance measurements.
#' Default: minimum (0), the first quantile (0.25), median (0.5), the third
#' quantile (0.75) and the maximum (1).
#' @param max.p (optional) The upper bound of P values to determine significance.
#' P <= max.p will be called significant.
#' @param max.range An integer specifying the maximum range (in bp) of in-group
#' allelic physical distances allowed to call consistent.
#' @param min.pIBD (optional) Minimum probability of the root of a minimum inclusive
#' clade to display a positive binary trait, such as having a pair of alleles co-occurring
#' or a specific allelic physical distance. Default: 0.9 (90\%).
#' @param output.dir (optional) Path of the output directory.
#' A relative path "output" under the current working directory is recommended
#' as GEMMA always create a directory named output to store its outputs.
#' Otherwise, you will end up with two output directories: one for yours, and
#' the other for GEMMA.
#' @param prefix (optional) For names of all output files
#' @param n.cores Number of cores used to run GEMMA in parallel where possible.
#' -1: automatically detect the number of available cores N, but use N - 1 cores (recommended)
#' 0: automatically detect the number of available cores and use all of them.
#' Be careful when the current R session is not running through SLURM.
#' >= 1: use the number of cores as specified. n.cores is reset to the maximal
#' number of available cores N when n.cores > N.
#' @param save.stages (optional) Whether to turn on stage control or not.
#' Recommend to turn it on when you are not sure whether the pipeline will be
#' finished smoothly.
#' @param del.temp (optional) Whether to delete temporary files or not.
#' @param skip (optional) Whether to avoid overwriting existing output files.
#'
#' @examples
#' time.start <- Sys.time()
#' assoc <- findPhysLink(snps = "input/noPhage_snps_1outgroup_var_regionFiltered_cons1.csv",
#' snps.delim = ",", pos.col = "Pos", allelic.pam = "input/allele_paMatrix_noHash_filtered.txt",
#' genetic.pam = "input/modified_allele_matrix_noHash_filtered.txt",
#' genes.excl = c("AMPH_Ecoli_Bla", "AmpC1_Ecoli_Bla", "AmpC2_Ecoli_Bla", "MrdA_Bla"),
#' phys.dists = "input/merged_dists_noHash.tsv", max.node.num = 2,
#' max.dist = 2.5e6, max.range = 2000, min.pIBD = 0.9, ingroup = NULL, outliers = "Outgroup",
#' min.mac = 1, min.co = 2, max.p = 0.05, output.dir = "output", prefix = "Ec",
#' gemma.path = "~/apps/gemma", n.cores = 8)
#' time.end <- Sys.time()
#' print(time.end - time.start)
#'
#' @author Yu Wan, \email{wanyuac@@gmail.com}
#' @export
#
#  Copyright 2017-2018 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 17 March 2017, the lastest edition: 10 Jan 2019

findPhysLink <- function(assoc.out = NULL,
                         snps = NULL, snps.delim = ",", pos.col = "Pos", min.mac = 1,
                         genetic.pam = NULL, genetic.pam.delim = "\t", genes.excl = NULL,
                         allelic.pam = NULL, allelic.pam.delim = "\t",
                         min.count = 2, min.co = 2, mapping = NULL,
                         phys.dists = NULL, dist.delim = "\t",
                         max.node.num = NULL, max.dist = NULL,
                         ingroup = NULL, outliers = NULL, ref = NULL,
                         tree = NULL, sample.dists = NULL,
                         d.qs = c(0, 0.25, 0.5, 0.75, 1),
                         max.p = 0.05, max.range = 2000, min.pIBD = 0.9,
                         output.dir = "output", prefix = NULL,
                         gemma.path = "gemma", n.cores = -1,
                         save.stages = TRUE, del.temp = TRUE, skip = TRUE) {

    # 1. Run the lmm function for association tests ###############
    print(paste0(Sys.time(), ": starting the findPhysLink function."))

    # Skip association tests if the output list of the function lmms is provided.
    # Notice it is not necessary to have the large element "snps" in the list assoc.out.
    if (! .checkAssocOut(assoc.out)) {
        print("The argument assoc.out will not be used by this function as it is not a valid result of association analysis.")
        stage.record <- paste0("temp/", prefix, "__assoc_out.rds")
        if (save.stages & file.exists(stage.record)) {
            assoc.out <- .recoverHistory(stage.record)
        } else {
            assoc.out <- lmm(snps = snps, snps.delim = snps.delim, pos.col = pos.col,
                             min.mac = min.mac, genetic.pam = genetic.pam,
                             genetic.pam.delim = genetic.pam.delim,
                             genes.excl = genes.excl, allelic.pam = allelic.pam,
                             allelic.pam.delim = allelic.pam.delim,
                             min.count = min.count, min.co = min.co,
                             mapping = mapping, ingroup = ingroup, outliers = outliers,
                             ref = ref, tree = tree, sample.dists = sample.dists,
                             output.dir = output.dir, prefix = prefix,
                             gemma.path = gemma.path, n.cores = n.cores,
                             save.stages = save.stages, skip = skip)
            if (save.stages) {
                saveRDS(assoc.out, file = stage.record)
            }
        }
    } else {
        print("The argument assoc.out will be used by this function as it is a valid result of association analysis.")
        .checkDir(output.dir)
    }

    # Cannot use attach(assoc.out) to unpack the list due to name clashes with this function's arguments.
    stage.outputs <- .initialiseStageRecords("findPhysLink", save.stages, prefix)  # ./temp/...

    # 2. Incorporate allelic physical distances into LMM results ###############
    stage.record <- stage.outputs[["ds"]]
    if (save.stages & file.exists(stage.record)) {
        records <- .recoverHistory(stage.record)
        ds <- records[["ds"]]
        lmms.ds <- records[["lmms.ds"]]
        score.dists <- TRUE
        rm(records)
    } else if (!is.null(phys.dists)) {
        ds <- .importPhysicalDists(dists = phys.dists, delim = dist.delim,
                                   ingroup = rownames(assoc.out[["alleles"]][["A"]]),
                                   outgroup = outliers)  # a data frame of original distance measurements
        lmms.ds <- .attachDistances(assoc.out[["lmms"]], ds)  # attach every distance measurement to allele pairs within the LMM outputs fitted under the alternative hypothesis H1
        lmms.ds <- .retrievePairID(lmms.ds, assoc.out[["lmms"]])  # attach a pair ID for every pair of queries
        if (save.stages) {
            saveRDS(list(ds = ds, lmms.ds = lmms.ds), file = stage.record)
        }
        score.dists <- TRUE
    } else {
        print("The function only perform LMM-based association tests because no allelic physical distance is provided.")
        score.dists <- FALSE
    }

    # 3. Scoring allelic physical distances ###############
    if (score.dists) {
        # 1) Infer identical-by-descent of measured distances and calculate slopes of the measurements ===============
        # return summary statistics of distance measurements per allele pair
        # Elements on the input list lmms will be concatenated into a single data frame as the output.
        stage.record <- stage.outputs[["ds.summary"]]
        if (save.stages & file.exists(stage.record)) {
            ds.stats <- .recoverHistory(stage.record)
        } else {
            ds.stats <- summariseDist(lmms = assoc.out[["lmms"]], lmms.ds = lmms.ds,
                                      tree = assoc.out[["struc"]][["tree"]],
                                      clades = assoc.out[["struc"]][["clades"]],
                                      allele.pam = assoc.out[["alleles"]][["A"]],
                                      max.nodes = max.node.num, max.dist = max.dist,
                                      colname.pair = "pair", colname.co = "n_xy",
                                      qs = d.qs, n.cores = n.cores)
            if (save.stages) {
                saveRDS(ds.stats, file = stage.record)
            }
        }


        # 2) Evaluate evidence of physical linkage for every pair of alleles ===============
        stage.record <- stage.outputs[["assoc"]]
        if (save.stages & file.exists(stage.record)) {
            assoc <- .recoverHistory(stage.record)  # In practice, it is very unlikely that a user wants to recover the result of this stage.
        } else {
            # Concatenate data frames and score evidence of physical linkage
            assoc <- evalPL(lmms = ds.stats, min.beta = 0, max.p = max.p,
                            max.range = max.range, min.pIBD = min.pIBD)

            # Merge structural random effects into the result table
            assoc <- merge(x = assoc,
                           y = assoc.out[["struc"]][["mc"]][, c("y", "x", "clade",
                                                                "size", "f_xy", "ds_max")],
                           by = c("y", "x"), all = TRUE, sort = FALSE)
            leading.columns <- c("pair", "y", "x", "y_pat", "x_pat", "dif", "n_y",
                                 "n_x", "n_xy", "m", "m_in", "score", "s_a", "w_d",
                                 "s_d", "beta", "p_adj", "l_remle", "clade", "size",
                                 "ds_max", "f_xy")
            assoc <- assoc[, c(leading.columns, setdiff(names(assoc), leading.columns))]  # rearrange columns of assoc
            assoc <- assoc[order(assoc$pair, assoc$y, decreasing = FALSE), ]
            if (save.stages) {
                saveRDS(assoc, file = stage.record)
            }
        }

        # 3) Expand the return object when physical distances are provided ===============
        assoc.out[["ds"]] <- ds
        assoc.out[["lmms.ds"]] <- lmms.ds
        assoc.out[["ds.stats"]] <- ds.stats
        assoc.out[["assoc"]] <- assoc
        assoc.out[["stage.outputs"]] <- append(assoc.out[["stage.outputs"]], stage.outputs)
    }

    # 4. Clearance ###############
    if (save.stages & del.temp) {
        print("Deleting temporary files.")
        system("rm -rf temp")
    }
    print(paste0(Sys.time(), ": the findPhysLink function has been implemented successfully."))

    return(assoc.out)
}
