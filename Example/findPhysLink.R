# Run function findPhysLink with GEMMA v0.96 for detection of physical linkage
# Copyright 2019 Yu Wan (wanyuac@gmail.com)
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 21 April 2019

Sys.setenv(LANG="en_US.UTF-8")
Sys.setenv(LC_ALL="en_US.UTF-8")

library(ape)
library(GeneMates)

tr <- read.tree("input/tr.tree")

# In theory, users can specify whatever name for the output directory, however, GEMME v0.96 does not
# support any user-specified name for its output directory and the directory will always be named
# "output". As a result, if output.dir = "Output_physLink", then two output directories will appear:
# Output_physLink and output. Therefore, it is recommended that users keep using "output" for the
# argument output.dir.
assoc <- findPhysLink(snps = "input/cgSNPs.csv",
                      snps.delim = ",", pos.col = "Pos", min.mac = 1,
                      ingroup = NULL, outliers = NULL, ref = NULL,
                      genetic.pam = "input/gpam.tsv", genetic.pam.delim = "\t",
                      genes.excl = c("AMPH_Ecoli_Bla", "AmpC2_Ecoli_Bla", "AmpC1_Ecoli_Bla", "MrdA_Bla"),
                      allelic.pam = "input/apam.tsv", allelic.pam.delim = "\t",
                      min.count = 2, min.co = 0, mapping = NULL,
                      tree = tr, sample.dists = cophenetic.phylo(tr),
                      phys.dists = "input/spds.tsv", dist.delim = "\t",
                      max.node.num = 2, max.dist = 250e3,
                      d.qs = c(0, 0.25, 0.5, 0.75, 1),
                      max.p = 0.05, max.range = 2000, min.pIBD = 0.9,
                      output.dir = "output", prefix = "demo", gemma.path = "~/Tools/gemma",
                      n.cores = 4, save.stages = TRUE, del.temp = FALSE, skip = TRUE)

snps <- assoc[["snps"]]  # a list
saveRDS(snps, file = "output/snps.rds")

assoc <- assoc[-which(names(assoc) == "snps")]
saveRDS(assoc, file = "output/physlink.rds")
