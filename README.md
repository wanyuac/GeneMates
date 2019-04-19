# GeneMates: an R package identifying horizontal gene co-transfer between bacteria

GeneMates is an R package implementing a network approach to identify horizontal gene co-transfer (HGcoT) between bacteria using whole-genome sequencing (WGS) data. It is particularly useful for investigating intra-species HGcoT, where presence-absence status of acquired genes is usually confounded by bacterial population structure due to clonal reproduction.

**Citation**

Since our article for this package is under preparation, the package and helper tools should be used with caution. Please cite this repository if you use our tool:

**Y. Wan, R.R. Wick, J. Zobel, M. Inouye, D. Ingle, K.E. Holt**, GeneMates: an R package detecting horizontal co-transfer of bacterial genes, GitHub: github.com/wanyuac/GeneMates, 2018.

This project is supported by the Department of Biochemistry and Molecular Biology, University of Melbourne.

## Table of Contents

1. [Installation](#installation)  
    * 1.1. [Dependencies](#dependencies)  
    * 1.2. [Helper scripts](#helpers)  
2. [Quick start](#quickStart)  

## <a name="installation">1. Installation</a>

There are two approaches to install this package. Take the GeneMates version 0.1.6 for example. Assuming that we are going to install the package under a user-specified directory Lib:

**R**

```R
install.packages(pkgs = "GeneMates_0.1.6.tar.gz", lib = "Lib")
```

**bash**

The program Rscript should be accessible as a command. Namely, the path of R should be added to $PATH before hand.

````bash
./install_GeneMates.sh GeneMates_0.1.6.tar.gz Lib 
````

## <a name="dependencies">1.1. Dependencies</a>

### Software

* [R](https://www.r-project.org) ≥ 3.3.3
* [GEMMA](https://github.com/genetics-statistics/GEMMA) 0.96
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) ≥ 2.2.30

### R packages

* parallel (≥ 3.3.3)
* data.table (≥ 1.10.4)
* ape (≥ 4.1)
* phytools (≥ 0.6-00)
* logistf (≥ 1.22)
* pheatmap (≥ 1.0.8)
* ggplot2 (≥ 2.2.1)
* [ggtree](https://github.com/GuangchuangYu/ggtree)(≥ 1.6.11)
* network (≥ 1.13.0.1)
* networkDynamic (≥ 0.9.0)

## <a name="helpers">1.2. Helper scripts</a>

We have developed the following scripts to help users to prepare input data for GeneMates:  

1. Screening genomic data for known genes  

   - [PAMmaker](https://github.com/wanyuac/PAMmaker), which compiles outputs of [SRST2](https://github.com/katholt/srst2) and geneDetector.  
   - [geneDetector](https://github.com/wanyuac/geneDetector) for genome assemblies.  

2. [readSimulator](https://github.com/wanyuac/readSimulator) for simulating short reads from template DNA sequences.  

3. [cgSNPs](https://github.com/wanyuac/cgSNPs "cgSNPs") for processing core-genome SNPs (cgSNPs).  

4. Measurement and compilation of allelic physical distances (APDs)  

    - [Bandage](https://github.com/rrwick/Bandage) ≥ 0.8.1 for distance measurement  
    - [physDist](https://github.com/wanyuac/physDist "physDist") for compiling APDs  

## <a name="quickStart">2. Quick start</a>

Function _findPhysLink_ is a pivotal function and a major user interface of this package. Here we show a typical procedure for the usage of GeneMates. We assume this procedure is run on a computing server.

```R
setwd("Analysis")

library(GeneMates)

tr <- read.tree("rooted.tree")  # a maximum-likelihood tree

assoc <- findPhysLink(snps = "Analysis/snps.csv",
                      snps.delim = ",", pos.col = "Pos", min.mac = 1,
                      genetic.pam = "gpam.tsv", genetic.pam.delim = "\t",
                      allelic.pam = "apam.tsv", allelic.pam.delim = "\t",
                      mapping = NULL, min.count = 15,
                      phys.dists = "prioritised_dists.tsv",
                      dist.delim = "\t", max.node.num = 2, max.dist = 250e3,
                      ref = "ref", tree = tr,
                      min.co = 2, d.qs = c(0, 0.25, 0.5, 0.75, 1), max.p = 0.05,
                      max.range = 2000, min.pIBD = 0.9,
                      output.dir = "Output", prefix = "demo", gemma.path = "~/apps/gemma",
                      n.cores = 8, save.stages = TRUE, del.temp = FALSE, skip = TRUE)

snps <- assoc[["snps"]]  # a large list
saveRDS(snps, file = "Out/snps.rds")  # Analysis/Out/snps.rds

assoc <- assoc[-which(names(assoc) == "snps")]
saveRDS(assoc, file = "Out/assoc.rds")  # Analysis/Out/assoc.rds
```

The element "snps" in the result list is usually too large to be loaded to an R session when the sample size or SNP number is large. Therefore we recommend to save the result list in two files.

