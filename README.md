# GeneMates: an R package performing association analysis for bacterial allelic presence-absence data and inferring horizontal gene co-transfer

This is an ongoing project looking for horizontal gene co-transfer (HGcoT) between bacteria of the same species using whole-genome sequencing (WGS) data. This repository provides an R package called GeneMates to perform the main analysis.

Since this project is in development and a paper is under preparation, the package and supplementary tools should be used with caution. Please cite this repository if you use our tool:

**Y. Wan, R.R. Wick, J. Zobel, M. Inouye, D. Ingle, K.E. Holt**, GeneMates: an R package detecting horizontal co-transfer of bacterial genes, GitHub: https://github.com/wanyuac/GeneMates, 2018.

This project is supported by the University of Melbourne.

## Installation

There are two approaches to install this package. Take the GeneMates version 0.1.6 for example.

### R

```
install.packages(pkgs = "GeneMates_0.1.6.tar.gz", lib = "Lib")
```

### bash
The R program should be accessible as a command. Namely, the path of R should be added to $PATH before hand. 

````bash
./install_GeneMates.sh GeneMates_0.1.6.tar.gz ~/R_lib
````

## Components

In a narrow sense, GeneMates is an R package; in a broad sense, GeneMates is a pack of the R package and supplementary tools. A complete GeneMates installation consists of the following components:  

1. R package GeneMates (this repository)  
2. Gene screen
	- [PAMmaker](https://github.com/wanyuac/PAMmaker), which compiles outputs of [SRST2](https://github.com/katholt/srst2) and geneDetector.
	- [geneDetector](https://github.com/wanyuac/geneDetector) for genome assemblies.
4. [readSimulator](https://github.com/wanyuac/readSimulator) for simulating short reads from template DNA sequences.
3. [cgSNPs](https://github.com/wanyuac/cgSNPs "cgSNPs") for processing core-genome SNPs (cgSNPs).
4. [physDist](https://github.com/wanyuac/physDist "physDist") for generating allelic physical distances (APDs).  

It is recommended to install all components under the same parental directory.

## Dependencies

### Software

* [R](https://www.r-project.org) >= 3.3.3
* [GEMMA](https://github.com/genetics-statistics/GEMMA) 0.96
* [Bandage](https://github.com/rrwick/Bandage) >= 0.8.1
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) >= 2.2.30

### R packages

* parallel (>= 3.3.3)
* data.table (>= 1.10.4)
* ape (>= 4.1)
* phytools (>= 0.6-00)
* logistf (>= 1.22)
* pheatmap (>= 1.0.8)
* ggplot2 (>= 2.2.1)
* [ggtree](https://github.com/GuangchuangYu/ggtree)(>= 1.6.11)
* network (>= 1.13.0.1)
* networkDynamic (>= 0.9.0)

## Scope
This tool was designed for identifying horizontally co-transferred accessory antimicrobial resistance genes (ARGs) in bacteria of the same species. To this end, it tests for associations between alleles of the ARGs. In theory, it is applicable to other kinds of bacterial genes when the following assumptions are satisfied. 

1. Bacteria are isolated within a short period, in which it is unlikely to have mutations in ARGs.
2. There is only a single allele per ARG per cell.

## Usage

Function _findPhysLink_ is the pivotal function of this package. It is an integration of several functions.

### Association analysis controlled for bacterial population structure

[To be continued]

### Inference of HGcoT

[To be continued]