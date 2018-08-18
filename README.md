# GeneMates: an R package performing association analysis for bacterial allelic presence-absence data and inferring horizontal gene co-transfer

This is an ongoing project looking for horizontal gene co-transfer (HGcoT) between bacteria of the same species using whole-genome sequencing (WGS) data. This repository provides an R package called GeneMates to perform the main analysis.

Since this project is in active development and a paper is under preparation, the package and supplementary tools should be used with caution. Please cite this repository if you use our tool:

**Y. Wan, R.R. Wick, J. Zobel, M. Inouye, D. Ingle, K.E. Holt**, GeneMates: an R package detecting horizontal co-transfer of bacterial genes, GitHub: https://github.com/wanyuac/GeneMates, 2018.

This project is supported by the University of Melbourne.

# Installation #

There are two approaches to install this package. Take the GeneMates version 0.1.6 for example.

### R ###

```
install.packages(pkgs = "GeneMates_0.1.6.tar.gz", lib = "Lib")
```

### bash ###
The R program should be accessible as a command. Namely, the path of R should be added to $PATH before hand. 

````bash
./install_GeneMates.sh GeneMates_0.1.6.tar.gz ~/R_lib
````

## Components ##

In a narrow sense, GeneMates is an R package; in a broad sense, GeneMates is a pack of the R package and supplementary tools. A complete GeneMates installation consists of the following components:  
  
1. R package GeneMates (this repository)  
2. Gene screen
  - [PAMmaker](https://github.com/wanyuac/PAMmaker "PAMmaker")
  - [screen_genes_in_assemblies](https://github.com/wanyuac/screen_genes_in_assemblies "screen_genes_in_assemblies")
3. [cgSNPs](https://github.com/wanyuac/cgSNPs "cgSNPs") for processing core-genome SNPs (cgSNPs)
4. [physDist](https://github.com/wanyuac/physDist "physDist") for generating allelic physical distances (APDs).  

It is recommended to install all components under the same parental directory.

## Dependencies ##

### Software ###

* [R](https://www.r-project.org) >= v3.3.3
* [GEMMA](https://github.com/genetics-statistics/GEMMA)
* [Bandage](https://github.com/rrwick/Bandage) >=v0.8.1
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) >= 2.2.30
 
### R packages ###

* parallel (>= 3.3.3)
* data.table (>= v1.10.4)
* ape (>= v4.1)
* phytools (>= v0.6-00)
* logistf (>= v1.22)
* pheatmap (>= v1.0.8)
* ggplot2 (>= v2.2.1)
* [ggtree](https://github.com/GuangchuangYu/ggtree)(>= v1.6.11)
* network (>= 1.13.0.1)
* networkDynamic (>= 0.9.0)

[To be continued]