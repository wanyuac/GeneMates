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
3. [Function hierarchy](#scriptOrganisation)  
    * 3.1. [Main function](#mainFunction)  
    * 3.2. [Network analysis](#networkAnalysis)  
    * 3.3. [Summary statistics of input data](#summaryStats)  
    * 3.4. [Analysis of population structural](#svd)  
    * 3.5. [Identification of mobile gene clusters](#MGEdiscovery)  
    * 3.6. [Spatiotemporal analysis](#spatiotemporal)  
    * 3.7. [Visualisation](#visualisation)  
    * 3.8. [Helper functions](#helpers)  
    * 3.9. [Alternative association analysis](#plr)  
4. [Inputs](#inputs)  
    * 4.1. [Import core-genome SNP matrix](#cgSNPs)  
5. [References](#references)

## <a name="installation">1. Installation</a>

There are two approaches to install this package. Take the GeneMates version 0.1.6 for example. Assuming that we are going to install the package under a user-specified directory Lib:

**R**

```R
install.packages(pkgs = "GeneMates_0.2.0.tar.gz", lib = "Lib")
```

**bash**

The program Rscript should be accessible as a command. Namely, the path of R should be added to $PATH before hand.

````bash
./install_GeneMates.sh GeneMates_0.2.0.tar.gz Lib 
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

3. [cgSNPs](https://github.com/wanyuac/cgSNPs) for processing core-genome SNPs (cgSNPs).  

4. Measurement and compilation of allelic physical distances (APDs)  

    - [Bandage](https://github.com/rrwick/Bandage) ≥ 0.8.1 for distance measurement  
    - [physDist](https://github.com/wanyuac/physDist) for compiling APDs  

## <a name="quickStart">2. Quick start</a>

Function _findPhysLink_ is a pivotal function and a major user interface of this package. Here we show a typical procedure for the usage of GeneMates. We assume this procedure is run on a computing server.

```R
setwd("Analysis")

library(GeneMates)

tr <- read.tree("rooted.tree")  # a maximum-likelihood tree

assoc <- findPhysLink(snps = "Analysis/snps.csv", snps.delim = ",", pos.col = "Pos", ref.pos = "Ref", min.mac = 1, genetic.pam = "gpam.tsv", genetic.pam.delim = "\t", allelic.pam = "apam.tsv", allelic.pam.delim = "\t", mapping = NULL, min.count = 15, phys.dists = "prioritised_dists.tsv", dist.delim = "\t", max.node.num = 2, max.dist = 250e3, ref = "ref", tree = tr, min.co = 2, d.qs = c(0, 0.25, 0.5, 0.75, 1), max.p = 0.05, max.range = 2000, min.pIBD = 0.9, output.dir = "Output", prefix = "demo", gemma.path = "~/Apps/gemma", n.cores = 8, save.stages = TRUE, del.temp = FALSE, skip = TRUE)

snps <- assoc[["snps"]]  # a large list
saveRDS(snps, file = "Out/snps.rds")  # Analysis/Out/snps.rds

assoc <- assoc[-which(names(assoc) == "snps")]
saveRDS(assoc, file = "Out/assoc.rds")  # Analysis/Out/assoc.rds
```

The element "snps" in the result list is usually too large to be loaded to an R session when the sample size or SNP number is large. Therefore we recommend to save the result list in two files.

## <a name="scriptOrganisation">3. Function hierarchy</a>

This section illustrates organisation of public functions in GeneMates.

### <a name = "mainFunction">3.1. Main function</a>

*findPhysLink*: the main function of GeneMates.  

- *lmm*: fits linear mixed models (LMMs) for allelic presence-absence status. It has the following subordinate functions:  
  - *importCoreGenomeSNPs*: reads and processes a core-genome SNP (cgSNP) table.
  - *importGeneticPAM*: reads and processes a genetic presence-absence matrix (PAM).
  - *importAllelicPAM*: reads and processes an allelic PAM.
  - *countAlleles*: counts the occurrence of each allele across bacterial isolates and summarises the number of alleles per gene.
  - *assignPairID*: assigns an identifier to each pair of linear models Y ~ X and X ~ Y.
  - *projectSamples*: performs sample projection via singular-value decomposition.
  - *tree2Clades*: converts a user-specified or projection tree into a presence-absence matrix of isolates in each clade.
  - *findMinIncClade*: determines the minimum inclusive clade (MIC) of each allele across the tree.
  - *testForStruEff*: hypothesis tests for structural random effects.
  - *corCladeProj*: estimates correlation between clades and sample projections.
- *summariseDist*: summarises APDs.  
- *evalPL*: evaluates evidence of physical linkage between alleles and scores edges in the output network.  

### <a name = "networkAnalysis">3.2. Network analysis</a>

The following functions analyse the output network of *findPhysLink*.  

- *mkNetwork*: creates a network object (Graph) from the result of *findPhysLink*.  
- *compileGraphs*: compiles graphs into a single network while retaining separation between individual graphs.  
- *summariseCliques*: summarises allele content per clique.  
- *summarisePhysDistForClusters*: summarises APDs for allele clusters.  
- *summariseDistsForEdges*: summarises APDs for each edge.  
- *betaSignCongruence*: compares signs of fixed effects estimated using LMMs and penalised logistic models (PLMs).  
- *mergeIddAlleles*: merges nodes representing identically distributed alleles into a single node in an output network of *findPhysLink*.  
- *alleleClusterDistr*: computes a matrix for presence of allele clusters based on allelic presence-absence status.  
- *allelicCoMatrix*: create a co-occurrence matrix from allelic PAM.  
- *compEdgeOccur*: compare similarity between edges' occurrence across bacterial isolates.  
- *countGeneClassLinks*: counts the number of links in a network at the level of gene classes.  
- *countNeighbours*: count neighbours of each node in a network.  
- *extractSubgraphs*: creates objects of the class "Graph" for subgraphs of a large graph.  
- *getClusterMemberCooccurrence*: retrieves member allele co-occurrence status and isolate distributions for clusters of alleles.  

### <a name = "smmaryStats">3.3. Summary statistics of input data</a>

- *calcGeneFreq*: calculates gene frequency using the output of *countAlleles*.  
- *corPatternProj*: estimates the correlation between allele distribution patterns and sample projections. This function is similar to *corCladeProj*.  

### <a name = "svd">3.4. Analysis of population structural</a>

- *screePlotPCs*: calculates proportion of genetic variation captured by each principal component.  
- *projectSamples*: see [Section 3.1](#mainFunction).  
- *testForStruEff*: as above.  
- *tree2Clades*: as above.  

### <a name = "MGEdiscovery">3.5. Identification of mobile gene clusters</a>

- *findSeq*: searches for a query sequence against a list of assemblies. This function requires accessibility to [Bandage](https://github.com/rrwick/Bandage).  

### <a name = "spatiotemporal">3.6. Spatiotemporal analysis</a>

GeneMates comprises functions analysing spatial and temporal distributions of bacterial isolates and alleles of genes of interest.

- Spatial analysis  
  - *countAllelesPerCountry*: counts alleles of interest in each country.  
  - *countAllelesPerGeneByCountry*: counts every allele of given genes in each country.  
  - *calcAllelicDiveristyPerCountry*: calculates Simpson's or Shannon's diversity index for alleles per country or any other kind of geographic regions. The function assumes that each strain must belong to only a single country or region.  
- Temporal analysis  
  - *countAllelesPerYear*: count alleles per year.  
  - *countAllelesPerGeneByYear*: counts every allele of given genes per year.  
  - *getAllelesEarliestAppearance*: identifies the earliest appearance of each allele specified by a vector.  
  - *mkCoocurNetwork*: converts a network produced by *findPhysLink* into a co-occurrence network.  
  - *tempNet*: creates a temporal network from graphs (such as a co-occurrence network).  

### <a name = "visualisation">3.7. Visualisation</a>

- *ringPlotPAM*: makes a ring plot to show presence-absence of genotypes and allelic co-occurrence.  
- *heatMapPAM*: an expansion of ggtree's gheatmap function for displaying an allelic PAM.  
- *drawHeatMap*: a generic function creating a heat map for a given variable.  
- *compPvalues*: draws a scatter plot and histograms to compare p-values from LMMs and PLMs.  
- *showGeneContent*: draws a bubble plot and two bar plots to summarise gene and allele frequencies.  

### <a name = "helpers">3.8. Helper functions</a>

Functions under this category are developed for helping users to extract and inspect specific aspects of results from *findPhysLink*.

- *getRowsXY*: a generic function used for retrieves rows in a data frame through a pair of keys.  
- extractPairedRows*: splits a data frame into two data frames for paired rows (that is, X ~ Y and Y ~ X) and unpaired rows (that is, only X ~ Y or Y ~ X exists in the input data frame), respectively.  
- *findMinIncCladeOfStrains*: finds out and summarises the minimal inclusive clade containing all given strains or isolates.  
- *getAllelesPerPattern*: returns a vector of allele names under each distribution pattern.  
- *getAssocGenePair*: searches the association table produced by *findPhysLink* for rows corresponding to alleles of a given pair of genes.  
- *getGeneClass*: extracts gene classe names from a vector of SRST2-formatted allele IDs.  
- *mkFilterTSV*: makes a guidance tab-delimited file as an input for the [physDist](https://github.com/wanyuac/physDist) pipeline.  
- *vertexAttr2Size*: maps a numeric vector to vertex sizes through a linear transformation.  
- *retrieveAlleleSetInfo*: retrives allelic presence-absence information given a vector of allele names.  

### <a name = "plr">3.9. Alternative association analysis</a>

- *plr*: uses Firth's penalised logistic regression rather than LMMs to model allelic presence-absence status. It works in a similar manner as *lmm*.  

## <a name = "inputs">4. Inputs</a>
GeneMates takes as input four kinds of data for detection of HGcoT. A function is created for importing each kind of data. This section explain the usage of these four functions.

### <a name = "cgSNPs">4.1. Import core-genome SNP matrix</a>

####  Input

Function *importCoreGenomeSNPs* reads a cgSNP table, encodes SNP genotypes, extracts biallelic SNPs and performs zero-centring of encoded SNP genotypes. The behaviour of this function is similar to the function *get_SNP_data* in R package [BugWAS](https://github.com/sgearle/bugwas) (see BUGWAS_modular.R). The function *importCoreGenomeSNPs* expects the SNP table to follow the output format of the script parseSNPtable.py in a read-mapping pipeline [RedDog](https://github.com/katholt/RedDog). Accordingly, the SNP table should be stored as a CSV file by default and an example of its structure is shown as follows.

| Pos,Ref,Isolate1,Isolate2,Isolate3,... |
| -------------------------------------- |
| 10,A,A,A,A,...                         |
| 21,C,C,C,C,...                         |
| 25,C,C,C,T,...                         |
| ... |

Here, the hyphen "-" denotes the SNP site at the 21<sup>st</sup> base (Pos) of the reference genome "Ref" is not present in Isolate2.

####  Arguments (excluding those specifying output filenames)

- snps.delim: a single character for the delimiter in the SNP table to be imported. Default: ",".
- pos.col: the name for the position column. Default: "Pos".
- ref.col: the name for the reference column. Default: "Ref".
- replace.ref: new name for the reference column when the column is present in the SNP table.
- ingroup: a character vector of isolate names to be included in the resulting SNP matrix, which may include reference SNP genotypes.
- outliers: a character vector of outlier isolate names, which may include the reference strain. These isolates will be excluded from the SNP matrix. Note that an error arises when ingroup and outliers overlap.
- min.mac: the minimum minor-allele count of each SNP across all isolates excluding outlier isolates. Default: 1, keep every SNP whose minor allele occurs at least once in the remaining isolates.

####  Limitation

Because of the code

```R
snps.var <- apply(snps.core, 2, function(x) length(unique(x)))
snps.bi <- snps.core[, as.integer(snps.var) == 2]
snp.alleles <- .getAllelePairs(snps.bi)
G <- .encodeAlleles(snps.bi, snp.alleles)
```

function *importCoreGenomeSNPs* only works correctly when all SNPs are detected in all ingroup genomes (namely, cgSNPs). In other words, missing SNP genotypes (each is denoted by a hyphen) create false genotype counts. It is not necessary to address this limitation for GeneMates because our estimation of bacterial population structure only relies on biallelic SNP sites that are found in all ingroup isolates.

#### Procedure

- Read the SNP table, replace the reference strain name when possible, exclude outlier isolates and/or keep ingroup isolates.
- Count the number of genotypes per SNP site, identify biallelic SNP sites in the ingroup genomes and determine minor (allele frequency < 0.5) and major alleles of each biallelic SNP site. Note that when both alleles of the same SNP site occur at the same frequency of 0.5, the alleles will be sorted alphabetically (namely, following the order of "A", "C", "G", "T" resulting from the behaviour of the R function *table*) and the first allele will be chosen as the minor allele.
- Encode alleles of each biallelic SNP site in accordance with a convention in genome-wide association analysis (GWAS)<sup>1, 2 (software manual)</sup>.
  - 1: minor allele
  - 0: major allele
- (Optional) filter out biallelic SNPs of insufficient minor-allele frequencies (MAFs).
- MAF-based zero-centring of the remaining biallelic SNP matrix.
- Save codes and positions of biallelic SNPs in the BIMBAM format for [GEMMA](https://github.com/genetics-statistics/GEMMA) (see the [manual](http://www.xzlab.org/software/GEMMAmanual.pdf) of GEMMA for details of this file format) and return a large list to the parental R session.

#### Output

Function *importCoreGenomeSNPs* returns a large list to R. Cautions must be taken to run this function due to the large size of the output list: users are advised to check their computer capacity in the first place. Elements of the output list are listed as follows.

- G: a binary matrix (rows: isolate names, columns: SNP sites) of biallelic SNP sites in the ingroup genomes. This matrix does not contain information of SNP positions.
- S: a column-wise zero-centred SNP matrix G.
- snp.alleles: a matrix of major and minor alleles of SNP sites in G.
- mac: a named integer vector of minor-allele counts of SNP sites in G.
- core: a matrix of unencoded cgSNP genotypes in ingroup genomes (rows denote isolate names and columns denote SNP sites). This matrix includes biallelic SNPs.
- var: a named integer vector of genotype counts across all cgSNPs in ingroup genomes.
- bi: a matrix of biallelic SNPs extracted from the matrix _core_.
- G.bimbam: a BIMBAM-formatted data frame directly converted from the matrix G.
- annots: a BIMBAM-formatted data frame of SNP positions.

## <a name = "references">References</a>

1. McVean, G. A Genealogical Interpretation of Principal Components Analysis. _PLOS Genet_. 5, e1000686 (2009).
2. Zhou, X. & Stephens, M. Genome-wide efficient mixed-model analysis for association studies. _Nat Genet_ 44, 821–824 (2012).