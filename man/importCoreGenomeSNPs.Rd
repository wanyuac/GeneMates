% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func__importCoreGenomeSNPs.R
\name{importCoreGenomeSNPs}
\alias{importCoreGenomeSNPs}
\title{Import genotypes of core-genome SNPs}
\usage{
importCoreGenomeSNPs(
  snps,
  snps.delim = ",",
  pos.col = "Pos",
  ref.col = "Ref",
  replace.ref = NULL,
  min.mac = 1,
  ingroup = NULL,
  outliers = NULL,
  G.file = "",
  annots.file = "",
  skip = TRUE
)
}
\arguments{
\item{snps}{either a path to the SNP table or a list generated by this function previously}

\item{snps.delim}{a single character for the delimiter in the SNP table to be imported}

\item{pos.col}{either a name or an index of the column for SNP positions}

\item{ref.col}{(optional) A string specifying the column for SNPs of the reference genome.}

\item{replace.ref}{Replace the column name specified by ref.col with this argument when ref.col is found amonst column names.
Such a column is seen in the SNP table created using RedDog. A user may want to substitute
it with a genuine strain name to match the reference isolate in the SNP matrix with that in the genetic and allelic matrices.}

\item{min.mac}{the minimal number of times that each minor allele occurs in the population (with outlier isolates excluded)
min.mac = 1: no filter of SNPs by the minor allele count (MAC) is applied.
min.mac = 2: the filter is applied to remove any isolate-specific SNPs.
min.mac > 2: the filter is applied for other purposes.}

\item{ingroup}{a character vector of isolate names to be included in the SNP matrix}

\item{outliers}{a character vector of names for outlier isolates (to be excluded from the SNP matrix)}

\item{G.file}{file path of a bimbam-formatted genotype matrix (G). Keep it an empty character ("" by default) to prevent the function from writing this genotype matrix to the hard disk.}

\item{annots.file}{file path of bimbam-formatted SNP annotations. Keep it "" to prevent the function from writing SNP annotations to the hard disk.}

\item{skip}{whether to skip overwriting existing output files.}
}
\description{
Read the SNP table usually produced by the parseSNPtable.py of RedDog (https://github.com/katholt/RedDog).
The function also compresses the SNP matrix into a pattern matrix. Notice it is the user's responsibility to ensure the
input SNP table only contains 100% conserved SNPs. This function does not check for the conservation status.
}
\author{
Yu Wan (\email{wanyuac@126.com})
}
