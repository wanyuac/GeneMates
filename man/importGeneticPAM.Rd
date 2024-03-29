% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func__importGeneticPAM.R
\name{importGeneticPAM}
\alias{importGeneticPAM}
\title{Import the genetic presence/absence (p/a) matrix}
\usage{
importGeneticPAM(
  pam,
  pam.delim = "\\t",
  outliers = NULL,
  min.count = 1,
  genes.rm = NULL,
  sample.order = NULL
)
}
\arguments{
\item{pam}{Expect pam to be "modified_allele_matrix.txt" produced by the cdhitFormatter pipeline.
First column of pam consists of isolate names. Other columns: presence/absence profile of every gene
The genetic PAM may contain more alleles than the allelic PAM does, because the latter applies a more stringent criterion to filter data.}

\item{pam.delim}{a single character for the delimiter in the text file of the genetic PAM}

\item{outliers}{a vector of isolate names to be excluded from rows of the PAM}

\item{min.count}{the minimal number of times that a gene occurs in all isolates excluding outliers
By default, only genes that do not occur at all are removed (min.count = 1).
When min.count = 2, isolate-specific genes are removed as well.}

\item{genes.rm}{a vector of gene names to be excluded from columns of the genetic PAM.}

\item{sample.order}{a vector of isolate names for the PAM to be matched with.
It can be used to filter the genetic PAM for in-group isolates.}
}
\description{
This function reads the genetic presence-absence matrix and returns a list of relevant matrices.
It automatically matches isolates to those in the SNP matrix.
}
\author{
Yu Wan (\email{wanyuac@126.com})
}
