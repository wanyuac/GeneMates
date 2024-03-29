% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func__distanceAnalyser__summariseDist.R
\name{summariseDist}
\alias{summariseDist}
\title{Attach new columns summarising distance measurements to LMM results (lmms)}
\usage{
summariseDist(
  lmms,
  lmms.ds,
  tree,
  clades,
  allele.pam,
  max.nodes = -1,
  max.dist = Inf,
  colname.pair = "pair",
  colname.co = "n_xy",
  qs = c(0, 0.25, 0.5, 0.75, 1),
  n.cores = n.cores
)
}
\arguments{
\item{lmms}{a data frame that must contain two columns x and y for allele names.}

\item{lmms.ds}{a data frame sharing the columns x and y of lmms and all of
their distance measurements.}

\item{tree}{A phylogenetic tree (an object of the phylo class defined in the
ape package) accompanying the clades argument.}

\item{clades}{A PAM of samples in every clade. It can be obtained using the
function tree2Clade.}

\item{allele.pam}{an allelic PAM of bacterial genes. Usually it is the element
A in the output list of the function importAllelicPAM.}

\item{max.nodes}{an integer specifying the maximal number of nodes to which
we can trust paths for distance measurements. Set max.nodes = -1 to skip this
fiter when the accuracy of measurements does not matter. Notice node number
equals zero for overlapping genes.}

\item{max.dist}{An inclusive upper bound for filterring distance measurements.
Measurements above this threshold will be ignored. Set it to Inf to turn off
this filter.}

\item{colname.pair}{a column name in the data frame lmms for indices of a pair
of alleles}

\item{colname.co}{a column name in the data frame lmms for the number of co-occurrence}

\item{qs}{quantile probabilities (between 0 and 1, inclusive) of physical
distances. The minimum, first quantile (25th percentile), median, third quantile
(75th percentile) and the maximum will always be calculated.}

\item{n.cores}{Number of computational cores that will be used in parallel for
this function.}
}
\description{
This function produces summary statistics of distance measurements
per pair of alleles.
}
\author{
Yu Wan (\email{wanyuac@126.com})
}
