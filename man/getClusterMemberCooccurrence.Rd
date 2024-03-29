% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/func__networkAnalyser__getClusterMemberCooccurrence.R
\name{getClusterMemberCooccurrence}
\alias{getClusterMemberCooccurrence}
\title{Retrieving member allele co-occurrence status and strain distributions
for clusters of alleles}
\usage{
getClusterMemberCooccurrence(
  clstr,
  cluster.colname = "cluster",
  allele.colname = "allele",
  pam,
  clade.pam = NULL,
  clade.sizes = NULL,
  sample.dists = NULL,
  n.cores = -1
)
}
\arguments{
\item{clstr}{a data frame consisting of at least two columns: one for cluster
IDs and the other for allele names (comma-delimited). I recommend to use
integers or characters as clster IDs. This argument can be generated by taking
the columns Clique and Alleles in the output of the function summariseCliques.}

\item{cluster.colname}{the column name for cluster IDs in clstr.}

\item{allele.colname}{the column name for allele IDs in clstr.}

\item{pam}{an allelic PAM with columns for allele names and rows for strain names.}

\item{clade.pam}{(optional) a matrix for the presence/absence of samples in
each clade of an input tree. It can be obtained using the function tree2Clades.}

\item{clade.sizes}{(optional) A named vector of integers for the number of
samples in each clade. It can be obtained from the element "sizes" in the
outputs of the function tree2Clades. Optional.}

\item{sample.dists}{(optional) A square matrix for distances between samples.
It can be acquired through the function projectSamples.}

\item{n.cores}{An integer determining the number of cores used for parallel
computing. Valid values are the same as those for the function findPhysLink.}
}
\value{
A data frame of 12 columns. Explanations to some columns: Column 1, named
by the argument cluster.colname, which contains cluster IDs; Column 2, named
by the argument allele.colname, which stores alleles per cluster; size, number
of alleles per cluster; co_max, maximum number of co-occurring alleles in
strains; co_perc, co_max / size * 100%.
}
\description{
This function returns a data frame given a data frame of cluster
IDs and alleles. In each row, a string of strain names (comma-delimited) are
returned if all of the alleles are co-occurring in these strains, otherwise,
an NA is returned because alleles may differ in strains when only a fixed number
(less than the number of all alleles) of alleles are co-occurring in the strains.
Moreover, the function reports the minimal inclusive clade of strains and the
frequency of co-occurrence events in this clade when all alleles are co-occurring.

Specifically, a cluster may refer to a community (also known as a module) or
a clique in a network.
}
\note{
Clade-wise summaries are not returned or only a few summary columns are
returned when the three optional arguments are incomplete (some may be NULL).
}
\examples{
assoc <- findPhysLink(...)
clusters <- ...  # from a network package of your preference
com.co <- getClusterMemberCooccurrence(com = clusters, pam = assoc[["alleles"]][["A"]], cluster.colname = "community",
clade.pam = assoc[["struc"]][["clades"]][["pam"]], clade.sizes = assoc[["struc"]][["clades"]][["sizes"]],
sample.dists = assoc[["struc"]][["C"]][["d"]], n.cores = 2)

Dependency: data.table, parallel

}
\author{
Yu Wan (\email{wanyuac@126.com})
}
