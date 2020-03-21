#' @title Project samples onto eigenvectors (of positive eigenvalues) of the semi-positive definite matrix t(S) * S
#'
#' @description First, this function performs a singular-value decomposition of the core-genome SNP (cgSNP)
#' matrix S through an eigen-decomposition of the relatedness matrix K = S * t(S) / L, where L is
#' the number of cgSNPs. It uses GEMMA to carry out the eigen-decompostion. It returns projections
#' of bacterial isolates on eigenvectors (also known as the right-singular vectors of S) of
#' V = t(S) * S. One may use it before the function testProjections of this package.
#'
#' It is mathematically equivalent to perform a singular-value decomposition of S through commands
#' svd(S) or svd(scale(G, center = TRUE)) to calculate projections of samples on singular vectors of
#' positive singular values (equivalently, eigenvalues). However, Phylix uses GEMMA to perform this
#' procedure because it is much faster than R, and the input matrix K is exactly the relatedness matrix
#' used by GEMMA to fit LMMs.
#'
#' We cannot specify the output directory for GEMMA yet. It stores every output under a directory
#' named "output". This directory appears under the current working directory of R. In addition, it is
#' recommended to use a different prefix for outputs of this calculation to that is used for calculating
#' the relatedness matrix. Otherwise, GEMMA will overwrite the previous log file.
#'
#' Output files: output/[prefix].eigenU.txt, output/[prefix].eigenD.txt and [prefix].log.txt
#' under the current working directory of R.
#'
#' @param K a string specifying the path to a file of the relatedness matrix produced using GEMMA
#' @param G path to the uncentred SNP genotype file G, such as [prefix]__G.txt, which has been used by GEMMA
#' to calculate the K matrix. Here, GEMMA uses it again to match samples between K, G and Y (see below), although
#' mathematically the eigen-decomposition does not need G and Y when samples are already matched, which is the
#' case when using Phylix.
#' @param Y path to the "phenotype" file, such as [prefix]__Y.txt under output/gene in Phylix's output
#' GEMMA requires this file to extract the same samples from the K matrix, although samples from K and Y
#' must match in a correct input data set for Phylix. Phylix keeps samples and their order the same with
#' the row names of the SNP matrix S (centred) or G (uncentred).
#' @param L number of core-genome SNPs, which equals the number of columns of SNP matrices G and S.
#' @param samples row names of cgSNP matrices S or G.
#' @param prefix output files are named [prefix]_eigenD.txt and [prefix]_eigenU.txt according to GEMMA's manual.
#' @param get.dists a logical option specifying whether to calculate a distance matrix between projections of samples
#' @param dist.method a character option specifying the method for calculating the sample distances. Valid values are the same
#' as the base::dist function.
#' @param get.tree a logical option specifying whether to return a midpoint-rooted neighbour-joining tree with the distance matrix.
#' The distance matrix is calculated when this option is set, regardless of the get.dists option.
#' @param gemma.path a path to an executable GEMMA.
#'
#' @examples
#' C <- projectSamples(K = "output/Kp_K.cXX.txt", G = "output/snp/Kp__G.txt", Y = "output/gene/Kp__Y.txt",
#' samples = rownames(G), L = ncol(G), prefix = "Kp_SVD", gemma.path = "~/apps/gemma")
#'
#' @author Yu Wan (\email{wanyuac@@126.com})
#' @export
#
#  Copyright 2017 Yu Wan
#  Licensed under the Apache License, Version 2.0
#  First edition: 16 July 2017, the lastest edition: 15 Jan 2018

projectSamples <- function(K, G, Y, L, samples, prefix = "S", get.dists = TRUE, dist.method = "euclidean",
                           get.tree = TRUE, gemma.path) {
    print(paste0(Sys.time(), ": projecting samples into a new space using singular-value decomposition."))

    # sanity check
    if (!.checkFiles(inputs = c(K, G, Y, gemma.path))) {  # check availability of four files
        stop("Input error: some files are not accessible.")
    }

    # set up output file names
    out <- c(U = paste0("output/", prefix, ".eigenU.txt"),
             D = paste0("output/", prefix, ".eigenD.txt"))  # Both files are created by GEMMA.

    # use GEMMA to perform an eigen-decomposition of the positive semi-definite K
    # The function findPhysLink calls GEMMA to compute K = S %*% t(S) / L
    # For Y having multiple columns, the "-n" argument does not change results as
    # GEMMA does not use this piece of information.
    cmd <- paste(gemma.path, "-g", G, "-p", Y, "-n 1 -k", K, "-eigen -o", prefix, sep = " ")
    print(paste0("Performing eigen-decomposition of K using the command ", cmd))
    system(cmd, wait = TRUE)

    # read the relateness matrix K
    # The function findPhysLink has ensured that sample names follow the same order between the SNP matrix (S) and the presence/absence matrix of resistance alleles (X & Y).
    K <- as.matrix(read.table(file = K, header = FALSE, sep = "\t"))
    K <- K[, -ncol(K)]  # remove the redundant last column of NA's (due to an extra tab at the end of each line)
    rownames(K) <- samples
    colnames(K) <- samples

    # read eigenvalues of K
    # eigen-decomposition of K: K = U * E * U^-1, where E is a diagonal matrix of eigenvalues, namely, E <- diag(ev).
    # The file [prefix].eigenD.txt consists of a single column, which corresponds to the principal diagonal of the matrix E
    ev <- as.numeric(read.table(file = out[["D"]], header = FALSE, sep = "\n")$V1)  # a vector of eigenvalues of K
    ev.ord <- order(ev, decreasing = TRUE)  # GEMMA orders eigenvalues in an increasing order
    ev <- ev[ev.ord]  # rearrange eigenvectors from the largest to the smallest following convention

    # convert eigenvalues of K into singular values of the SNP matrix S
    # Eigenvalues of K are non-negative because K is a positive semi-definite square matrix.
    # We are not interested in naught eigenvalues as they do not capture any variance in data.
    ev.pve <- which(ev > 0)
    num.ev0 <- length(ev) - length(ev.pve)  # number of zero eigenvalues
    if (num.ev0 > 0) {
        print(paste("Remove", num.ev0, "eigenvectors, which are related to naught eigenvalues.", sep = " "))
    }
    ev <- ev[ev.pve]  # only retain positive eigenvalues
    sv <- round(sqrt(L * ev), digits = 10)  # singular values of L * K = S * t(S)
    D <- diag(sv)  # the diagonal matrix of positive singular values used in the prcomp function (S = U * D * t(V))

    # read eigenvectors of K
    # Matrix L * K shares the same eigenvectors as K, but differs in eigenvalues and singular values.
    U <- as.matrix(read.table(file = out[["U"]], header = FALSE, sep = "\t"))  # unit eigenvectors of both K and L * K = S * t(S)
    U <- U[, -ncol(U)]  # remove the redundant NA column (the last column)
    U <- U[, ev.ord]  # rearrange eigenvectors following eigenvalues
    U <- U[, ev.pve]  # eigenvectors corresponding to positive eigenvalues - capturing variance in the cgSNP data
    rownames(U) <- samples
    colnames(U) <- paste("e", 1 : ncol(U), sep = "")  # assign eigenvector IDs to columns

    # We do not need to calculat anther matrix V (of orthogonal columns) because C = S * V = U * D.
    C <- U %*% D  # projections of samples into a new space of ncol(U) dimensions; "C" stands for "casts".
    rownames(C) <- samples
    colnames(C) <- paste("c", 1 : ncol(C), sep = "")  # IDs of cast vectors - projections

    # get pairwise projection distances of samples
    if (get.dists | get.tree) {
        print("Calculating pairwise distances between projections of samples.")
        sample.ds <- dist(x = C, method = dist.method, diag = FALSE, upper = FALSE)
        if (get.tree) {
            require(ape)
            require(phytools)
            print("Constructing a midpoint-rooted neighbour-joining tree for samples using the distance matrix.")
            tr.proj <- midpoint.root(nj(sample.ds))  # the tree, midpoint rooted
        } else {  # Do not build the tree but return the distances.
            tr.proj <- NA
        }
        sample.ds <- as.matrix(sample.ds)  # convert a dist object into a matrix
    } else {
        sample.ds <- NA
        tr.proj <- NA
    }

    # return eigenvalues (ev) and eigenvectors (U) of K, singular values (sv) and projections (C) of S
    return(list(C = C, sv = sv, ev = ev, U = U, K = K, d = sample.ds, tr = tr.proj))
}
