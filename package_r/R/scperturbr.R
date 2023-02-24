# Title: scperturbR (E-statistics for Seurat objects)
# Author: Stefan Peidli
# Maintainer: Stefan Peidli <stefan.peidli at gmail.com>
# Description: Computes E-test statistics for each group in a Seurat object,
# using the E-distance in any space given (e.g. pca).
# Depends: R (>= 3.5.0), Seurat, dplyr, rdist, energy
# Date: 2023-02-22
# Publication: https://doi.org/10.1101/2022.08.20.504663
# URL: https://github.com/sanderlab/scPerturb
# BugReports: https://github.com/sanderlab/scPerturb/issues
# Encoding: UTF-8
# References: Rizzo & Székely (2013) https://doi.org/10.1016/j.jspi.2013.03.018

library(Seurat)
library(dplyr)
library(rdist)
library(energy)


#' @title etest
#'
#' @description Performs E-testing on a Seurat object.
#'     Computes E-test statistics for each group in a Seurat object,
#'     using the E-distance in space given by reduction to the group defined
#'     by control.
#' @param seurat_object An object of class Seurat.
#' @param groupy An object of class character. Points to the column in the
#'     Seurat object's meta data that contains the group labels.
#' @param control An object of class character. The group that is used as the
#'     control.
#' @param reduction An object of class character. The reduction / embedding in
#'     seurat_object that is used to compute the E-distance in.
#' @return Returns an object of class data.frame. For each group contains the
#'     E-test p-value and the E-distance to control group.
#' @examples
#'     # Add some code illustrating how to use the function
#' @importFrom Seurat Embeddings VariableFeatures RunPCA
#' @importFrom energy eqdist.etest
#' @importFrom stats dist
#' @importFrom dplyr select
#' @export
etest <- function(seurat_object, groupby = "perturbation", control = "control",
                  reduction = "pca", verbose = TRUE, permutations = 1000) {
    if (class(seurat_object) != "Seurat") {
        stop("The first argument must be a Seurat object.")
    }
    if (!(reduction %in% names(seurat_object@reductions))) {
        if (reduction == "pca") {
            if (verbose) {
                print("No PCA found in the Seurat object. Computing PCA now.")
            }
            var_features <- Seurat::VariableFeatures(seurat_object)
            seurat_object <- Seurat::RunPCA(seurat_object,
                                            features = var_features,
                                            verbose = FALSE)
        } else {
            stop("The specified reduction was not found in the Seurat object.")
        }
    }
    labels <- seurat_object[[]][[groupby]]
    groups <- unique(labels)
    if (!(control %in% groups)) {
        stop("The specified control group was not found in the groupby column 
             in the seurat_object metadata.")
    }
    emb <- Seurat::Embeddings(seurat_object, reduction = reduction)

    df <- data.frame(row.names = groups, pval = rep(NA, length(groups)))
    if (verbose) {
        print("Computing E-test statistics for each group.")
        }
    for (group in groups) {
        x <- as.matrix(emb)[labels == control, ]
        y <- as.matrix(emb)[labels == group, ]
        d <- stats::dist(rbind(x, y))
        res <- energy::eqdist.etest(d, sizes = c(nrow(x), nrow(y)),
                                    distance = TRUE, R = permutations)
        df[group, "pval"] <- res$"p.value"
        df[group, "edist"] <- res$statistic
    }
    return(df)
}

#' @title edist
#'
#' @description Computes pairwise E-distances on a Seurat object.
#'     Computes E-distance between all groups in a Seurat object in space
#'     given by reduction.
#' @param seurat_object An object of class Seurat.
#' @param groupy An object of class character. Points to the column in the
#'     Seurat object's meta data that contains the group labels.
#' @param reduction An object of class character. The reduction / embedding in
#'     seurat_object that is used to compute the E-distance in.
#' @param sample_correction An object of class logical. If TRUE, the
#'     E-distances are corrected for sample size. Will make it not a proper
#'     distance, leads to negative values.
#' @return Returns an object of class data.frame. For each group contains the
#'     E-test p-value and the E-distance to control group.
#' @examples
#'     # Add some code illustrating how to use the function
#' @importFrom Seurat Embeddings VariableFeatures RunPCA
#' @importFrom rdist cdist pdist
#' @importFrom dplyr select
#' @export
edist <- function(seurat_object, groupby = "perturbation", reduction = "pca",
                  sample_correction = FALSE, verbose = TRUE) {
    if (class(seurat_object)!="Seurat") {
        stop("The first argument must be a Seurat object.")
    }
    if (!(reduction %in% names(seurat_object@reductions))) {
        if (reduction == "pca") {
            if (verbose) {
                print("No PCA found in the Seurat object. Computing PCA now.")
            }
            var_features <- Seurat::VariableFeatures(seurat_object)
            seurat_object <- Seurat::RunPCA(seurat_object,
                                            features = var_features,
                                            verbose = FALSE)
        } else {
            stop("The specified reduction was not found in the Seurat object.")
        }
    }
    labels <- seurat_object[[]][[groupby]]
    groups <- unique(labels)
    emb <- Seurat::Embeddings(seurat_object, reduction = reduction)

    df <- setNames(data.frame(matrix(ncol = length(groups),
                                     nrow = length(groups)),
                              row.names = groups),
                   groups)
    if (verbose) {
        print("Computing E-test statistics for each group.")
    }
    completed_groups <- c()
    for (groupx in groups) {
        for (groupy in groups){
            if (groupy %in% completed_groups) {
                next
            }  # skip if already computed
            x <- as.matrix(emb)[labels == groupx, ]
            y <- as.matrix(emb)[labels == groupy, ]
            # this is the original edist function by Rizzo et al.
            # res <- energy::edist(c(x,y), s=c(50,50), distance = FALSE)  

            N <- nrow(x)
            M <- nrow(y)

            dist_xy <- rdist::cdist(x,y)
            dist_x <- rdist::pdist(x)
            dist_y <- rdist::pdist(y)

            if (sample_correction) {
                ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * (N - 1))) - (sum(dist_y) / (M * (M - 1)))
            } else {
                ed <- 2 * (sum(dist_xy) / (N * M)) - (sum(dist_x) / (N * N)) - (sum(dist_y) / (M * M))
            }

            df[groupx, groupy] <- df[groupy, groupx] <- ed
        }
        completed_groups <- c(completed_groups, groupx)
    }
    return(df)
}