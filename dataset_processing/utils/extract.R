#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2023-11-23
# @Author  : Stefan Peidli
# @Version : $Id$
# @Description: extract data from seurat object
# Takes a Seurat rds file and extracts the data into a format that can be used
# in other programming languages. Outputs the following files to the output
# folder:
# - matrix.mtx: sparse matrix in matrix market format
# - genes.csv: gene names
# - barcodes.csv: cell barcodes
# - seurat_EMB.csv: Embedding coordinates for any reduction in the data object
# - seurat_metadata.csv: metadata (or replaces "seurat" by "name" option)

library("optparse")
library(Seurat)
library(Matrix)

# parse command line arguments
option_list <- list(
  make_option(c("-f", "--file"), type = "character", default = NULL,
              help = "rds file path", metavar = "character"),
  make_option(c("-n", "--name"), type = "character", default = "seurat",
              help = "output file name prefix [default=%default]",
              metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "./",
              help = "output folder path [default= %default]",
              metavar = "character")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$file)) {
  print_help(opt_parser)
  stop("At least one argument must be supplied (rds file path).n",
    call. = FALSE
  )
}


print("Reading metadata...")
seurat_data <- readRDS(file = opt$file)
metadata <- seurat_data[[]]

print("Reading GEX data...")
mat <- GetAssayData(object = seurat_data, slot = "counts")

print("Writing...")
target <- paste0(opt$out, "/", opt$name)
Matrix::writeMM(mat, paste0(target, "_matrix.mtx"))  # then gzip it optionally
write.csv(rownames(seurat_data), paste0(target, "_genes.csv"))
write.csv(colnames(seurat_data), paste0(target, "_barcodes.csv"))
write.csv(metadata, paste0(target, "_metadata.csv"))
# Embeddings
for (embedding in names(seurat_data@reductions)) {
  df <- seurat_data@reductions[[embedding]][[]]
  write.csv(df, paste0(target, "_", embedding, ".csv"))
}

print("Finished extraction from rds file.")
