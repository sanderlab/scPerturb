#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))

# CLI arguments
option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL,
                help="rds file path", metavar="character"),
    make_option(c("-n", "--name"), type="character", default="seurat",
                help="prefix for exported file names [default=%default]", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="./",
                help="output folder path [default=%default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (rds file path).n", call.=FALSE)
}
file <- opt$file
prefix <- opt$name
out <- opt$out

print('Reading metadata...')
seurat_data <- readRDS(file = file)
pca <- Seurat::Embeddings(seurat_data, reduction = "pca")
umap <- Seurat::Embeddings(seurat_data, reduction = "umap")
metadata <- seurat_data[[]]

# GEX
print('Reading GEX data...')
mat = Seurat::GetAssayData(object = seurat_data, slot = "counts")

# Export
print('Writing...')
Matrix::writeMM(mat, paste0(out, "/", prefix, '_matrix.mtx'))  # then could gzip it
write.csv(rownames(seurat_data), paste0(out, "/", prefix, '_genes.csv'))
write.csv(colnames(seurat_data), paste0(out, "/", prefix, '_barcodes.csv'))
write.csv(pca, paste0(out, "/", prefix, "_pca.csv"))
write.csv(umap, paste0(out, "/", prefix, "_umap.csv"))
write.csv(metadata, paste0(out, "/", prefix, "_metadata.csv"))

print('Done!')