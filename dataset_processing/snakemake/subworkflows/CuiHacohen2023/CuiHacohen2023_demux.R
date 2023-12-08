library('Seurat')

# Load the data
data <- Seurat::ReadMtx(snakemake@input[['sample_mtx']], snakemake@input[['sample_barcodes']], snakemake@input[['sample_features']], strip.suffix=T)
tags <- Seurat::ReadMtx(snakemake@input[['tags_mtx']], snakemake@input[['tags_barcodes']], snakemake@input[['tags_features']], feature.column=1)
# data = Seurat::ReadMtx('cytokine-samples17-matrix.mtx', 'cytokine-samples17-barcodes.tsv', 'cytokine-samples17-features.tsv', strip.suffix=T)
# tags = Seurat::ReadMtx('cytokine-hashtags17-matrix.mtx', 'cytokine-hashtags17-barcodes.tsv', 'cytokine-hashtags17-features.tsv', feature.column=1)

# Subset RNA and HTO counts by joint cell barcodes
joint.bcs <- intersect(colnames(data), colnames(tags))
print(length(joint.bcs))
data <- data[, joint.bcs]
tags <- tags[, joint.bcs]
print(rowSums(tags))
tags <- tags[rowSums(tags) > 10000,] # Remove tags with <10k hashtag counts

# Create Seurat object and demultiplex
seurat_object <- Seurat::CreateSeuratObject(counts = data)
seurat_object[["HTO"]] <- Seurat::CreateAssayObject(counts = tags)
seurat_object <- Seurat::NormalizeData(seurat_object, assay = "HTO", normalization.method = "CLR")
seurat_object <- Seurat::MULTIseqDemux(seurat_object, assay = "HTO")
print(head(seurat_object[[]]))

# Save the demultiplexed metadata
write.csv(seurat_object[[]], snakemake@output[['demux']])