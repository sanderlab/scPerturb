path <- "inst/"
expression_matrix <- Seurat::ReadMtx(
  paste0(path, 'matrix.mtx'),
  paste0(path, 'barcodes.csv'),
  paste0(path, 'genes.csv'),
  cell.column = 1,
  feature.column = 2,
  cell.sep = ",",
  feature.sep = ",",
  skip.cell = 1,
  skip.feature = 1,
)
latent <- read.csv(paste0(path, "/seurat_pca.csv"), row.names = 1)
metadata <- read.csv(paste0(path, "/seurat_metadata.csv"), row.names = 1)
colnames(expression_matrix) <- rownames(metadata)
seurat_object_reloaded <- CreateSeuratObject(counts = expression_matrix,
                                             meta.data = metadata)
seurat_object_reloaded[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(latent), key = "pca_", assay = DefaultAssay(seurat_object_reloaded))

print(getwd())
et <- etest(seurat_object = seurat_object_reloaded, groupby = 'perturbation', 
            control = 'control', reduction = 'pca')
