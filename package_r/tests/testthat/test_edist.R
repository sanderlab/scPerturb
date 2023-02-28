path <- "testdata"
expression_matrix <- Seurat::ReadMtx(
  testthat::test_path(path, "matrix.mtx"),
  testthat::test_path(path, "barcodes.csv"),
  testthat::test_path(path, "genes.csv"),
  cell.column = 1,
  feature.column = 2,
  cell.sep = ",",
  feature.sep = ",",
  skip.cell = 1,
  skip.feature = 1,
)
latent <- read.csv(testthat::test_path(path, "seurat_pca.csv"), row.names = 1)
metadata <- read.csv(testthat::test_path(path, "seurat_metadata.csv"),
                     row.names = 1)
colnames(expression_matrix) <- rownames(metadata)
seurat_object <- Seurat::CreateSeuratObject(counts = expression_matrix,
                                             meta.data = metadata)
seurat_object[["pca"]] <- Seurat::CreateDimReducObject(embeddings = as.matrix(latent),
                                               key = "pca_",
                                               assay = Seurat::DefaultAssay(seurat_object))

df <- edist(seurat_object, groupby = "perturbation", reduction = "pca",
            sample_correction = FALSE, verbose = TRUE)
