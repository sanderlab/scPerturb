# This script is in parts taken from 
# https://github.com/GreenleafLab/SpearATAC_MS_2021/blob/main/AnalyzeSpearATAC/Scripts/Analysis-of-Spear-ATAC-Screen.R
# which was published under the MIT licence:
 
# MIT License

# Copyright (c) 2020 GreenleafLab

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



#todo:
# change sgsgNT to sgNT
# Use mat2triplet instead of summary for sparse matrices. Even though: 
# "The summary() method for "sparseMatrix", summary,sparseMatrix-method."" 
# Double check that i is row and j is col index. The R documentation has an error.



library(ArchR)
library(Biostrings)
library(here)
library(tidyverse)
library(rhdf5)
library(SummarizedExperiment)

#Helpful functions designed for SpearATAC Analysis
source(paste0(here(),"/scATAC/scATAC_utility_functions.R")) #R is beautiful
#Functions to create plots
source(paste0(here(),"/scATAC/plotting_functions.R"))

# dataset <- "Spear_ATAC"
# cell_line <-  "GM12878" #for spearATAC "K562", "MCF7", "GM12878"

# dataset <- "CRISPR_sciATAC"
# cell_line <-  "K562_2" #for CRISPR_sciATAC "K562_1", "K562_2"

dataset <- "ASAP_seq"
cell_line <-  "CD4+_T_cells" #for ASAP_seq "CD4+_T_cells"


#QC paras
minTSS = 4 #TSS enrichment score (signal to noise ratio)
minFrags = 1000 #min number mapped frags per cell
maxFrags = 1e+05 #max number mapped frags per cell
nSg <- 20 #min number of sgRNA counts for max Assignment
Spec <- 0.8 #min fraction of number max Assignment / total sgRNA count
# min_PurityRatio <- 0.9 #Spear-ATAC method to disregard cells whose neighborhood does not have mostly the same sgRNAs

output_dir = paste0(here(), "/scATAC/output/", dataset, "/", cell_line)
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

if(dataset == "Spear_ATAC"){
  addArchRGenome("hg38")
} 
if(dataset == "CRISPR_sciATAC"){
  addArchRGenome("hg19")
}
if(dataset == "ASAP_seq"){
  addArchRGenome("hg38")
}
addArchRThreads(100)


#Make Sure Data is present
#source("Scripts/Download-Test-SpearATAC-Data.R")

#Input Fragments
fragmentFiles <- getFileNames(dataset, cell_line, "fragments")
sgRNAFiles <- getFileNames(dataset, cell_line, "sgRNA_assignment") #sgRNA Assignments Files whose names match the sample names in the ArchRProject

if(dataset == "Spear_ATAC"){
  barcodeFiles <- getFileNames(dataset, cell_line, "barcodes")
  #Get Valid Barcodes
  validBC <- getValidBarcodes(
      csvFiles = unname(barcodeFiles),
      sampleNames = names(barcodeFiles)
    )
  } else if( (dataset == "CRISPR_sciATAC") | (dataset == "ASAP_seq")){
    # 10x barcode output file not provided so I cannot run this check and assume all bcs to be valid.
    validBC = NULL
  }


#Create Arrow Files
setwd(output_dir) #because there doesn't seem to be a way to give createArrowFiles an output folder for all files
ArrowFiles <- createArrowFiles(
	inputFiles=fragmentFiles, 
	# validBarcodes = validBC,
	minTSS = minTSS,
  minFrags = minFrags,
  maxFrags = maxFrags,
  )
#setwd(here()) #I better go back

#Make an ArchR Project
proj <- ArchRProject(ArrowFiles)

#Create an sgRNA assignment matrix
#Iterate over each sgRNA aligned file
sgAssign <- assign_sgRNA(sgRNAFiles, proj, dataset)

#str_replace(sgAssign$sgAssign, "sgsgNT", "control")
if(dataset == "CRISPR_sciATAC"){
  #adapt to the nomenclature of Spear-Atac so that the rest of the scripts can just run
  if( cell_line == 'K562_1'){
    sgAssign$sgAssign <- str_replace(sgAssign$sgAssign, "non_targeting", "sgsgNT")
  }
  if( cell_line == 'K562_2'){
    sgAssign$sgAssign <- str_replace(sgAssign$sgAssign, "NT", "sgsgNT")
  }
}
if(dataset == "ASAP_seq"){
  #adapt to the nomenclature of Spear-Atac so that the rest of the scripts can just run
  sgAssign$sgAssign <- str_replace(sgAssign$sgAssign, "sgNTC", "sgsgNT")
}
#Add all this information to your ArchRProject
proj <- addCellColData(proj, data = sgAssign$sgAssign, name = "sgAssign", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgCounts, name = "sgCounts", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgTotal, name = "sgTotal", cells = rownames(sgAssign), force = TRUE)
proj <- addCellColData(proj, data = sgAssign$sgSpec, name = "sgSpec", cells = rownames(sgAssign), force = TRUE)

#Identify those sgRNA Assignemnets that passed your cutoffs
proj$sgAssign2 <- NA
proj$sgAssign2[which(proj$sgCounts > nSg & proj$sgSpec > Spec)] <- proj$sgAssign[which(proj$sgCounts > nSg & proj$sgSpec > Spec)]
if(dataset == "Spear_ATAC"){
  proj$sgAssign3 <- stringr::str_split(proj$sgAssign2,pattern="\\-",simplify=TRUE)[,1]
}
if( (dataset == "CRISPR_sciATAC") | (dataset == "ASAP_seq")){
  proj$sgAssign3 <- stringr::str_split(proj$sgAssign2,pattern="\\_",simplify=TRUE)[,1]
}

#Print Numbers
#table(proj$sgAssign3)
# saveRDS(proj, paste0(output_dir, "/ArchRProject-W-sgRNA-Assignments.rds"))

#We now want to clean our sgAssignments based on the homogeneity of the sgRNA signal. This analysis shouldnt filter more than 5-10%
#of sgRNA assignments. This will help resolve your differential analyses but is not crucial for downstream analysis.

# proj <- cleanSgRNA(ArchRProj = proj, groupSg = "sgAssign3", individualSg = "sgAssign2", nonTarget = "sgsgNT")

#Clean Up sgRNA assignments based on Purity Ratio (NOTE you may need to adjust this cutoff)
# proj$sgAssignFinal <- "UNK"
# proj$sgAssignFinal[proj$PurityRatio >= min_PurityRatio] <- proj$sgAssign3[proj$PurityRatio >= min_PurityRatio]
# proj$sgIndividual <- ifelse(proj$sgAssignFinal=="UNK", "UNK", proj$sgAssign2)


#I do not want to run the purity analysis, so I need to do this:
proj$sgAssignFinal<- proj$sgAssign3
proj$sgIndividual <- proj$sgAssign2

#Lets create an unbiased LSI-UMAP to sgRNA assignments by creating an iterativeLSI reduction + UMAP
proj <- addIterativeLSI(proj)
proj <- addUMAP(proj)


#Call Peaks using sgAssignments
proj <- addGroupCoverages(proj, groupBy = "sgAssignFinal", force = TRUE)
proj <- addReproduciblePeakSet(proj, 
	groupBy = "sgAssignFinal", force = TRUE, 
)

#Create PeakMatrix from Peak Set
proj <- addPeakMatrix(proj, force = TRUE)

#Create Motif Deviations Matrix using Vierstra Motifs
motifPWMs <- readRDS(paste0(here(),"/scATAC/data/Vierstra-Human-Motifs.rds"))

#these steps take an hour:
proj <- addMotifAnnotations(proj, motifPWMs = motifPWMs, name = "Vierstra")
proj <- addDeviationsMatrix(proj, peakAnnotation = "Vierstra", force = TRUE, matrixName = 'ChromVar')

saveRDS(proj, paste0(output_dir, "/ArchRProject.rds"))
# proj <- readRDS(paste0(output_dir, "/ArchRProject.rds"))

####################################
# Compute Differential Peaks
####################################

bgd <- grep("sgNT", proj$sgIndividual, value=TRUE) %>% unique
# bgd <- bgd[!grepl("-11|-12|-8|-5|-6", bgd)]
proj$sgAssignClean <- proj$sgAssignFinal
# proj$sgAssignClean[grepl("-11|-12|-8|-5|-6", proj$sgIndividual)] <- "UNK"

#Sort sgRNA Targets so Results are in alphabetical order
useGroups <- sort(unique(proj$sgAssignClean)[unique(proj$sgAssignClean) %ni% c("UNK")])

if((dataset == "CRISPR_sciATAC") & (cell_line == 'K562_2')){
  useGroups <- useGroups[!(useGroups == "SUPT16H")] #no cells in that group
}

bgdGroups <- unique(grep("sgNT", proj$sgAssignClean,value=TRUE))

#Differential Peaks
diffPeaks <- getMarkerFeatures(
	ArchRProj = proj, 
	testMethod = "binomial",
	binarize = TRUE,
	useMatrix = "PeakMatrix",
	useGroups = useGroups, 
	bgdGroups = bgdGroups, 
	groupBy = "sgAssignClean",
	bufferRatio = 0.95,
	maxCells = 250,
)


###plotting
#I need to test if this works here, or if I have to plot earlier.
# plot_QC(proj, sgAssign, nSg, sgSpec)


#can't find marker genes

#  diff_genes <- getMarkerFeatures(
# 	ArchRProj = proj, 
#     testMethod = "wilcoxon",
#     binarize = FALSE,
#     useMatrix = "GeneScoreMatrix",
#     useGroups = useGroups, 
#     bgdGroups = bgdGroups, 
#     groupBy = "sgAssignClean",
#     # bufferRatio = 0.95,
#     # maxCells = 250,
#     )

# gene_markers = getMarkers(diff_genes)

# gene_marker_matrix_z <- plotMarkerHeatmap(
#   seMarker = diff_genes, 
#   cutOff = "FDR <= 0.01 & Log2FC >= 0.5", 
#   transpose = TRUE,
#   returnMatrix = TRUE
# )


#I get nans when I run this. Not sure why, so I will leave it
# diff_TF <- getMarkerFeatures(
# 	ArchRProj = proj, 
#     testMethod = "wilcoxon",
#     binarize = FALSE,
#     useMatrix = "VierstraMatrix",
#     useGroups = useGroups, 
#     bgdGroups = bgdGroups, 
#     groupBy = "sgAssignClean",
# 	useSeqnames = 'z',
#     bufferRatio = 0.95,
#     maxCells = 250,
#     )

# plotMarkerHeatmap(
#   seMarker = diff_TF, 
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.1", 
#   transpose = TRUE,
#   returnMatrix = FALSE
# )


#gather everything I want to export#
####################################

# getAvailableMatrices(proj)


LSI <- getReducedDims(
      ArchRProj = proj,
      reducedDims = "IterativeLSI",
      returnMatrix = TRUE,
      dimsToUse = NULL,
      scaleDims = NULL,
      corCutOff = 0.75
    )

gene_scores <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
)
rowData_gene_scores =rowData(gene_scores)
colData_gene_scores =colData(gene_scores)
gene_scores <- assay(gene_scores)

peak_bc <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
)
rowData_peak_bc =rowData(peak_bc)
colData_peak_bc =colData(peak_bc)
peak_bc <- assay(peak_bc)

marker_peak_bc <- plotMarkerHeatmap(
  seMarker = diffPeaks, 
  cutOff = "FDR <= 0.05 & Log2FC >= 2",  #these are the Spear-ATAC thresholds
  transpose = TRUE,
  returnMatrix = TRUE
)

ChromVar <- getMatrixFromProject(
  ArchRProj = proj,
  useMatrix = "ChromVar",
)
rowData_ChromVar =rowData(ChromVar)
colData_ChromVar =colData(ChromVar)
ChromVar <- assays(ChromVar)$z

cell_meta_data <- getCellColData(proj) #I need to check that nucleosome ratio + TSSEnrichment are in there even when I don't use the plotting function eg. plotTSSEnrichment.
PeakSet <- getPeakSet(proj) #this is needed to annotate the peak_bc matrix

class(LSI)
class(gene_scores)
class(peak_bc)
class(marker_peak_bc)
class(ChromVar)
class(cell_meta_data)
class(PeakSet)

#export
export_dir <- paste0(output_dir, "/export/")
dir.create(export_dir, showWarnings = TRUE, recursive = TRUE)

write.csv(LSI, paste0(export_dir, "LSI.csv"))

write.csv(rowData_gene_scores, paste0(export_dir, "gene_scores_rowData.csv"))
write.csv(colData_gene_scores, paste0(export_dir, "gene_scores_colData.csv"))
write.csv(summary(gene_scores), file = gzfile(paste0(export_dir, "gene_scores.csv.gz")))

write.csv(rowData_peak_bc, paste0(export_dir, "peak_bc_rowData.csv"))
write.csv(colData_peak_bc, paste0(export_dir, "peak_bc_colData.csv"))
write.csv(summary(peak_bc), file = gzfile(paste0(export_dir, "peak_bc.csv.gz")))

write.csv(marker_peak_bc, paste0(export_dir, "marker_peak_bc.csv"))

write.csv(rowData_ChromVar, paste0(export_dir, "ChromVar_rowData.csv"))
write.csv(colData_ChromVar, paste0(export_dir, "ChromVar_colData.csv"))
write.csv(summary(ChromVar), file = gzfile(paste0(export_dir, "ChromVar.csv.gz")))

write.csv(cell_meta_data, paste0(export_dir, "cell_meta_data.csv"))
write.csv(PeakSet, paste0(export_dir, "PeakSet.csv"))

# matrix <- assays(diff_genes)[[2]][1:100,]
# matrix_long <- matrix %>% 
# 	as_tibble(rownames = 'gene') %>% 
# 	mutate(gene = as.integer(gene)) %>% 
# 	pivot_longer(cols = !gene, names_to = "guide")

# gg <- matrix_long %>% 
# 	ggplot(aes(x=gene, y=guide, fill=value)) + 
# 	geom_tile()

# plotPDF(gg, name = "test", width = 6, height = 6, addDOC=FALSE)
