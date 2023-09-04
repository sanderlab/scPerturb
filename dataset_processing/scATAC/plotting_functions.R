

plot_QC <- function(proj, sgAssign, nSg, sgSpec){
    #Plot QC of scATAC-seq data
    p1 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment", plotAs = "violin", alpha = 0.4)
    p2 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "PromoterRatio", plotAs = "violin", alpha = 0.4)
    p3 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "NucleosomeRatio", plotAs = "violin", alpha = 0.4)
    p4 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", name = "BlacklistRatio", plotAs = "violin", alpha = 0.4)
    plotPDF(p1, p2, p3, p4, name = "QC-singles-scATAC", addDOC = FALSE, width = 4, height = 4)

    #Plot Pseudo-Bulk Average QC of scATAC-seq data
    p1 <- plotTSSEnrichment(proj)
    p2 <- plotFragmentSizes(proj)
    plotPDF(p1, p2, name = "QC-pseudobulk-scATAC", addDOC = FALSE, width = 4, height = 4)

    #Plot Cutoffs of sgRNA Assignment (NOTE: You may want to adjust these cutoffs based on your results!)
    p <- ggPoint(log10(sgAssign$sgCounts), sgAssign$sgSpec, colorDensity = TRUE) +
    geom_vline(xintercept = log10(nSg), lty = "dashed") + 
    geom_hline(yintercept = Spec, lty = "dashed") +
    xlab("log10(sgCounts)") + ylab("Specificity")
    plotPDF(p, name = "Plot-Assignment-Density", addDOC = FALSE, width = 5, height = 5)
}

plot_UMAP <- function(proj){
    #Create Color Palettes
    pal1 <- paletteDiscrete(paste0(unique(proj$sgAssign2)))
    pal1["NA"] <- "lightgrey"

    pal2 <- paletteDiscrete(paste0(unique(proj$sgAssign3)))
    pal2["NA"] <- "lightgrey"

    pal3 <- paletteDiscrete(paste0(unique(proj$sgAssignFinal)))
    pal3["UNK"] <- "lightgrey"


    p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample")
    p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssign2", pal = pal1)
    p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssign3", pal = pal2)
    p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "sgAssignFinal", pal = pal3)
    plotPDF(p1, p2, p3, p4, name = "Plot-UMAP", addDOC=FALSE)
}
