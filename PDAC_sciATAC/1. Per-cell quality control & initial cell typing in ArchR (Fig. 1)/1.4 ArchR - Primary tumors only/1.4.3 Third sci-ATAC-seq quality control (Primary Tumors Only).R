library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)

PlotPrefix = "Primary_TSSMixed_nFrags1k_FRIS0.4"



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSS2.5_nFrags1k_FRIS0.4/AllCells")

# Filter for TSS>=4, but only for non-epithelial clusters
proj <- proj[(proj$Clusters %in% c("C1","C2","C3","C4","C5","C6","C7","C9","C10","C17","C18")) | (proj$TSSEnrichment>=4)]

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells", dropCells = T)



############################################################
# iLSI and UMAP clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",force=TRUE,iterations=3, clusterParams = list(resolution=2),varFeatures = 100000)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force=TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE,minDist = 0.1)
proj <- addImputeWeights(proj)



############################################################
# Export UMAP plots and outputs

p <- list()
p[["UMAP"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p[["AllSamples"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample_ID", embedding = "UMAP", labelMeans = F) + ggtitle("All Samples")
for(sample in sort(unique(proj$Sample_ID))){
    p[[sample]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", highlightCells = proj$cellNames[proj$Sample_ID==sample], labelMeans = F) + ggtitle(paste("Sample ",sample, sep = ""))}

plotPDF(p, name = paste(PlotPrefix, "_UMAP_Clusters_SampleID.pdf", sep = ""),
        ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

# Get UMAP Dimensions File
dims <- data.frame(Barcode = proj$cellNames,
                   Dimension_1 = proj@embeddings@listData[["UMAP"]]@listData[["df"]][[1]],
                   Dimension_2 = proj@embeddings@listData[["UMAP"]]@listData[["df"]][[2]])
dims <- data.frame(Barcode = gsub("^([^#])+.", "", dims$Barcode),
                   Dimension_1 = dims$Dimension_1,
                   Dimension_2 = dims$Dimension_2,
                   Sample = gsub(".[^#]*$", "", dims$Barcode))
write.table(dims, file = paste(getOutputDirectory(proj),"/Plots/", PlotPrefix, "_UMAPDimensions.dims", sep = ""), row.names = F, col.names = F, sep = "\t", quote = F)

# Get Cluster annotation file
clusterAnnot <- data.frame(Barcode = proj$cellNames,
                           Cluster = proj$Clusters)
clusterAnnot <- data.frame(Barcode = gsub("^([^#])+.", "", clusterAnnot$Barcode),
                           Cluster = proj$Clusters,
                           Sample = gsub(".[^#]*$", "", clusterAnnot$Barcode))
write.table(clusterAnnot, file = paste(getOutputDirectory(proj),"/Plots/", PlotPrefix, "_ClusterAnnotation.annot", sep = ""), row.names = F, col.names = F, sep = "\t", quote = F)



######################################################################
# Save ArchRProject

proj <- saveArchRProject(proj)


