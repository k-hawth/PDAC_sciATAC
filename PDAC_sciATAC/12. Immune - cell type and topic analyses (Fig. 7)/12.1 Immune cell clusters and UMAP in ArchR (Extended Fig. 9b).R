library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)

PlotPrefix = "ImmuneTSSMixnFrags1kFRIS0.4"



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells")

# Filter for fibroblast clusters
proj <- proj[proj$Clusters %in% c("C2","C3","C4","C5","C6")]

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune", dropCells = T)



############################################################
# iLSI and UMAP clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",force=TRUE,iterations=3, clusterParams = list(resolution=2),varFeatures = 50000)
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



############################################################
# Prepare new ArchRProject that will have peak set added

proj <- saveArchRProject(proj)

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune_MACS2ArchR")



############################################################
# Add background peaks

peaks <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells_MACS2ArchR/Plots/Peaks_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)

proj <- addPeakMatrix(proj)




############################################################
# Add chromVAR TF motif accessibility matrix

# Add BSgenome.Hsapiens.UCSC.hg38
library(BSgenome.Hsapiens.UCSC.hg38)

# Add cisbp motif annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)

# Call background peaks
proj <- addBgdPeaks(proj,method = "ArchR", force=T)

# Run chromVAR
proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Motif",
    force = TRUE
)




############################################################
# Prepare new ArchRProject that will have cisTopic-based UMAP and clustering

proj <- saveArchRProject(proj)

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Immune_MACS2ArchR_cisTopic43")


PlotPrefix = "ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_cisTopic43"



############################################################
# Add cisTopic data

# UMAP
cisTopic_dims <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/cisTopic_50/Perplexity_topics/cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Scaled.UMAP.dims")

# UMAP name in ArchR
umap_name = "cisTopic_43"



# PhenoGraph clusters
cisTopic_annot <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Immune/cisTopic_50/Perplexity_topics/cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Scaled.k15.pg.annot")

# pg annotation name in ArchR
pg_name = "cistopic43_pg_k15"

# pg to ArchRProject
cisTopic_annot$V1 <- paste("sears_tumors_hg38#",cisTopic_annot$V1, sep = "")
proj <- addCellColData(proj, cells = cisTopic_annot$V1, data = cisTopic_annot$V2, name = pg_name)



# UMAP to ArchRProject
proj <- addUMAP(proj, name = umap_name, saveModel = F,force = T)
newdims <- data.frame(cell = proj$cellNames)
cisTopic_dims$V1 <- paste("sears_tumors_hg38#",cisTopic_dims$V1,sep="")
newdims$`IterativeLSI#UMAP_Dimension_1` <- cisTopic_dims$V2[match(newdims$cell, cisTopic_dims$V1)]
newdims$`IterativeLSI#UMAP_Dimension_2` <- cisTopic_dims$V3[match(newdims$cell, cisTopic_dims$V1)]
rownames(newdims) <- newdims$cell
newdims <- newdims[c(2,3)]
proj@embeddings@listData[[umap_name]]@listData[["df"]] <- newdims
rm(cisTopic_dims,newdims)



############################################################
# Plot cisTopic UMAP and clusters

p <- list()
p[["UMAP"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = pg_name, embedding = umap_name)
p[["AllSamples"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample_ID", embedding = umap_name, labelMeans = F) + ggtitle("All Samples")
for(sample in sort(unique(proj$Sample_ID))){
    p[[sample]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = pg_name, embedding = umap_name, highlightCells = proj$cellNames[proj$Sample_ID==sample], labelMeans = F) + ggtitle(paste("Sample ",sample, sep = ""))}

plotPDF(p, name = paste(PlotPrefix, "_UMAP_Clusters_SampleID.pdf", sep = ""),
        ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)



######################################################################
# Save ArchRProject

proj <- saveArchRProject(proj)


