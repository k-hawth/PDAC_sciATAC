library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/AllSamples_TSSMixed_nFrags1k_FRIS0.4")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3")



############################################################
# Apply filters

# Filter to include epithelial/malignant clusters only
proj <- proj[proj$Clusters %in% c("C2","C3","C4","C5","C6","C7","C8","C9","C11","C12")==F]



############################################################
# Save filtered ArchRProject

proj <- saveArchRProject(proj, dropCells = T)



############################################################
# iLSI and UMAP clustering

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",force=TRUE,iterations=3, clusterParams = list(resolution=2),varFeatures = 100000)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force=TRUE)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=TRUE,minDist = 0.2)
proj <- addImputeWeights(proj)



############################################################
# Add peak matrix

peaks <- read.table("~/PeakCalling/MACS3/AllSamples_TSSMixed_nFrags1k_FRIS0.4/AllSamplesTSSMixednFrags1kFRIS0.4_MACS3.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)

proj <- addPeakMatrix(proj)



############################################################
# Run chromVAR

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



######################################################################
# Save ArchRProject

proj <- saveArchRProject(proj)

# Save new ArchRProject that will have cisTopic-based UMAP and PhenoGraph clusters added
proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# Add cisTopic data

# UMAP
cisTopic_dims <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_ScaledMatrix_UMAP_n30_minDist0.1.dims")

# UMAP name in ArchR
umap_name = "cisTopic_50"

# PhenoGraph clusters
cisTopic_annot <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Scaled.k15.pg.annot")

# pg annotation name in ArchR
pg_name = "cistopic50_pg_k15"



# Add cisTopic-based UMAP to ArchRProject
proj <- addUMAP(proj, name = umap_name, saveModel = F,force = T)
newdims <- data.frame(cell = proj$cellNames)
cisTopic_dims$V1 <- paste("sears_tumors_hg38#",cisTopic_dims$V1,sep="")
newdims$`IterativeLSI#UMAP_Dimension_1` <- cisTopic_dims$V2[match(newdims$cell, cisTopic_dims$V1)]
newdims$`IterativeLSI#UMAP_Dimension_2` <- cisTopic_dims$V3[match(newdims$cell, cisTopic_dims$V1)]
rownames(newdims) <- newdims$cell
newdims <- newdims[c(2,3)]
proj@embeddings@listData[[umap_name]]@listData[["df"]] <- newdims
rm(cisTopic_dims,newdims)

# Add cisTopic-based PhenoGraph clusters to  ArchRProject
cisTopic_annot$V1 <- paste("sears_tumors_hg38#",cisTopic_annot$V1, sep = "")
proj <- addCellColData(proj, cells = cisTopic_annot$V1, data = cisTopic_annot$V2, name = pg_name)



############################################################
# Plot cisTopic UMAP and clusters

p <- list()
p[["UMAP"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = pg_name, embedding = umap_name)
p[["AllSamples"]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample_ID", embedding = umap_name, labelMeans = F) + ggtitle("All Samples")
for(sample in sort(unique(proj$Sample_ID))){
    p[[sample]] <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = pg_name, embedding = umap_name, highlightCells = proj$cellNames[proj$Sample_ID==sample], labelMeans = F) + ggtitle(paste("Sample ",sample, sep = ""))}

plotPDF(p, name = paste(PlotPrefix, "_UMAP_Clusters_SampleID.pdf", sep = ""),
        ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)