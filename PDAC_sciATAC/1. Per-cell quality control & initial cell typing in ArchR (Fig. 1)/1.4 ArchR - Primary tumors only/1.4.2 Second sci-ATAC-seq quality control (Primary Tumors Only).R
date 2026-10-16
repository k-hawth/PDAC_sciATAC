library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)

PlotPrefix = "Primary_TSS2.5_nFrags1k_FRIS0.4"



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/bam_q30/PrimarySamples_TSS2.5_nFrags1k/AllCells")

# Read in per-cell fraction of reads overlapping Meuleman W, et al. (2020) DHS sites
qualityMetrics <- read.table("~/ArchR/AllSamples_nFrags500_QualityMetrics.txt",header = T, comment.char="")

# Remove any cell with less than 40% fraction of reads overlapping lineage DHS sites
proj <- proj[proj$cellNames %in% qualityMetrics$cell[qualityMetrics$FRIS_IndexDHS>=0.4]]

# Save filtered ArchRProject
proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSS2.5_nFrags1k_FRIS0.4/AllCells", dropCells = T)



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



############################################################
# Prepare new ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSS2.5_nFrags1k_FRIS0.4/AllCells")

proj <- saveArchRProject(proj, outputDirectory = "~/ArchR/Primary_TSS2.5_nFrags1k_FRIS0.4/AllCells_IndexDHS")



############################################################
# Add Index DHS peak set to ArchRProject

# Add peak set - DHS index
# Peaks obtained from Meuleman W, et al. (2020) (https://doi.org/10.1038/s41586-020-2559-3))
peaks <- read.table("~/peak_sets/DHS_Index_Meuleman_2020/DHS_Index_and_Vocabulary_hg38_WM20190703.non_overlap.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)

# Annotate pan-human peak set by component
Index_DHS_components <- c(
    Tissueinvariant = '~/peak_sets/DHS_Index_Meuleman_2020/Tissueinvariant.txt',
    Placentaltrophoblast = '~/peak_sets/DHS_Index_Meuleman_2020/Placentaltrophoblast.txt',
    Neural = '~/peak_sets/DHS_Index_Meuleman_2020/Neural.txt',
    Primitiveembryonic = '~/peak_sets/DHS_Index_Meuleman_2020/Primitiveembryonic.txt',
    Digestive = '~/peak_sets/DHS_Index_Meuleman_2020/Digestive.txt',
    Lymphoid = '~/peak_sets/DHS_Index_Meuleman_2020/Lymphoid.txt',
    Musculoskeletal = '~/peak_sets/DHS_Index_Meuleman_2020/Musculoskeletal.txt',
    StromalB = '~/peak_sets/DHS_Index_Meuleman_2020/StromalB.txt',
    Myeloiderythroid = '~/peak_sets/DHS_Index_Meuleman_2020/Myeloiderythroid.txt',
    Cancerepithelial = '~/peak_sets/DHS_Index_Meuleman_2020/Cancerepithelial.txt',
    Pulmonarydevel = '~/peak_sets/DHS_Index_Meuleman_2020/Pulmonarydevel.txt',
    Renalcancer = '~/peak_sets/DHS_Index_Meuleman_2020/Renalcancer.txt',
    StromalA = '~/peak_sets/DHS_Index_Meuleman_2020/StromalA.txt',
    Organdevel.renal = '~/peak_sets/DHS_Index_Meuleman_2020/Organdevel.renal.txt',
    Vascularendothelial = '~/peak_sets/DHS_Index_Meuleman_2020/Vascularendothelial.txt',
    Cardiac = '~/peak_sets/DHS_Index_Meuleman_2020/Cardiac.txt'
)

proj <- addPeakAnnotations(ArchRProj = proj, regions = Index_DHS_components, name = "Index_DHS",force=T)

proj <- addBgdPeaks(proj,method = "ArchR",force=T)

proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = "Index_DHS",
    force = TRUE
)



############################################################
# Plot DHS Index Components

markerDHS <- getFeatures(proj, useMatrix = "Index_DHSMatrix")
markerDHS <- sort(grep("z:", markerDHS, value = TRUE))


proj <- addImputeWeights(proj)
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "Index_DHSMatrix", 
    name = markerDHS, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)



p <- plotGroups(proj,
                groupBy = "Clusters",
                colorBy = "Index_DHSMatrix",
                name = markerDHS,
                imputeWeights = getImputeWeights(proj),
                maxCells = 10000000,
                quantCut = c(0,1),
                plotAs = "violin"
)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)