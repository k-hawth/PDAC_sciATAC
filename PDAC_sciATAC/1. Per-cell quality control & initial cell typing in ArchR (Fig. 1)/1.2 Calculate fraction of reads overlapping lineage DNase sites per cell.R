library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(32)
set.seed(1)



############################################################
# Prepare new ArchRProject

# Load root Arrow File containing all filtered and unfiltered cells
proj <- ArchRProject(
  ArrowFiles = "~/ArchR/bam_q30/ArrowFile_Root/sears_tumors_hg38.arrow", 
  outputDirectory = "~/ArchR/QualityMetricsTest",
  copyArrows = TRUE
)


# Read in sample annotation
annot <- read.table("~/sample_id.annot")
annot$V1 <- paste("sears_tumors_hg38#", annot$V1, sep = "")
annot <- annot[annot$V1 %in% proj$cellNames,]
proj <- addCellColData(proj, data = annot$V2, cells = annot$V1, name = "Sample_ID")



# Add peak set - DHS index
# Peaks obtained from Meuleman W, et al. (2020) (https://doi.org/10.1038/s41586-020-2559-3))
peaks <- read.table("~/peak_sets/DHS_Index_Mueleman_2020/DHS_Index_and_Vocabulary_hg38_WM20190703.non_overlap.bed")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))

proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)



# Remove cells with under 500 reads
proj <- proj[proj$nFrags>=500]



# Record quality metrics
df <- data.frame(cell = proj$cellNames, Sample_ID = proj$Sample_ID, TSSEnrichment = proj$TSSEnrichment, nFrags = proj$nFrags, PromoterRatio = proj$PromoterRatio, FRIS_IndexDHS = proj$FRIP)
write.table(df, file = "~/ArchR/AllSamples_nFrags500_QualityMetrics.txt", row.names = F, col.names = T, sep = '\t', quote = F)