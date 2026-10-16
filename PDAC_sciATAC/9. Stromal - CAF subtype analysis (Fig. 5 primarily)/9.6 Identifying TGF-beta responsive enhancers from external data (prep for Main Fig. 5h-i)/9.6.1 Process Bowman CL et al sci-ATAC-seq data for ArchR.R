########################################################################################################################
# Analysis of external sci-ATAC-seq dataset from Bowman CL and Daniel CJ, et al. (2025) (https://doi.org/10.1158/0008-5472.can-24-3437)

# This section creates an ArrowFile and ArchRProject using the data from Bowman CL and Daniel CJ, et al. (2025)
# This external sci-ATAC-seq dataset profiles the chromatin response of pancreatic stellate cells (PSCs) to TGF-beta

# In these data, there is a shPIN1 condition (PIN1 knockout) and a shSCR condition (control).
# Here, I compare only the shSCR control conditions (shSCR untreated PSCs vs. shSCR TGFb-treated PSCs)

########################################################################################################################

library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)




ArrowFiles <- createArrowFiles( inputFiles ="~/PIN1_ArchR/HPSC_PIN1.bbrd.q10.filt.bam",
                                outputNames = "archr_pin1",
                                sampleNames = "archr_pin1",
                                minTSS = 0,
                                minFrags = 0,
                                addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                maxFrags = 100000000000,
                                gsubExpression= ":.*")

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10,
    knnMethod = "UMAP",
    LSIMethod = 1
)


proj <- ArchRProject(
    ArrowFiles = "~/PIN1_ArchR/archr_pin1.arrow", 
    outputDirectory = "~/PIN1_ArchR/Project_PDAC_sciATAC",
    copyArrows = TRUE
)



############################################################
# Load the annotation file designating the treatment group of each cell

condition <- read.table("~/PIN1_ArchR/Dims and Annots/celltype_cond.annot")
condition$V1 <- paste("archr_pin1#",condition$V1,sep="")

annot <- data.frame(cell = proj$cellNames)
annot$group <- condition$V2[match(annot$cell,condition$V1)]

proj <- addCellColData(proj, data = annot$group, cells = annot$cell, name = "sample_ID")



############################################################
# Add background set of peaks 
# This is the set of peaks called in the original paper by Bowman CL and Daniel CJ, et al. (2025)

peaks <- read.table("~/peak_sets/LangerE_TGFB_PSC/background_called_peaks_BowmanCL_DanielCJ_2025.txt")
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))
proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)