library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(50)
set.seed(1)



############################################################
# Create root unfiltered ArchR Arrow file from aligned sequencing data

ArrowFiles <- createArrowFiles( inputFiles ="sears_tumors_hg38na.bbrd.q30.bam",
                                sampleNames = "sears_tumors_hg38",
                                minTSS = 0, 
                                minFrags = 0,
                                maxFrags = 10000000000000000000000,
                                addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                gsubExpression= ":.*")

# moved arrow file



############################################################
# Create ArchRProject and filter for TSS>=2.5 and nFrags>=1000

proj <- ArchRProject(
    ArrowFiles = "~/ArchR/bam_q30/ArrowFile_Root/sears_tumors_hg38.arrow", 
    outputDirectory = "~/ArchR/bam_q30/DoubletCalling/AllSamples_TSS2.5_nFrags1k_WithDoublets",
    copyArrows = TRUE
)

# Filter cells that do not meet cutoffs
proj <- proj[proj$TSSEnrichment>=2.5 & proj$nFrags>=1000]

# Save ArchRProject
proj <- saveArchRProject(proj, dropCells = TRUE)



############################################################
# Add sample ID to each cell

# Read in sample annotation
annot <- read.table("~/sample_id.annot")
annot$V1 <- paste("sears_tumors_hg38#", annot$V1, sep = "")
annot <- annot[annot$V1 %in% proj$cellNames,]
proj <- addCellColData(proj, data = annot$V2, cells = annot$V1, name = "Sample_ID")

# Save ArchRProject after adding sample ID
proj <- saveArchRProject(proj)



############################################################
# Split ArchRProject up based on which samples were multiplexed together

# Read in library prep annot
annot <- read.table("~/sears_disaggregates_allPreps.annot")
annot$V1 <- paste("sears_tumors_hg38#", annot$V1, sep = "")
annot <- annot[annot$V1 %in% proj$cellNames,]



# Create ArchRProjects split by samples that were multiplexed together
prep1 <- saveArchRProject(proj[proj$cellNames %in% annot$V1[annot$V2=="sears_disaggregates_library1"]],
                         outputDirectory = "~/ArchR/bam_q30/DoubletCalling/TSS2.5_nFrags1k/Prep1",
                         dropCells = TRUE)

prep2 <- saveArchRProject(proj[proj$cellNames %in% annot$V1[annot$V2=="sears_disaggregates_library2"]],
                          outputDirectory = "~/ArchR/bam_q30/DoubletCalling/TSS2.5_nFrags1k/Prep2",
                          dropCells = TRUE)

prep3a <- saveArchRProject(proj[proj$cellNames %in% annot$V1[annot$V2=="sears_disaggregates_library3_group1"]],
                          outputDirectory = "~/ArchR/bam_q30/DoubletCalling/TSS2.5_nFrags1k/Prep3a",
                          dropCells = TRUE)

prep3b <- saveArchRProject(proj[proj$cellNames %in% annot$V1[annot$V2=="sears_disaggregates_library3_group2"]],
                           outputDirectory = "~/ArchR/bam_q30/DoubletCalling/TSS2.5_nFrags1k/Prep3b",
                           dropCells = TRUE)

prep3c <- saveArchRProject(proj[proj$cellNames %in% annot$V1[annot$V2=="sears_disaggregates_library3_group3"]],
                           outputDirectory = "~/ArchR/bam_q30/DoubletCalling/TSS2.5_nFrags1k/Prep3c",
                           dropCells = TRUE)



# Add doublet scores to each prep project
prep1 <- addDoubletScores(prep1) # UMAP Projection R^2 = 0.92132
prep2 <- addDoubletScores(prep2) # UMAP Projection R^2 = 0.98675
prep3a <- addDoubletScores(prep3a) # UMAP Projection R^2 = 0.99135
prep3b <- addDoubletScores(prep3b) # UMAP Projection R^2 = 0.98747
prep3c <- addDoubletScores(prep3c) # UMAP Projection R^2 = 0.9954



############################################################
# Check the amount of doublets filtered by default filterRatio

filterDoublets(prep1) # sears_tumors_hg38 : 80 of 2840 (2.8%)
filterDoublets(prep2) # sears_tumors_hg38 : 817 of 9042 (9%)
filterDoublets(prep3a) # sears_tumors_hg38 : 583 of 7636 (7.6%)
filterDoublets(prep3b) # sears_tumors_hg38 : 1193 of 10925 (10.9%)
filterDoublets(prep3c) # sears_tumors_hg38 : 832 of 9122 (9.1%)



############################################################
# Increase filterRatio to filter >=10% doublets

filterDoublets(prep1, filterRatio = 3.51)
filterDoublets(prep2, filterRatio = 1.11)
filterDoublets(prep3a, filterRatio = 1.31)
filterDoublets(prep3b, filterRatio = 0.912)
filterDoublets(prep3c, filterRatio = 1.10)

# Compile list of doublets
doublets <- data.frame(V1 = proj$cellNames, V2 = "Doublet")
doublets$V2[doublets$V1 %in% filterDoublets(prep1, filterRatio = 3.51)$cellNames] = "Singlet"
doublets$V2[doublets$V1 %in% filterDoublets(prep2, filterRatio = 1.11)$cellNames] = "Singlet"
doublets$V2[doublets$V1 %in% filterDoublets(prep3a, filterRatio = 1.31)$cellNames] = "Singlet"
doublets$V2[doublets$V1 %in% filterDoublets(prep3b, filterRatio = 0.912)$cellNames] = "Singlet"
doublets$V2[doublets$V1 %in% filterDoublets(prep3c, filterRatio = 1.10)$cellNames] = "Singlet"

doublets$V1 <- gsub("sears_tumors_hg38#", "", doublets$V1)

# Output doublet list
write.table(doublets, file = "~/ArchR/bam_q30/DoubletCalling/DoubletCalls_10Percent_AllSamplesTSS2.5nFrags1k_q30bam.txt", row.names = F, col.names = F, sep = '\t', quote = F)



############################################################
# Save ArchRProjects

prep1 <- saveArchRProject(prep1)
prep2 <- saveArchRProject(prep2)
prep3a <- saveArchRProject(prep3a)
prep3b <- saveArchRProject(prep3b)
prep3c <- saveArchRProject(prep3c)