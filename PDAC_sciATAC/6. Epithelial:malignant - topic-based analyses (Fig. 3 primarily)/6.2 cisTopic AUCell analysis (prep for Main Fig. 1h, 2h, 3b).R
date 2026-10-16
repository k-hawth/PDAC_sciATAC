library(cisTopic)
library(Matrix)
library(ggplot2)
set.seed(1)



########################################################################################################################
# SET VARIABLES
########################################################################################################################

PlotPrefix <- "Epithelial_cisTopic50"
PlotSuffix <- "cisTopicAUC"

# Specify folder where outputs will go. Folder will be created if one does not exist
outputDirectory <- "~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/AUC_Epithelial_Signatures"


cisTopicObject <- readRDS("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Derivative_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3.cistopicObject_Binarized.rds")
cisTopic_modelNumber <- 50

peaks <- read.table("~/PeakCalling/MACS3/AllSamples_TSSMixed_nFrags1k_FRIS0.4/AllSamplesTSSMixednFrags1kFRIS0.4_MACS3.bed")

ChIPSeeker_max_distToTSS <- 50000



########################################################################################################################
# PREPARE GENE SET LIST
########################################################################################################################

signatures <- as.data.frame(readxl::read_xlsx("~/gene_sets/PDAC/PDAC_Gene_Signatures.xlsx", sheet = 2))

# Arrange gene sets into a named list
gene_sets <- list()
for(gene_set in colnames(signatures)){
    gene_sets[[gene_set]] <- signatures[[gene_set]][is.na(signatures[[gene_set]])==F]
    rm(gene_set)
}



########################################################################################################################
# GENERATE "PEAK SETS" FROM GENE SET LIST
########################################################################################################################

# This code converts gene signatures into pseudo-peak signatures
# It does so by collecting all peaks whose nearest gene within 50kb corresponds to a gene in a gene signature

# Annotate peaks with nearest gene
chipseek <- function(peakset,
                     chr_col = "V1",
                     start_col = "V2",
                     end_col = "V3",
                     tss_range = c(-3000,3000)){
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    peak <- peakset[c(chr_col,start_col,end_col)]
    peak <- GenomicRanges::GRanges(seqnames = peak[,chr_col], ranges = IRanges::IRanges(start = peak[,start_col], end = peak[,end_col]))
    peakAnno <- ChIPseeker::annotatePeak(peak, tssRegion=tss_range,TxDb=txdb, annoDb="org.Hs.eg.db")
    peakAnno <- as.data.frame(peakAnno)
    
    if(identical(paste(peakset[,chr_col],peakset[,start_col],peakset[,end_col],sep = ""), paste(peakAnno$seqnames,peakAnno$start,peakAnno$end,sep=""))){
        peakAnno <- peakAnno[colnames(peakAnno) %in% c("seqnames","start","end")==F]
        peakset <- cbind(peakset,peakAnno)
        colnames(peakset) <- make.unique(colnames(peakset))
        return(peakset)
    } else{
        print("ChIPSeeker Error: Unequal peakAnno and peakset")
    }
}

peaks <- chipseek(peaks)



# Create output folder (if it does not exist)
system(paste("mkdir -p",outputDirectory)) 


# List of "peak set" file paths
peak_set_paths <- list()


# Generate text files for each "peak set" gene set
for(gene_set in names(gene_sets)){
    gene_set_peaks <- peaks
    gene_set_peaks <- gene_set_peaks[gene_set_peaks$SYMBOL %in% gene_sets[[gene_set]],]
    gene_set_peaks <- gene_set_peaks[abs(gene_set_peaks$distanceToTSS) <= ChIPSeeker_max_distToTSS,]
    gene_set_peaks <- gene_set_peaks[c(1:3)]
    geneset_path <- paste(outputDirectory,"/",gene_set,".txt",sep = "")
    write.table(gene_set_peaks, file = geneset_path, row.names = F, col.names = F, sep = '\t', quote = F)
    peak_set_paths[[gene_set]] <- geneset_path
    rm(gene_set_peaks, geneset_path)
}



########################################################################################################################
# PREP CISTOPIC OBJECT
########################################################################################################################

# Ensure the cisTopic object is set to use the correct topic model (50 in this case)
cisTopicObject <- selectModel(cisTopicObject, select = cisTopic_modelNumber)
cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=F)

dev.off()
system("rm Rplots.pdf")



########################################################################################################################
# AUC ANALYSIS: per-topic scores
########################################################################################################################

# Extract all peak sets
PlotSuffixGeneSet <- ".*"
signatures <- peak_set_paths[grep(paste("^",PlotSuffixGeneSet,"_",sep=""),names(peak_set_paths))]


# Change plot output name
PlotSuffixGeneSet <- "cisTopicAUC_SubtypeSignatures"


# Output AUC scores
auc <- signaturesHeatmap(getSignaturesRegions(cisTopicObject, signatures = signatures, labels=names(signatures), minOverlap = 0.4))
auc <- auc@ht_list[["Normalised AUC score"]]@matrix
outputFile <- paste(outputDirectory,"/",PlotPrefix,"_",PlotSuffix,"_", PlotSuffixGeneSet,"_AUCScoresForHeatmap.txt",sep="")
write.table(auc, file = outputFile, row.names = T, col.names = T, sep = '\t', quote = F)



########################################################################################################################
# AUC ANALYSIS: per-cell scores
########################################################################################################################

# Change plot output name
PlotSuffixGeneSet <- "AllMalignantSignatures"

# AUC threshold: this number will be multiplied by the number of total peaks
AUC_MaxRank = 0.03  #default is 0.03

cisTopicObject <- getSignaturesRegions(cisTopicObject, signatures = signatures, labels=names(signatures), minOverlap = 0.4)

# Build predictive matrix
pred.matrix <- predictiveDistribution(cisTopicObject)

# AUC cell rankings
aucellRankings <- AUCell_buildRankings(pred.matrix, plot=FALSE, verbose=FALSE)

# Calculate signature enrichment (cell version)
cisTopicObject <- signatureCellEnrichment(cisTopicObject,
                                          aucellRankings,
                                          selected.signatures='all',
                                          aucMaxRank = AUC_MaxRank*nrow(aucellRankings),
                                          plot=FALSE)

# Save AUC scores to file
outputFile <- paste(outputDirectory,"/",PlotPrefix,"_",PlotSuffix,"_", PlotSuffixGeneSet,"_AUCMaxRank",AUC_MaxRank,"_AUCScores",".txt",sep="")

cell_signatures <- cisTopicObject@cell.data
cell_signatures <- cell_signatures[-c(1,2)]
write.table(cell_signatures, file = outputFile, row.names = T, col.names = T, sep = '\t', quote = F)
