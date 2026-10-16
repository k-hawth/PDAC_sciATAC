library(cisTopic)
library(Matrix)
library(ggplot2)
set.seed(1)

library(pheatmap)


########################################################################################################################
# SET VARIABLES
########################################################################################################################

PlotPrefix <- "Fibroblast_cisTopic27"
PlotSuffix <- "WangY2021_AUC"


# Specify folder where outputs will go. Folder will be created if one does not exist
outputDirectory <- "~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/WangY2021_AUC"


cisTopicObject <- readRDS("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/cisTopic_alpha50_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR.cistopicObject_Binarized.rds")
cisTopic_modelNumber <- 27


peaks <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/AllCells_MACS2ArchR/Plots/Peaks_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4.bed")


ChIPSeeker_max_distToTSS <- 50000



########################################################################################################################
# PREPARE GENE SET LIST
########################################################################################################################

# Stromal gene signatures
signatures <- as.data.frame(readxl::read_xlsx("~/gene_sets/PDAC/PDAC_Gene_Signatures.xlsx", sheet = 3))
signatures <- signatures[,grep("WangY2021",colnames(signatures)),drop=F]

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


# Output AUC scores
auc <- signaturesHeatmap(getSignaturesRegions(cisTopicObject, signatures = signatures, labels=names(signatures), minOverlap = 0.4))
auc <- auc@ht_list[["Normalised AUC score"]]@matrix
outputFile <- paste(outputDirectory,"/",PlotPrefix,"_",PlotSuffix,"_AUCScoresForHeatmap.txt",sep="")
write.table(auc, file = outputFile, row.names = T, col.names = T, sep = '\t', quote = F)



########################################################################################################################
# HEATMAP: Plot per-topic enrichment scores
########################################################################################################################


heatmap <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/WangY2021_AUC/Fibroblast_cisTopic27_WangY2021_AUC_AUCScoresForHeatmap.txt", header = T)


pheatmap(heatmap, color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000),
         clustering_method = "ward",
         treeheight_row = 10,
         treeheight_col = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         labels_col = gsub("Topic","T",colnames(heatmap)),
         angle_col = 90)