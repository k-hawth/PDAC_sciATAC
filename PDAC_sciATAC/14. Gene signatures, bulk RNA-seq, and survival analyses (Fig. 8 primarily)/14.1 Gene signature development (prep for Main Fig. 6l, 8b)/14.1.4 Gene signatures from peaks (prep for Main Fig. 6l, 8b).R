# Gene signature development

# In this section, peak sets are translated into proxy "gene signatures" for transcriptomic analysis. This is done by linking peaks to their nearest gene.

# Because some of the peak sets contain thousands of peaks, we use only up to 200 genes from these linkages. The top 200 are prioritized based on the cisTopic probabilistic model data we exported earlier.



library(rtracklayer)


############################################################
# Load cisTopic data

# FOR MALIGNANT: Load mean peak probability by primary tumor subtype. This will be used to prioritize the top 200 genes based on the score of their linked peak
region_scores <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic50_Epithelial_RegionCellScores_MeanByCelltypeSep.txt", sep = '\t', header = T, comment.char = "", quote = "")


# FOR STROMAL: Load probability of peak belonging to each stromal topic. This will be used to prioritize the top 200 genes based on the score of their linked peak
caf_topicPeaks <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/Stromal_cisTopic_PeakTopicProbability.txt", sep = '\t', header = T, comment.char = "", quote = "")



############################################################
# Function to link peaks with nearest genes

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
        print("ChIPseeker Error: Unequal peakAnno and peakset")
    }
}



############################################################
# Load gene annotation
# This will be used to exclude "pseudogene" gene types

# GENCODE V41
annotation <- import("~/reference_genomes/gencode/gencode.v41.annotation.gtf.gz")
annotation <- as.data.frame(annotation@elementMetadata@listData)

# Prepare annotation object for converting ensembl ID to gene name
annotation$gene_id_short <- gsub("\\..*","",annotation$gene_id)  # Remove Ensembl ID end tags
annotation <- annotation[c("gene_id","gene_id_short","gene_type","gene_name")]  # Keep only the info of interest
annotation <- unique(annotation)



############################################################
# Gene signature filters

# Exclude peak-gene links exceeding 50kb in distance
max_distance_to_TSS <- 50000

# Exclude pseudogenes
gene_types_to_keep <- unique(annotation$gene_type)
gene_types_to_keep <- grep("pseudogene", gene_types_to_keep, value = T, invert = T)




########################################################################################################################
# PEAK-TO-GENE SIGNATURES: PDAC epigenetic subtypes (SIP, classical-NRP, and basal committed)
########################################################################################################################

table <- data.frame()



############################################################
# Convert SIP peaks to SIP gene signature

feature_name <- "SIP_Signature"

peakmarkers <- read.table("~/ArchR/Epithelial_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50_PrimaryPDACOnly/Plots/DifferentialAccessibility_PrimaryPDAC_Specific_Intermediate.txt", header = T, sep = '\t', comment.char = "", quote = "")

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "seqnames", "start", "end")

# Annotate with mean probability of peak in primary SIP
peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_SIP <- region_scores$Subtype.Intermediate.Progenitor[match(peakmarkers$peak, rownames(region_scores))]
peakmarkers <- peakmarkers[order(peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_SIP, decreasing = T),]

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL






############################################################
# Convert basal committed peaks to BC gene signature

feature_name <- "BC_Signature"

peakmarkers <- read.table("~/ArchR/Epithelial_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50_PrimaryPDACOnly/Plots/DifferentialAccessibility_PrimaryPDAC_Specific_Basal.txt", header = T, sep = '\t', comment.char = "", quote = "")

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "seqnames", "start", "end")

# Annotate with mean probability of peak in primary basal committed
peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_BasalCommitted <- region_scores$Basal.Committed[match(peakmarkers$peak, rownames(region_scores))]
peakmarkers <- peakmarkers[order(peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_BasalCommitted, decreasing = T),]

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL





############################################################
# Convert classical-NRP peaks to C-NRP gene signature

feature_name <- "C_NRP_Signature"

peakmarkers <- read.table("~/ArchR/Epithelial_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50_PrimaryPDACOnly/Plots/DifferentialAccessibility_PrimaryPDAC_Specific_ClassicalNRP.txt", header = T, sep = '\t', comment.char = "", quote = "")

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "seqnames", "start", "end")

# Annotate with mean probability of peak in primary classical-NRP
peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_ClassicalNRP <- region_scores$Classical.NRP[match(peakmarkers$peak, rownames(region_scores))]
peakmarkers <- peakmarkers[order(peakmarkers$cisTopic_MeanProbabilityOfPeak_in_PrimaryTumor_ClassicalNRP, decreasing = T),]

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL






########################################################################################################################
# PEAK-TO-GENE SIGNATURES: Stromal topics associated with CAF subtypes (LRRC15+ myCAFs, NGFR+ myCAFs, and iCAFs)
########################################################################################################################


############################################################
# LRRC15+ myCAF epigenetic subtype signature


feature_name <- "LRRC15_myCAF_Signature_Topic15"

peakmarkers <- caf_topicPeaks[caf_topicPeaks$Topic=="Topic 15",]
peakmarkers <- peakmarkers[order(peakmarkers$Peak_Topic_Probability, decreasing = T),]

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "peakChromosome", "peakStart", "peakEnd")

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL



############################################################
# iCAF epigenetic subtype signature

feature_name <- "iCAF_Signature_Topic11"

peakmarkers <- caf_topicPeaks[caf_topicPeaks$Topic=="Topic 11",]
peakmarkers <- peakmarkers[order(peakmarkers$Peak_Topic_Probability, decreasing = T),]

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "peakChromosome", "peakStart", "peakEnd")

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL



############################################################
# NGFR+ myCAF epigenetic subtype signature

feature_name <- "NGFR_myCAF_Signature_Topic22"

peakmarkers <- caf_topicPeaks[caf_topicPeaks$Topic=="Topic 22",]
peakmarkers <- peakmarkers[order(peakmarkers$Peak_Topic_Probability, decreasing = T),]

# Annotate with nearest gene info
peakmarkers <- chipseek(peakmarkers, "peakChromosome", "peakStart", "peakEnd")

# Remove peaks that do not meet inclusion criterion
peakmarkers <- peakmarkers[abs(peakmarkers$distanceToTSS)<=max_distance_to_TSS,]
peakmarkers$gene_type <- annotation$gene_type[match(peakmarkers$ENSEMBL,annotation$gene_id_short)]
peakmarkers <- peakmarkers[peakmarkers$gene_type %in% gene_types_to_keep,]

peakmarkers <- unique(peakmarkers[,c("ENSEMBL","SYMBOL"),drop=F])

# Limit to top 200 genes
if(length(rownames(peakmarkers))>200){
    peakmarkers <- peakmarkers[1:200,]
}

# Add to gene signatures table
table[1:length(rownames(peakmarkers)),paste(feature_name,"Ensembl",sep="_")] <- peakmarkers$ENSEMBL
table[1:length(rownames(peakmarkers)),paste(feature_name,"GeneSymbol",sep="_")] <- peakmarkers$SYMBOL




############################################################
# Export results

write.table(table, file = "~/gene_sets/PDAC/sciATAC_Gene_Signatures.txt", row.names = F, col.names = T, sep = '\t', quote = F)


