library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

PlotPrefix = "cisTopic50_Epithelial_PrimaryPDAC_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3"



############################################################
# Create ArchRProject with only primary PDAC samples, excluding PAV tumors

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")

# Remove PAV tumors
proj <- proj[proj$Sample_ID %in% c("ST-00017310-T","ST-00019171-T")==F]

# Remove metastases
proj <- proj[grep("-T",proj$Sample_ID)]

# Save new ArchRProject
proj <- saveArchRProject(proj, drop = T, outputDirectory = "~/ArchR/Epithelial_PrimaryPDACOnly_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



########################################################################################################################
# DIFFERENTIAL ACCESSIBILITY: Malignant subtypes vs. acinar/ductal
# In this section, differential accessibility tests will be performed between acinar/ductal cells and each epigenetic subtype 
########################################################################################################################

############################################################
# Differential accessibility: primary classical-NRP vs. acinar/ductal

group1 = "Classical-NRP"
group2 = "Acinar/Ductal"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

classicalNRP_vs_acinar <- peakmarkers

# 7668 acinar
# 12533 classical



############################################################
# Differential accessibility: primary basal vs. acinar/ductal

group1 = "Basal Committed"
group2 = "Acinar/Ductal"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

basal_vs_acinar <- peakmarkers

# 3437 acinar
# 42298 basal



############################################################
# Differential accessibility: primary SIP vs. acinar/ductal

group1 = "Subtype-Intermediate Progenitor"
group2 = "Acinar/Ductal"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.2 & abs(Log2FC) >= 0.5")
peakmarkers <- as.data.frame(markerList)
rownames(peakmarkers) <- 1:length(rownames(peakmarkers))
peakmarkers$group_name[peakmarkers$Log2FC>0] <- group1
peakmarkers$group_name[peakmarkers$Log2FC<0] <- group2
peakmarkers$peak <- paste(peakmarkers$seqnames,peakmarkers$start,peakmarkers$end, sep="_")

intermediate_vs_acinar <- peakmarkers

# 3076 acinar
# 15114 intermediate



########################################################################################################################
# PEAK INTERSECTIONS: Get consensus peak sets
# This section will identify peaks that are more accessible/less accessible in acinar/ductal cells compared to each of the other subtypes
########################################################################################################################

# Peak intersections function
intersectPeaks <- function(peakObjList){
    output <- peakObjList[[1]]
    for(peakset in peakObjList[2:length(peakObjList)]){
        output <- output[output$peak %in% peakset$peak,]
    }
    return(output)
}



# Peak quality filters
filter_fdr = 0.01
filter_log2fc = 1



############################################################
# Intersect peaks: malignant

peakSetList <- list()

filter_group = "Classical-NRP"
peakmarkers <- classicalNRP_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Basal Committed"
peakmarkers <- basal_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Subtype-Intermediate Progenitor"
peakmarkers <- intermediate_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

shared_malignant <- intersectPeaks(peakSetList)
dim(intersectPeaks(peakSetList))

write.table(peakmarkers, file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Consensus_PDAC.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)

# 4178 consensus PDAC



############################################################
# Intersect peaks: acinar/ductal

peakSetList <- list()

filter_group = "Acinar/Ductal"
peakmarkers <- classicalNRP_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Acinar/Ductal"
peakmarkers <- basal_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

filter_group = "Acinar/Ductal"
peakmarkers <- intermediate_vs_acinar
peakmarkers <- peakmarkers[peakmarkers$group_name==filter_group & peakmarkers$FDR<=filter_fdr & abs(peakmarkers$Log2FC)>=filter_log2fc,]
peakSetList[[length(peakSetList)+1]] <- peakmarkers

shared_acinar <- intersectPeaks(peakSetList)
dim(intersectPeaks(peakSetList))

write.table(peakmarkers, file = paste(getOutputDirectory(proj),"/Plots/","DifferentialAccessibility_PrimaryPDAC_Consensus_Acinar.txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)

# 571 consensus acinar



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)



########################################################################################################################
# TF Motif Enrichment
# This section uses ArchR functions to count the number of peaks (in a given peak set) that have a TF motif
# First, we have ArchR identify TF motifs present in the background peak set (all MACS v3 called peaks)
# Second, we use ArchR's peakAnnoEnrichment to count the number of peaks in each consensus peak set with each TF motif
# Third, we perform a hypergeometric test for each TF motif examined and correct the p-values
########################################################################################################################


############################################################
# Function: use ArchR to count the number of query peaks and background peaks containing a specific motif

count_motif_hits <- function(peaks_to_test = peaks_to_test,
                             peakset_name = "Peak_Set",
                             name_annot = "Motif",
                             ArchRProject = proj,
                             SE_object = markersPeaks){
    
    # Flag all peaks in peaks_to_test, unflag all other background peaks
    SE_peaks <- as.data.frame(rowData(SE_object))
    keep <- assay(SE_object,"FDR")
    colnames(keep) <- "keep"
    keep$keep <- 0
    keep$peak <- paste(SE_peaks$seqnames,SE_peaks$start,SE_peaks$end,sep="_")
    keep$keep[keep$peak %in% peaks_to_test$peak] = 1
    keep <- keep[,"keep",drop=F]
    assay(SE_object,"keep", withDimnames = F) <- keep
    
    # Run ArchR peakAnnoEnrichment
    enrichEncode <- peakAnnoEnrichment(
        seMarker = SE_object,
        ArchRProj = ArchRProject,
        peakAnnotation = name_annot,
        cutOff = "keep == 1", # keep only peaks in peaks_to_test
        background = "all"
    )
    
    results <- data.frame(TF = rownames(enrichEncode))
    for(assay in names(enrichEncode@assays)){
        temp1 <- as.data.frame(assay(enrichEncode,assay))
        colnames(temp1) <- "assay"
        results[assay] <- temp1[match(results$TF,rownames(temp1)),"assay"]
    }
    rm(assay,temp1)
    results <- results[order(results$Enrichment,decreasing = T),]
    results <- data.frame(Query_Peak_Set = peakset_name,
                          TF_Motif = results$TF,
                          Enrichment = results$Enrichment,
                          Query_Hits = results$CompareFrequency,
                          Query_PeakCount = results$nCompare,
                          Query_Proportion_Hits = results$CompareProportion,
                          Background_Hits = results$BackgroundFrequency,
                          Background_PeakCount = results$nBackground)
    return(results)
}



############################################################
# Initiate random marker peaks test
# the results of this test do not matter; this is only to generate the SummarizedExperiment object that the function above will take to perform the custom enrichment

group1 = "pg_1"
group2 = "pg_2"

set.seed(1)
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix", 
    groupBy = "cistopic50_pg_k15",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = group1,
    bgdGroups = group2,
    testMethod = "binomial",
    binarize = T
)



############################################################
# Count overlap of TF motifs with consensus peak sets

peaks_to_test <- shared_malignant
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Consensus PDAC Peaks",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
motif_consensusPDAC <- results



peaks_to_test <- shared_acinar
results <- count_motif_hits(peaks_to_test = peaks_to_test,
                            peakset_name = "Consensus Acinar/Ductal Peaks",
                            name_annot = "Motif",
                            SE_object = markersPeaks)
results$TF_Motif <- gsub("_.*","",results$TF_Motif)
motif_consensusAcinar <- results



############################################################
# Perform hypergeometric tests
# Correct p-values using all results

results_all <- rbind(motif_consensusPDAC, motif_consensusAcinar)

balls_drawn_white = results_all$Query_Hits
total_white_balls = results_all$Background_Hits
total_black_balls = results_all$Background_PeakCount - results_all$Background_Hits
total_balls_drawn = results_all$Query_PeakCount

p_ln <- phyper(q = balls_drawn_white - 1,
               m = total_white_balls,
               n = total_black_balls,
               k = total_balls_drawn,
               lower.tail = F,
               log.p = T)
p_log10 <- -(p_ln/log(10))
p_log10_bonferroni <- p_log10 - log10( length(p_log10) )

# Convert negative values to 0
p_log10_bonferroni[p_log10_bonferroni<0] <- 0

results_all$pval_log10 <- p_log10
results_all$pval_log10_Bonferroni <- p_log10_bonferroni



############################################################
# Plot barplot of TF motifs

table <- results_all

features_to_use <- c()

# Get top significant acinar/ductal features
plot <- table[table$Query_Peak_Set=="Consensus Acinar/Ductal Peaks",]
plot$Enrichment[plot$pval_log10_Bonferroni<(-log10(0.05))] = 0
plot <- plot[order(plot$pval_log10, decreasing = T),]
features_to_use <- c(features_to_use, plot$TF_Motif[1:5])

# Get top enriched PDAC features
plot <- table[table$Query_Peak_Set=="Consensus PDAC Peaks",]
plot$Enrichment[plot$pval_log10_Bonferroni<(-log10(0.05))] = 0
plot <- plot[order(plot$Enrichment, decreasing = T),]
features_to_use <- c(features_to_use,plot$TF_Motif[1:5])

plot <- table
plot$significance <- plot$pval_log10_Bonferroni
plot$significance[plot$Query_Peak_Set=="Consensus PDAC Peaks"] <- -plot$significance[plot$Query_Peak_Set=="Consensus PDAC Peaks"]
plot <- plot[plot$TF_Motif %in% features_to_use,]
plot$TF_Motif <- factor(plot$TF_Motif, levels = rev(c("HNF1B","HNF1A","FOXD2","FOXJ2","NFIC","PLAGL1","MYC","MXI1","MYCN","CREB3L2")))

# Set x-axis and enrichment maximum
plot$Enrichment <- log2(plot$Enrichment)
plot$Enrichment[plot$Enrichment<0] <- 0
enrichment_max <- max(plot$Enrichment)
padj_max <- max(abs(plot$significance))

ggplot(plot, aes(x = significance, y = TF_Motif)) +
    geom_bar(stat = "identity", aes(fill = Enrichment)) +
    default_theme +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
        # axis.text.x = element_text(vjust = 3)
    ) +
    
    scale_fill_gradient(low = "white", high = "red",
                        limit = c(0,enrichment_max),
                        guide = guide_colourbar(
                            direction = "vertical",
                            title.theme = element_blank(),
                            label.theme = element_blank(),
                            barwidth=0.75, 
                            barheight = 3,
                            ticks.linewidth = NA)
    ) +
    
    xlim(-padj_max,padj_max) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    geom_hline(yintercept = 5.5, linewidth = 0.5, color = "gray30", linetype = 2)






########################################################################################################################
# Module cCRE enrichment in consensus peak sets

# This section counts the overlap between consensus peak sets and 150 module cCRE peak sets
# The 150 module cCRE peak sets come from the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)
# Hypergeometric tests then determine significance
########################################################################################################################

############################################################
# Query:
# Consensus peak sets

topic_list <- list()

peakmarkers <- shared_malignant
peakmarkers <- GRanges(seqnames = peakmarkers$seqnames, IRanges(start = peakmarkers$start, end = peakmarkers$end))
topic_list[["Consensus PDAC Peaks"]] <- peakmarkers

peakmarkers <- shared_acinar
peakmarkers <- GRanges(seqnames = peakmarkers$seqnames, IRanges(start = peakmarkers$start, end = peakmarkers$end))
topic_list[["Consensus Acinar/Ductal Peaks"]] <- peakmarkers



############################################################
# Enrichment:
# Load 150 module cCRE peak sets in the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)

peakset_file_path <- "~/peak_sets/FetalAdultAtlas_Zhang_2021/module_cCREs"
peakset_files <- list.files(path = peakset_file_path, full.names = F)
peakset_list <- list()
for(unique_peakset in peakset_files){
    peakmarkers <- read.table(paste(peakset_file_path, unique_peakset, sep = "/"), header = F, sep = '\t')
    peakmarkers <- GRanges(seqnames = peakmarkers$V1, IRanges(start = peakmarkers$V2, end = peakmarkers$V3))
    peakset_list[[unique_peakset]] <- peakmarkers
    rm(peakmarkers)
}



############################################################
# Background:
# Load background MACS v3 called peak set

peakmarkers <- read.table("~/PeakCalling/MACS3/AllSamples_TSSMixed_nFrags1k_FRIS0.4/AllSamplesTSSMixednFrags1kFRIS0.4_MACS3.bed")
peakmarkers <- GRanges(seqnames = peakmarkers$V1, IRanges(start = peakmarkers$V2, end = peakmarkers$V3))
background <- peakmarkers



#########################################################
# Function to perform hypergeometric tests on list of peak sets

peakSetEnrichment <- function(topic_list,
                              peakset_list,
                              bgd){
    output <- data.frame()
    
    for(unique_topic in names(topic_list)){
        for(unique_peakset in names(peakset_list)){
            peakset_topic <- topic_list[[unique_topic]]
            peakset_cCRE <- peakset_list[[unique_peakset]]
            peakset_all_topics <- bgd
            topic_length = length(peakset_topic)
            
            total_balls_drawn = topic_length
            balls_drawn_white = length(peakset_topic[overlapsAny(query = peakset_topic, subject = peakset_cCRE)])
            total_white_balls = length(peakset_all_topics[overlapsAny(query = peakset_all_topics, subject = peakset_cCRE)])
            total_black_balls = length(peakset_all_topics) - total_white_balls
            
            p_ln <- phyper(q = balls_drawn_white - 1,
                           m = total_white_balls,
                           n = total_black_balls,
                           k = total_balls_drawn,
                           lower.tail = F,
                           log.p = T)
            p_log10 <- -(p_ln/log(10))
            p_archr <- p_log10
            enrichment <- (balls_drawn_white/total_balls_drawn)/(total_white_balls/length(peakset_all_topics))
            p_archr
            enrichment
            
            output[length(rownames(output))+1,"Test_Set"] <- unique_topic
            output[length(rownames(output)),"Test_Peaks"] <- total_balls_drawn
            output[length(rownames(output)),"Test_Peak_Hits"] <- balls_drawn_white
            output[length(rownames(output)),"Bgd_Peaks"] <- total_white_balls + total_black_balls
            output[length(rownames(output)),"Bgd_Peak_Hits"] <- total_white_balls
            
            output[length(rownames(output)),"Peak_Set"] <- unique_peakset
            output[length(rownames(output)),"Enrichment"] <- enrichment
            output[length(rownames(output)),"mlog10p"] <- p_archr
            output[length(rownames(output)),"p_val"] <- exp(p_ln)
        }
    }
    return(output)
}



#########################################################
# Perform enrichments

test <- peakSetEnrichment(topic_list = topic_list,
                          peakset_list = peakset_list,
                          bgd = background)



#########################################################
# Calculate p-values

results_all <- test

# Correct p-values using all results
balls_drawn_white = results_all$Test_Peak_Hits
total_white_balls = results_all$Bgd_Peak_Hits
total_black_balls = results_all$Bgd_Peaks - results_all$Bgd_Peak_Hits
total_balls_drawn = results_all$Test_Peaks

p_ln <- phyper(q = balls_drawn_white - 1,
               m = total_white_balls,
               n = total_black_balls,
               k = total_balls_drawn,
               lower.tail = F,
               log.p = T)
p_log10 <- -(p_ln/log(10))
p_log10_bonferroni <- p_log10 - log10( length(p_log10) )

# Convert negative values to 0
p_log10_bonferroni[p_log10_bonferroni<0] <- 0

results_all$mlog10p <- p_log10
results_all$mlog10Padj <- p_log10_bonferroni



#########################################################
# Tabulate results

table <- data.frame(Query_Peak_Set = results_all$Test_Set,
                    CisRegulatoryModule_Peak_Set = results_all$Peak_Set,
                    pval_log10 = results_all$mlog10p,
                    pval_log10_Bonferroni = results_all$mlog10Padj,
                    Enrichment = results_all$Enrichment,
                    Background_PeakCount = results_all$Bgd_Peaks,
                    Background_Hits = results_all$Bgd_Peak_Hits,
                    Background_Proportion_Hits = results_all$Bgd_Peak_Hits/results_all$Bgd_Peaks,
                    Query_PeakCount = results_all$Test_Peaks,
                    Query_Hits = results_all$Test_Peak_Hits,
                    Query_Proportion_Hits = results_all$Test_Peak_Hits/results_all$Test_Peaks)
table$Query_Peak_Set <- gsub("Topic_","Topic ", table$Query_Peak_Set)
table$Query_Peak_Set <- gsub("\\.bed","", table$Query_Peak_Set)

table$CisRegulatoryModule_Peak_Set <- gsub("\\.txt","",table$CisRegulatoryModule_Peak_Set)
table$CisRegulatoryModule_Peak_Set <- paste("module",table$CisRegulatoryModule_Peak_Set,sep="")

table <- table[order(table$Enrichment, decreasing = T),]



#########################################################
# Plot dotplot of module cCREs (Fig. 3d)

plot <- table
plot$Module <- plot$CisRegulatoryModule_Peak_Set
plot$Module <- gsub("^module[0-9]*_","",plot$Module)
plot$Module <- gsub("_"," ",plot$Module)
plot$Module <- gsub("and","&",plot$Module)
plot$Module <- gsub(" or ","/",plot$Module)
plot$Module[plot$Module=="Pulm Epithelial"] <- "Pulmonary Epithelial"


ggplot(plot[plot$Query_Peak_Set=="Consensus Acinar/Ductal Peaks",],
       aes(x = pval_log10_Bonferroni, y = Module)) +
    geom_point(aes(color = Enrichment), size = 1.60) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        axis.ticks.x = element_line(linewidth = line_width*0.5, color = "black"),
        axis.ticks.y = element_line(linewidth = line_width*0.5, color = "black"),
        
        axis.text.x = text_settings,
        axis.text.y = text_settings,
        
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        
        # panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linewidth = line_width, fill = NA),
        
        plot.background = element_blank(),
        plot.title = title_settings,
        
        legend.text = text_settings,
        legend.title = title_settings
    ) +
    scale_color_gradientn(colours = colorRampPalette(ArchR::ArchRPalettes$coolwarm)(1000),
                          # limit = c(0, max_enrichment),
                          guide = guide_colourbar(
                              direction = "vertical",
                              title.theme = element_blank(),
                              # label.theme = element_blank(),
                              barwidth=0.75, 
                              barheight = 3,
                              ticks.linewidth = NA)) +
    xlim(0,42)



ggplot(plot[plot$Query_Peak_Set=="Consensus PDAC Peaks",],
       aes(x = pval_log10_Bonferroni, y = Module)) +
    geom_point(aes(color = Enrichment), size = 1.60) +
    theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        axis.ticks.x = element_line(linewidth = line_width*0.5, color = "black"),
        axis.ticks.y = element_line(linewidth = line_width*0.5, color = "black"),
        
        axis.text.x = text_settings,
        axis.text.y = text_settings,
        
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        
        # panel.grid = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(linewidth = line_width, fill = NA),
        
        plot.background = element_blank(),
        plot.title = title_settings,
        
        legend.text = text_settings,
        legend.title = title_settings
    ) +
    scale_color_gradientn(colours = colorRampPalette(ArchR::ArchRPalettes$coolwarm)(1000),
                          # limit = c(0, max_enrichment),
                          guide = guide_colourbar(
                              direction = "vertical",
                              title.theme = element_blank(),
                              # label.theme = element_blank(),
                              barwidth=0.75, 
                              barheight = 3,
                              ticks.linewidth = NA)) +
    xlim(0,155)
