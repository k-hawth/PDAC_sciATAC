library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")

# Create new ArchRProject
proj <- saveArchRProject(proj, "~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50_TopicPeaksBgd")



############################################################
# Load all peaks from all topics

read_topics <- function(topic_file_directory,
                        topic_file_count,
                        topic_file_prefix = "Topic_",
                        topic_file_suffix = ".bed"){
    topic_list <- list()
    for(topic_file in paste(topic_file_prefix,1:topic_file_count,topic_file_suffix, sep = "")){
        topic_list[[topic_file]] <- read.table(file = paste(topic_file_directory,topic_file,sep = ""), header = F, sep = '\t', quote = "", comment.char = "")
    }
    return(topic_list)
}

cat_topics <- function(topic_list = topic_list){
    for(topic in names(topic_list)){
        topic_df <- topic_list[[topic]]
        topic_df$Topic = topic
        topic_list[[topic]] <- topic_df
    }
    cat <- topic_list[[1]]
    for(topic in names(topic_list)){
        cat <- rbind(cat,topic_list[[topic]])
    }
    cat <- unique(cat)
    cols <- colnames(cat)
    cols <- cols[cols!="Topic"]
    cat <- cat[c("Topic",cols)]
    return(cat)
}


catTopics <- read_topics("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/topics50_bedfiles_binarized/",
                         topic_file_count = 50,
                         topic_file_prefix = "Topic_",
                         topic_file_suffix = ".bed")
catTopics <- cat_topics(catTopics)






########################################################################################################################
# TF Motif Enrichment
# This section uses ArchR functions to count the number of peaks (in a given peak set) that have a TF motif
# First, we have ArchR identify TF motifs present in the background peak set (all unique peaks across all topics)
# Second, we use ArchR's peakAnnoEnrichment to count the number of peaks in each topic with each TF motif
# Third, we perform a hypergeometric test for each topic and compare topic TF motif hits to background TF motif hits
########################################################################################################################



############################################################
# Add topic peaks as background peak set

peaks <- unique(catTopics[,2:4])
peaks <- GRanges(seqnames = peaks$V1, ranges = IRanges(start = peaks$V2, end = peaks$V3))
proj <- addPeakSet(ArchRProj = proj, peakSet = peaks,force=T)
proj <- addPeakMatrix(proj,force=T)

# Call background peaks
proj <- addBgdPeaks(proj,method = "ArchR", force=T)

# Add cisbp motif annotations
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=T)



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
# Count TF motifs in each topic

topic_TF_motifs <- data.frame()
for(unique_topic in unique(catTopics$Topic)){
    peaks_to_test <- catTopics[catTopics$Topic==unique_topic,]
    peaks_to_test <- peaks_to_test[,2:4]
    colnames(peaks_to_test) <- c("seqnames","start","end")
    peaks_to_test$peak <- paste(peaks_to_test$seqnames,peaks_to_test$start,peaks_to_test$end,sep="_")
    
    results <- count_motif_hits(peaks_to_test = peaks_to_test,
                                peakset_name = unique_topic,
                                name_annot = "Motif",
                                SE_object = markersPeaks)
    topic_TF_motifs <- rbind(topic_TF_motifs, results)
}



############################################################
# Perform hypergeometric tests

results_all <- topic_TF_motifs
results_all$Query_Peak_Set <- gsub("_|\\.bed$"," ",results_all$Query_Peak_Set)
results_all$TF_Motif <- gsub("_.*","",results_all$TF_Motif)                       



# Correct p-values using all results
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
# Save ArchRProject

proj <- saveArchRProject(proj)






########################################################################################################################
# Plot heatmap of curated TF motifs
########################################################################################################################

# Make matrix of corrected p-values
matrix_mlog10Padj <- data.frame()
for(unique_feature in unique(results_all$TF_Motif)){
    for(unique_peakset in unique(results_all$Query_Peak_Set)){
        matrix_mlog10Padj[unique_feature,unique_peakset] <- results_all$pval_log10_Bonferroni[results_all$TF_Motif==unique_feature & results_all$Query_Peak_Set==unique_peakset]
    }
    rm(unique_feature, unique_peakset)
}



# Scale mlog10Padj matrix manually
colScaleMinMax <- function(x){
    (x - min(x))/(max(x) - min(x))
}

# Set TFs to columns
mat <- matrix_mlog10Padj
mat <- as.data.frame(t(mat))

# Scale columns min-max (0-100)
for(unique_colname in colnames(mat)){
    mat[[unique_colname]] <- colScaleMinMax(mat[[unique_colname]])
}
mat <- mat * 100
mat <- as.data.frame(t(mat))



############################################################
# Set scaled 0-100 mlog10Padj matrix values to 0, if non-significant

for(unique_feature in unique(rownames(matrix_mlog10Padj))){
    for(unique_peakset in unique(colnames(matrix_mlog10Padj))){
        if(matrix_mlog10Padj[unique_feature,unique_peakset] < (-log10(0.05)) ){
            mat[unique_feature,unique_peakset] = 0
        }
    }
    rm(unique_feature, unique_peakset)
}



############################################################
# Limit to subtype-associated topics

colnames(mat) <- gsub(" ","",colnames(mat))
colnames(mat) <- gsub("^Topic","",colnames(mat))

# Order matrix by topic subtype
mat <- mat[,c("23","9","22","47","42","20","6","34","18","41","10","44","30","45","46","29","3","37","40")]



############################################################
# ANNOTATION DF: subtype-associated topics

heatmap_annot_col <- data.frame(col = colnames(mat))
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% c("3","10","18","29","30","34","37","40","41","44","45","46")] <- "Classical-NRP"
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% c("6","20","42")] <- "Subtype-Intermediate Progenitor"
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% c("9","22","47","23")] <- "Basal Committed"
rownames(heatmap_annot_col) <- heatmap_annot_col$col
heatmap_annot_col <- heatmap_annot_col[,c("Subtype"),drop=F]



############################################################
# ANNOTATION COLORS: PDAC signatures

annot_colors_col <- list(
    "Subtype" = c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a")
)



############################################################
# Plot heatmap

# List of features to plot
tfs_to_show <- c("ZEB1","SNAI1","TWIST1", # EMT
                 "NEUROD1","ASCL1","MAFA", # Endocrine
                 "PTF1A","NR5A2", # Exocrine
                 "HNF4A","HNF4G","HES1","PDX1","MNX1","FOXA2","HNF1B", # Progenitor
                 "JUN","BACH1","FOS", # AP-1
                 "GATA6","GATA4","ONECUT2","ELF3") # Epithelial

cairo_pdf(width = 4.6, height = 2.17,
          file = "~/cistopic_heatmap_TF_Motif.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(as.matrix(mat[tfs_to_show,,drop=F]),
                         angle_col = "90",
                         show_colnames = T,
                         show_rownames = T,
                         cluster_cols = F,
                         cluster_rows = F,
                         color = colorRampPalette(ArchR::ArchRPalettes$horizon)(10000),
                         border_color = NA,
                         annotation_col = heatmap_annot_col,
                         annotation_colors = annot_colors_col,
                         cellwidth = 7,
                         cellheight = 7,
                         fontsize = 6,
                         annotation_legend = F,
                         
                         heatmap_legend_param = list(border = FALSE),
                         row_names_side="left",
                         row_dend_side = "right",
                         column_names_side = "bottom"
)

dev.off()