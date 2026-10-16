library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)

library(IRanges)
library(GenomicRanges)
library(ggrepel)
library(pheatmap)

proj <- loadArchRProject("~/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/")



########################################################################################################################
# Cell type enhancer enrichment in topics
# This section counts the overlap between 222 cell type peak sets with topic peaks and background peaks
# The 222 cell type peak sets come from the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)
# Hypergeometric tests then determine significance
########################################################################################################################

############################################################
# Query:
# Load all 50 epithelial/malignant topics

# Topic Peak Sets
topic_list <- list()
for(unique_topic in paste("Topic_",1:50,".bed",sep="")){
    peakmarkers <- read.table(paste("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/topics50_bedfiles_binarized/", unique_topic, sep=""), header = F, sep = '\t')
    peakmarkers <- GRanges(seqnames = peakmarkers$V1, IRanges(start = peakmarkers$V2, end = peakmarkers$V3))
    topic_list[[unique_topic]] <- peakmarkers
    rm(peakmarkers)
}



############################################################
# Enrichment:
# Load 222 cell type peak sets in the chromatin atlas by Zhang K and Hocker JD, et al. (2021) (https://doi.org/10.1016/j.cell.2021.10.024)

peakset_file_path <- "~/peak_sets/FetalAdultAtlas_Zhang_2021/celltype_cCREs"
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
# Load background set of all unique peaks across all 50 topics

peakmarkers <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/topics50_bedfiles_binarized/EpithelialAllSamples_cisTopic50_TopicPeaksOnly.bed", header = F, sep = '\t')
peakmarkers <- GRanges(seqnames = peakmarkers$V1, IRanges(start = peakmarkers$V2, end = peakmarkers$V3))
peakset_all_topics <- peakmarkers



#########################################################
# Function to perform hypergeometric tests on list of topics

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
                          bgd = peakset_all_topics)



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
                    Celltype_cCRE_Peak_Set = results_all$Peak_Set,
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

table$Celltype_cCRE_Peak_Set <- gsub("\\.txt","",table$Celltype_cCRE_Peak_Set)
table$Celltype_cCRE_Peak_Set <- gsub("_"," ",table$Celltype_cCRE_Peak_Set)





########################################################################################################################
# Example scatterplot of topic enrichment of 222 cell type peak sets
########################################################################################################################

#########################################################
# Theme Settings

line_width <- 0.75
point_size <- 0.65
text_size <- 6
title_size <- 6

line_settings <- element_line(linewidth = line_width, color = "black")
text_settings <- element_text(colour = "black", size = text_size)
title_settings <- element_text(colour = "black", size = title_size, vjust = -3)


# Default Theme
default_theme <- theme(
    axis.title.x = element_blank(),
    axis.title.y = text_settings,
    
    axis.ticks.x = element_line(linewidth = line_width*0.5, color = "black"),
    axis.ticks.y = element_line(linewidth = line_width*0.5, color = "black"),
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(linewidth = line_width, fill = NA),
    
    plot.background = element_blank(),
    plot.title = title_settings,
    
    legend.text = text_settings,
    legend.title = title_settings
)



#########################################################
# Plot Topic 42 results

plot <- table
plot <- plot[plot$Query_Peak_Set=="Topic 42",]

label_tfs <- NA

plot$label = plot$Celltype_cCRE_Peak_Set
plot$label[ (plot$pval_log10_Bonferroni < 600) & plot$Celltype_cCRE_Peak_Set %in% label_tfs==F ] <- NA

ggplot(plot, aes(y = pval_log10_Bonferroni, x = Enrichment)) +
    geom_point(color = "#e1468a", stroke = 0, size = 1) +
    # geom_text(aes(label = label), size = 6/.pt, nudge_y = 0.1) +
    default_theme +
    # geom_vline(xintercept = -log10(0.05),linetype = 2, linewidth = 0.2) +
    geom_text_repel(aes(label = label),
                    size = 6/.pt,
                    min.segment.length = 0.0,
                    max.overlaps = 200,
                    segment.size = 0.5/.pt,
                    force = 2,
                    box.padding = 0.2,
                    max.time = 3,
                    max.iter = 100000,
                    nudge_x = -0.1) +
    ylab("-log10(Adjusted P)") + theme(aspect.ratio = 1.3) 






########################################################################################################################
# Plot p-value rank of mucinous GI, esophageal epithelial, and tuft lineages per topic
########################################################################################################################

classicalNRP_topics <- c(3,10,18,29,30,34,37,40,41,44,45,46)
classicalNRP_topics <- paste("Topic ", classicalNRP_topics, sep = "")

basal_topics <- c(9,22,23,47)
basal_topics <- paste("Topic ", basal_topics, sep = "")

intermed_topics <- c(6,20,42)
intermed_topics <- paste("Topic ", intermed_topics, sep = "")



#########################################################
# Remove spurious topics with fewer than 5 significant results

results <- data.frame(
    TF = table$Celltype_cCRE_Peak_Set,
    Peakset = table$Query_Peak_Set,
    Enrichment = table$Enrichment,
    mlog10p = table$pval_log10,
    mlog10Padj = table$pval_log10_Bonferroni
)

min_significant_hits <- 5
bad_topics <- as.data.frame(table(results$Peakset[results$mlog10Padj<=(-log10(0.05))]))
bad_topics$good_topics <- length(unique(results$TF))-bad_topics$Freq
bad_topics <- bad_topics$Var1[bad_topics$good_topics<min_significant_hits]

results <- results[results$Peakset %in% bad_topics==F,]



#########################################################
# Within each topic, rank cell type peak sets by p-value

output <- data.frame(celltype = unique(results$TF))
for(unique_topic in unique(results$Peakset)){
    test <- results[results$Peakset==unique_topic,]
    test$p_rank <- rank(test$mlog10p)
    output[,unique_topic] <- test$p_rank[match(output$celltype,test$TF)]
}
rownames(output) <- output$celltype
output <- output[,-1,drop=F]



#########################################################
# Get a list of mucinous GI lineages

GI_lineages <- unique(
    grep("Goblet|Foveolar|Enterocyte|Colon Epithelial",results$TF,value = T,ignore.case = T)
)

# Exclude airway goblet and fetal cell types
GI_lineages <- grep("Sm Ms|Neuroendocrine|Airway|Fetal",GI_lineages, invert = T, value = T)



#########################################################
# Calculate average rank of mucinous GI lineages for each topic

plot <- data.frame()
for(unique_topic in unique(colnames(output))){
    plot <- rbind(plot,
                  data.frame(lineage = rownames(output),
                             topic = unique_topic,
                             rank = output[[unique_topic]])
    )
}


plot <- plot[plot$lineage %in% c("Tuft","Esophageal Epithelial",GI_lineages),]
plot <- plot[plot$topic %in% c(basal_topics,classicalNRP_topics,intermed_topics),]


plot$subtype[plot$topic %in% basal_topics] <- "BC"
plot$subtype[plot$topic %in% classicalNRP_topics] <- "C-NRP"
plot$subtype[plot$topic %in% intermed_topics] <- "SIP"


plot$lineage_type <- plot$lineage
plot$lineage_type[plot$lineage_type %in% GI_lineages] <- "Mucinous GI (Avg.)"


plot <- plot %>%
    group_by(topic, lineage_type) %>%
    mutate(avg_rank = mean(rank)) %>%
    select(-lineage,-rank) %>%
    unique()



#########################################################
# Plot p-value rank as a heatmap

heatmap <- data.frame()
for(unique_topic in unique(plot$topic)){
    for(unique_lineage in unique(plot$lineage_type))
        heatmap[unique_lineage,unique_topic] <- plot$avg_rank[plot$topic==unique_topic & plot$lineage_type==unique_lineage]
}


topic_order <- c(
    "Topic 23",
    "Topic 9",
    "Topic 22",
    "Topic 47",
    "Topic 42",
    "Topic 20",
    "Topic 6",
    "Topic 41",
    "Topic 18",
    "Topic 40",
    "Topic 44",
    "Topic 45",
    "Topic 10",
    "Topic 34")
topic_order <- topic_order[topic_order %in% colnames(heatmap)]

heatmap <- heatmap[c("Esophageal Epithelial","Mucinous GI (Avg.)","Tuft"),topic_order,]

heatmap_annot <- data.frame(Subtype = colnames(heatmap))
rownames(heatmap_annot) <- heatmap_annot$Subtype
heatmap_annot$Subtype[rownames(heatmap_annot) %in% classicalNRP_topics] <- "Classical-NRP"
heatmap_annot$Subtype[rownames(heatmap_annot) %in% basal_topics] <- "Basal Committed"
heatmap_annot$Subtype[rownames(heatmap_annot) %in% intermed_topics] <- "Subtype-Intermediate Progenitor"

colnames(heatmap) <- gsub("Topic ", "T", colnames(heatmap))
rownames(heatmap_annot) <- gsub("Topic ", "T", rownames(heatmap_annot))

pheatmap(heatmap,
         cluster_cols = F,
         cluster_rows = F,
         annotation_col = heatmap_annot,
         annotation_colors = list(Subtype = c("Basal Committed" = "#bd3531",
                                              "Subtype-Intermediate Progenitor" = "#e1468a",
                                              "Classical-NRP" = "#f79747")),
         cellwidth = 8,
         cellheight = 8,fontsize = 6,
         border_color = NA,angle_col = 90,
         color = colorRampPalette( ArchR::ArchRPalettes$greenBlue )(1000))


