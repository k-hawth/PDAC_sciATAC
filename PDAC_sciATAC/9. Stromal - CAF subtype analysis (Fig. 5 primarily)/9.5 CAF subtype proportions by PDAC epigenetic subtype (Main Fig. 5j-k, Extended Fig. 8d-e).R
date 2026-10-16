library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

library(ggsignif)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/")



############################################################
# Load table of primary tumors and their epigenetic subtype

subtype <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_PrimaryTumorEpigeneticSubtype.txt", header = T, sep = '\t')



############################################################
# Calculate CAF cluster proportions

# Get proportions by this cellColData column
clusterName <- "celltype"

clusterProp <- data.frame(Cell = proj$cellNames,
                          Sample = proj$Sample_ID,
                          Cluster = proj@cellColData[,clusterName]
)

# Retain CAFs only
clusterProp <- clusterProp[clusterProp$Cluster!="Endothelial",]

# Add total cell count per sample
clusterProp <- clusterProp %>%
    group_by(Sample) %>%
    mutate(Total_Count = n())

# Add cluster counts per sample
output <- data.frame()
for(unique_sample in unique(clusterProp$Sample)){
    for(unique_clust in unique(clusterProp$Cluster)){
        output[length(rownames(output))+1,"Sample"] <- unique_sample
        output[length(rownames(output)),"Cluster"] <- unique_clust
        output[length(rownames(output)),"Cluster_Count"] <- length(rownames(clusterProp[clusterProp$Sample==unique_sample & clusterProp$Cluster==unique_clust,]))
    }
    rm(unique_sample,unique_clust)
}
output$Total_Count <- clusterProp$Total_Count[match(output$Sample,clusterProp$Sample)]
clusterProp <- output
clusterProp$Proportion <- clusterProp$Cluster_Count/clusterProp$Total_Count
rm(output)



############################################################
# Filter for samples that have an epigenetic subtype and minimum cell count

# Filter for samples with min CAF cell count
clusterProp <- clusterProp[clusterProp$Total_Count>=10,]

# Filter for samples with an epigenetic subtype
clusterProp$Subtype <- subtype$Subtype[match(clusterProp$Sample,subtype$Sample)]
clusterProp <- clusterProp[clusterProp$Subtype %in% c("Basal Committed","Classical-NRP","Subtype-Intermediate Progenitor"),]



############################################################
# ggplot theme

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


colorPalette <- c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a"
)




############################################################
# Plot LRRC15+ myCAF proportions

plot <- clusterProp
plot <- plot[plot$Cluster=="LRRC15+ myCAF",]

ylim_adjustment <- 1.75  # set to 1 for no change
step_adjustment <- c(0,0.25,0.18)

ggplot(plot, aes(x = Subtype, y = Proportion)) +
    geom_boxplot(outlier.alpha = 0, aes(fill = Subtype), color = "black", lwd = line_width*0.5, fatten = line_width*2) +
    geom_point(size = 0.65) +
    geom_signif(comparisons = list(
        c("Basal Committed","Subtype-Intermediate Progenitor"),
        c("Classical-NRP","Subtype-Intermediate Progenitor"),
        c("Basal Committed","Classical-NRP")),
        vjust = -0.1,
        step_increase = step_adjustment,
        test = "wilcox.test",
        textsize = 6/.pt,
        margin_top = 0.15) +
    scale_fill_manual(values = colorPalette) +
    # ylim(min(plot$feature),abs(max(plot$feature)-min(plot$feature))*ylim_adjustment) +
    ylab("Count / Total CAFs") +
    default_theme +
    theme(
        legend.position = "none"
    )



############################################################
# Plot NGFR+ myCAF proportions

plot <- clusterProp
plot <- plot[plot$Cluster=="NGFR+ myCAF",]

ggplot(plot, aes(x = Subtype, y = Proportion)) +
    geom_boxplot(outlier.alpha = 0, aes(fill = Subtype), color = "black", lwd = line_width*0.5, fatten = line_width*2) +
    geom_point(size = 0.65) +
    geom_signif(comparisons = list(
        c("Basal Committed","Subtype-Intermediate Progenitor"),
        c("Classical-NRP","Subtype-Intermediate Progenitor"),
        c("Basal Committed","Classical-NRP")),
        vjust = -0.1,
        step_increase = step_adjustment,
        test = "wilcox.test",
        textsize = 6/.pt,
        margin_top = 0.15) +
    scale_fill_manual(values = colorPalette) +
    # ylim(min(plot$feature),abs(max(plot$feature)-min(plot$feature))*ylim_adjustment) +
    ylab("Count / Total CAFs") +
    default_theme +
    theme(
        legend.position = "none"
    )



############################################################
# Plot iCAF proportions

plot <- clusterProp
plot <- plot[plot$Cluster=="iCAF",]

ggplot(plot, aes(x = Subtype, y = Proportion)) +
    geom_boxplot(outlier.alpha = 0, aes(fill = Subtype), color = "black", lwd = line_width*0.5, fatten = line_width*2) +
    geom_point(size = 0.65) +
    geom_signif(comparisons = list(
        c("Basal Committed","Subtype-Intermediate Progenitor"),
        c("Classical-NRP","Subtype-Intermediate Progenitor"),
        c("Basal Committed","Classical-NRP")),
        vjust = -0.1,
        step_increase = step_adjustment,
        test = "wilcox.test",
        textsize = 6/.pt,
        margin_top = 0.15) +
    scale_fill_manual(values = colorPalette) +
    # ylim(min(plot$feature),abs(max(plot$feature)-min(plot$feature))*ylim_adjustment) +
    ylab("Count / Total CAFs") +
    default_theme +
    theme(
        legend.position = "none"
    )



############################################################
# Plot PSC proportions

plot <- clusterProp
plot <- plot[plot$Cluster=="PSC",]

ggplot(plot, aes(x = Subtype, y = Proportion)) +
    geom_boxplot(outlier.alpha = 0, aes(fill = Subtype), color = "black", lwd = line_width*0.5, fatten = line_width*2) +
    geom_point(size = 0.65) +
    geom_signif(comparisons = list(
        c("Basal Committed","Subtype-Intermediate Progenitor"),
        c("Classical-NRP","Subtype-Intermediate Progenitor"),
        c("Basal Committed","Classical-NRP")),
        vjust = -0.1,
        step_increase = step_adjustment,
        test = "wilcox.test",
        textsize = 6/.pt,
        margin_top = 0.15) +
    scale_fill_manual(values = colorPalette) +
    # ylim(min(plot$feature),abs(max(plot$feature)-min(plot$feature))*ylim_adjustment) +
    ylab("Count / Total CAFs") +
    default_theme +
    theme(
        legend.position = "none"
    )