library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



colorPalette <- c("Basal Committed" = "#bd3531",
                  "Classical-NRP" = "#f79747",
                  "Subtype-Intermediate Progenitor" = "#e1468a",
                  "Malignant PAV" = "#77a00a",
                  "Acinar/Ductal" = "#174a7e",
                  "Metaplastic NRP" = "#7c177e",
                  "Islet" = "#c48eed",
                  "Intestinal Epithelial" = "#177e53",
                  "Other PAV" = "#77a00a",
                  "Metaplastic Gastric" = "#ffd898",
                  "Pulmonary Epithelial" = "#729292",                 
                  "Non.Highlighted" = "#D3D3D3"
)



########################################################################################################################
# Prepare cisTopic topic scores for linear mixed effects analysis
# In this section, each sample will have their epithelial/malignant cells separated by cell type.
# The mean topic score of each cell type subpopulation will then be calculated for the sample.
########################################################################################################################

# Load cisTopic data
cistopic <- as.data.frame(read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix.matrix"))

# Scale matrix
cistopic <- scale(cistopic, center = TRUE, scale = TRUE)
rownames(cistopic) <- paste(rownames(cistopic), "_Scaled", sep = "")

# Transpose
cistopic <- as.data.frame(t(cistopic))

rownames(cistopic) <- paste("sears_tumors_hg38",rownames(cistopic), sep = "#")



############################################################
# Function to collect mean topic value of each cell type within each sample

sampleGroup_Mean <- function(
    sample_group_cell_annot,
    cell_x_value_matrix
){
    
    final_output <- sample_group_cell_annot
    final_output <- final_output %>%
        group_by(sample) %>%
        mutate(count_total = n()) %>%
        group_by(sample,group) %>%
        mutate(count_group = n()) %>%
        select(-cell) %>%
        unique()
    
    for(feature in colnames(cell_x_value_matrix)){
        output <- sample_group_cell_annot
        output$feature <- cell_x_value_matrix[match(output$cell, rownames(cell_x_value_matrix)),feature]
        
        output <- output %>%
            group_by(sample,group) %>%
            mutate(mean_feature_score = mean(feature)) %>%
            select(-cell) %>%
            unique()
        
        final_output[[paste(feature,"_Mean",sep="")]] <- output$mean_feature_score[match(paste(final_output$sample,final_output$group),paste(output$sample,output$group))]
    }
    
    
    
    return(final_output)
    
}



############################################################
# Run function and filter for cell types

df <- data.frame(sample = proj$Sample_ID,
                 group = proj$celltype,
                 cell = proj$cellNames)

df <- sampleGroup_Mean(df, cistopic)


# Discard tumor subpopulation if less than min cell threshold
df <- df[df$count_group>=10,]

# Remove unwanted cell types
df <- df[grep("Other|Pulmonary",df$group, invert = T),]

# Discard metastases
df <- df[grep("-T",df$sample),]

# Convert groups to factors
df$group <- factor(df$group)






########################################################################################################################
# Plot topics with all cell types included
########################################################################################################################


############################################################
# ggplot theme

# Theme Settings
text_settings <- element_text(colour = "black", size = 6)
line_settings <- element_line(linewidth = 0.5)

# Default Theme
default_theme <- theme(
    axis.title.x = text_settings,
    axis.title.y = text_settings,
    
    axis.ticks.x = line_settings,
    axis.ticks.y = line_settings,
    
    axis.text.x = text_settings,
    axis.text.y = text_settings,
    
    #axis.line.x = element_line(linewidth = 0.5),
    #axis.line.y = element_line(linewidth = 0.5),
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    
    legend.text = text_settings,
    legend.title = text_settings
)





plot <- df

# Convert groups to factors
plot$group <- factor(plot$group, levels = c("Basal Committed","Classical-NRP", "Subtype-Intermediate Progenitor","Acinar/Ductal","Islet","Metaplastic NRP","Metaplastic Gastric","Intestinal Epithelial","Malignant PAV"))



pdf(width = 4.3, height = 1.7, file = "topic_scaled_avg.pdf")
for(unique_topic in grep("^Topic",colnames(plot),value = T)){
    p <- ggplot(plot, aes(x = group, y = !!sym(unique_topic))) +
        geom_boxplot(outlier.alpha = 0, aes(fill = group)) +
        geom_point(size = 1) +
        default_theme +
        theme(
            axis.title.x = element_blank(),
            plot.title = element_text(size = 6,vjust = -2),
            legend.position = "none"
        ) +
        scale_x_discrete(labels = scales::label_wrap(10)) +
        scale_fill_manual(values = colorPalette) +
        ylab("Avg Scaled Topic Score\nof Cell Type Per Sample") +
        ggtitle(gsub(" Scaled Mean","", gsub("_"," ",unique_topic)))
    print(p)
    rm(unique_topic)
}
dev.off()
