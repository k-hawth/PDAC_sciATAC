# In this section, we perform two survival analyses:
# 1) Calculating univariate hazard ratios for GSVA scores corresponding to PDAC and CAF subtype signatures
# 2) Kaplan-Meier curves and log-rank test comparing survival of SIP tumors in the sci-ATAC-seq cohort to basal committed and classical-NRP tumors.


library(survival)
library(survminer)
library(ggplot2)



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



########################################################################################################################
# Load data
########################################################################################################################


############################################################
# Load per-sample GSVA scores

gsva <- read.table("~/BulkRNA_GSVAScores.txt", header = T, sep = '\t', quote = "", comment.char = "")



############################################################
# Load survival information for the bulk RNA-seq cohort from Link JM, et al (2025) (https://doi.org/10.1038/s43018-024-00881-3)

# Supplemental_Dataset_1 obtained from https://github.com/engjen/Liver_Lung_PDAC
bulkRNA_cohort_info <- readxl::read_xlsx("~/Downloads/Supplemental_Dataset_1.xlsx", sheet = 2)

# Filter for primary tumors only
bulkRNA_cohort_info <- bulkRNA_cohort_info[bulkRNA_cohort_info$Tumor_Type=="Primary",]



############################################################
# Prepare survival data for bulk RNA-seq survival analysis

survival <- data.frame(sample = bulkRNA_cohort_info$Public_Specimen_ID,
                       survival = bulkRNA_cohort_info$Days_from_Diagnosis_to_FU,
                       censor = bulkRNA_cohort_info$Vital_Status_at_FU
)

# Censoring
survival$censor[survival$censor=="Alive"] <- "1"
survival$censor[survival$censor=="Dead"] <- "2"
survival$censor <- as.numeric(survival$censor)

# Remove -T2 samples to prevent double-counting the same patient
survival <- survival[grep("-T2",survival$sample, invert = T),]



########################################################################################################################
# Calculate univariate hazard ratio using GSVA scores for each signature
########################################################################################################################


############################################################
# PDAC gene signatures

feature_list <- rev(c("SIP_Signature",
                      "MoffittRA2015_PDAC_Classical",
                      "RaghavanS2021_scClassical",
                      "C_NRP_Signature",
                      "HwangWL2022_MalignantLineage_NeuralLikeProgenitor",
                      "HwangWL2022_MalignantLineage_Neuroendocrine_Like",
                      "BC_Signature",
                      "MoffittRA2015_PDAC_Basal_Like",
                      "RaghavanS2021_scBasal"))


plot <- data.frame()
for(feature in feature_list){
    
    cox <- survival
    cox$feature <- gsva[match(cox$sample, gsva$Public_Specimen_ID), feature]
    cox <- coxph(Surv(survival, censor) ~ feature, data = cox)
    
    pval <- summary(cox)$coefficients[,"Pr(>|z|)"]
    coef <- summary(cox)$coefficients[,"exp(coef)"]
    upper <- summary(cox)$conf.int[,"upper .95"]
    lower <- summary(cox)$conf.int[,"lower .95"]
    
    plot[length(rownames(plot))+1,"Feature"] <- feature
    plot[length(rownames(plot)),"HR"] <- coef
    plot[length(rownames(plot)),"pval"] <- pval
    plot[length(rownames(plot)),"Upper_CI"] <- upper
    plot[length(rownames(plot)),"Lower_CI"] <- lower
}


# FDR column
plot$FDR_number <- p.adjust(plot$pval, method = "fdr")
plot$FDR <- plot$FDR_number
plot$FDR <- signif(plot$FDR, digits = 3)
plot$FDR <- paste("FDR=",plot$FDR,sep="")

plot$FDR_formatC <- plot$FDR_number
plot$FDR_formatC[grep("e",plot$FDR_formatC)] <- NA
plot$FDR_formatC <- formatC(signif(plot$FDR_formatC, digits = 3), format = "fg", flag = "#", digits = 3)
plot$FDR_formatC <- paste("FDR=",plot$FDR_formatC,sep="")
plot$FDR_formatC[grep("NA", plot$FDR_formatC)] <- NA
plot$FDR[is.na(plot$FDR_formatC)==F] <- plot$FDR_formatC[is.na(plot$FDR_formatC)==F]

# Color of point
plot$point_color <- "black"
plot$point_color[p.adjust(plot$pval, method = "fdr")<=0.05] = "red"

plot$Feature <- factor(plot$Feature, levels = unique(feature_list))

p_cox <- ggplot(plot, aes(x = log2(HR), y = Feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    geom_linerange(aes(xmin=log2(Lower_CI), xmax=log2(Upper_CI))) +
    geom_point(shape = 15, size = 2, color = plot$point_color) +
    geom_text(aes(label = FDR), nudge_y = 0.30, size = 6/.pt) +
    default_theme +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = text_settings
    ) +
    xlab("log2(Univariate HR)") +
    xlim(-(max(abs(log2(c(plot$Upper_CI,plot$Lower_CI))))), max(abs(log2(c(plot$Upper_CI,plot$Lower_CI)))))
p_cox



############################################################
# CAF gene signatures

feature_list <- rev(c("NGFR_myCAF_Signature_Topic22",
                      "LRRC15_myCAF_Signature_Topic15",
                      "iCAF_Signature_Topic11"))


plot <- data.frame()
for(feature in feature_list){
    
    cox <- survival
    cox$feature <- gsva[match(cox$sample, gsva$Public_Specimen_ID), feature]
    cox <- coxph(Surv(survival, censor) ~ feature, data = cox)
    
    pval <- summary(cox)$coefficients[,"Pr(>|z|)"]
    coef <- summary(cox)$coefficients[,"exp(coef)"]
    upper <- summary(cox)$conf.int[,"upper .95"]
    lower <- summary(cox)$conf.int[,"lower .95"]
    
    plot[length(rownames(plot))+1,"Feature"] <- feature
    plot[length(rownames(plot)),"HR"] <- coef
    plot[length(rownames(plot)),"pval"] <- pval
    plot[length(rownames(plot)),"Upper_CI"] <- upper
    plot[length(rownames(plot)),"Lower_CI"] <- lower
}


# FDR column
plot$FDR_number <- p.adjust(plot$pval, method = "fdr")
plot$FDR <- plot$FDR_number
plot$FDR <- signif(plot$FDR, digits = 3)
plot$FDR <- paste("FDR=",plot$FDR,sep="")

plot$FDR_formatC <- plot$FDR_number
plot$FDR_formatC[grep("e",plot$FDR_formatC)] <- NA
plot$FDR_formatC <- formatC(signif(plot$FDR_formatC, digits = 3), format = "fg", flag = "#", digits = 3)
plot$FDR_formatC <- paste("FDR=",plot$FDR_formatC,sep="")
plot$FDR_formatC[grep("NA", plot$FDR_formatC)] <- NA
plot$FDR[is.na(plot$FDR_formatC)==F] <- plot$FDR_formatC[is.na(plot$FDR_formatC)==F]

# Color of point
plot$point_color <- "black"
plot$point_color[p.adjust(plot$pval, method = "fdr")<=0.05] = "red"

plot$Feature <- factor(plot$Feature, levels = unique(feature_list))

p_cox <- ggplot(plot, aes(x = log2(HR), y = Feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    geom_linerange(aes(xmin=log2(Lower_CI), xmax=log2(Upper_CI))) +
    geom_point(shape = 15, size = 2, color = plot$point_color) +
    geom_text(aes(label = FDR), nudge_y = 0.30, size = 6/.pt) +
    default_theme +
    theme(
        axis.title.y = element_blank(),
        axis.title.x = text_settings
    ) +
    xlab("log2(Univariate HR)") +
    xlim(-(max(abs(log2(c(plot$Upper_CI,plot$Lower_CI))))), max(abs(log2(c(plot$Upper_CI,plot$Lower_CI)))))
p_cox



########################################################################################################################
# Kaplan-Meier curve of sci-ATAC-seq cohort
########################################################################################################################


############################################################
# Load metadata for the sci-ATAC-seq cohort

# Sample info on the sci-ATAC-seq cohort in Supp. Table 1
sciATAC_cohort_info <- readxl::read_xlsx("~/PDAC_sciATAC/Supplementary_Tables/Supplementary_Table_1.xlsx", sheet = 1)


# Prepare survival info for primary PDAC tumors with an epigenetic subtype assigned
survival_sciATAC <- sciATAC_cohort_info
survival_sciATAC <- survival_sciATAC[survival_sciATAC$Specimen_Type=="Primary Tumor" & survival_sciATAC$sciATAC_Subtype %in% c("Basal Committed","Classical-NRP", "Subtype Intermediate Progenitor"),]

survival_sciATAC <- data.frame(sample = survival_sciATAC$Public_Specimen_ID,
                               survival = survival_sciATAC$Days_from_Diagnosis_to_FU,
                               censor = survival_sciATAC$Vital_Status_at_FU,
                               group = survival_sciATAC$sciATAC_Subtype
)

# Censoring
survival_sciATAC$censor[survival_sciATAC$censor=="Alive"] <- "1"
survival_sciATAC$censor[survival_sciATAC$censor=="Dead"] <- "2"
survival_sciATAC$censor <- as.numeric(survival_sciATAC$censor)



############################################################
# Perform log-rank test comparing SIP vs. classical-NRP, and SIP vs. basal committed

reference_group <- "Subtype Intermediate Progenitor"
survival_sciATAC$group <- factor(survival_sciATAC$group, levels = unique(survival_sciATAC$group))
survival_sciATAC$group <- relevel(survival_sciATAC$group, ref = reference_group)


# Perform log-rank tests
pairwise_results <- pairwise_survdiff(Surv(survival, censor) ~ group, 
                                      data = survival_sciATAC,
                                      p.adjust.method = "none")
pairwise_results <- as.data.frame(pairwise_results$p.value)


# Get FDR-corrected p-values 
p.adjust(
    c(pairwise_results["Classical-NRP","Subtype Intermediate Progenitor"],
      pairwise_results["Basal Committed","Subtype Intermediate Progenitor"]),
    method = "fdr")


# Visualize as KM plot
p_km <- ggsurvplot(survfit(Surv(survival, censor) ~ group, data = survival_sciATAC),
                   pval = F,
                   conf.int = F,
                   palette = c("#e1468a","#bd3531","#f79747"),
                   risk.table = F,
                   ggtheme = 
                       theme_bw() +
                       default_theme,
                   xlim = c(0,1750))

p_km