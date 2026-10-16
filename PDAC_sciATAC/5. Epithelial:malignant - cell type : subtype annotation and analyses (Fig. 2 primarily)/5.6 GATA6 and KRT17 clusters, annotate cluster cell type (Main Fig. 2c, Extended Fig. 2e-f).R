library(ArchR)
set.seed(1)
addArchRGenome("hg38")
library(parallel)
library(dplyr)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50")



############################################################
# Extract gene scores from GeneScoreMatrix

counts <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", logFile = NULL, binarize = FALSE)
gene_index <- as.data.frame(counts@elementMetadata@listData)
gene_index$idx <- as.numeric(rownames(gene_index))
counts <- assay(counts)

df <- data.frame(cell = proj$cellNames,
                 pg = proj$cistopic50_pg_k15,
                 gata6 = counts[gene_index$idx[gene_index$name=="GATA6"],proj$cellNames],
                 krt17 = counts[gene_index$idx[gene_index$name=="KRT17"],proj$cellNames])



############################################################
# Compare GATA6 gene scores between acinar ductal cluster and each malignant cluster

p1 <- data.frame(x = df$pg,
                 y = df$gata6)

# All malignant clusters
overall <- data.frame(group = c("pg_3","pg_17","pg_1","pg_13","pg_14","pg_19","pg_9","pg_22","pg_6","pg_10","pg_5","pg_12","pg_4","pg_16"))

for(group in overall$group){
    overall$wilcox[overall$group==group] <- wilcox.test(
        p1$y[p1$x==group],
        p1$y[p1$x=="pg_11"],paired = F
    )[[3]]
    
    overall$avg_log2FC[overall$group==group] <- log2(mean(p1$y[p1$x==group]) / mean(p1$y[p1$x=="pg_11"]))
}

overall_gata6 <- overall



############################################################
# Compare KRT17 gene scores between acinar ductal cluster and each malignant cluster

p1 <- data.frame(x = df$pg,
                 y = df$krt17)

# All malignant clusters
overall <- data.frame(group = c("pg_3","pg_17","pg_1","pg_13","pg_14","pg_19","pg_9","pg_22","pg_6","pg_10","pg_5","pg_12","pg_4","pg_16"))

for(group in overall$group){
    overall$wilcox[overall$group==group] <- wilcox.test(
        p1$y[p1$x==group],
        p1$y[p1$x=="pg_11"],paired = F
    )[[3]]
    
    overall$avg_log2FC[overall$group==group] <- log2(mean(p1$y[p1$x==group]) / mean(p1$y[p1$x=="pg_11"]))
}

overall_krt17 <- overall



############################################################
# Apply BH FDR correction to all p-values

overall_gata6$gene = "GATA6"
overall_krt17$gene = "KRT17"

overall <- rbind(overall_gata6,overall_krt17)
overall$padj <- p.adjust(overall$wilcox,method = "fdr")

overall_krt17 <- overall[overall$gene=="KRT17",]
overall_gata6 <- overall[overall$gene=="GATA6",]



############################################################
# Classify clusters as KRT17+, GATA6+, or KRT17+ GATA6+

basal <- overall_krt17$group[overall_krt17$padj<=0.05 & overall_krt17$avg_log2FC >= 0.5]
classical <- overall_gata6$group[!( overall_gata6$padj<=0.05 & overall_gata6$avg_log2FC<= (-0.5) )]
hybrid <- basal[basal %in% classical==T]

classical <- classical[classical %in% hybrid==F]
basal <- basal[basal %in% hybrid==F]


overall$Subtype[overall$group %in% basal] <- "KRT17+"
overall$Subtype[overall$group %in% classical] <- "GATA6+"
overall$Subtype[overall$group %in% hybrid] <- "KRT17+ GATA6+"



############################################################
# Plot data
 
overall$point_color[overall$group %in% basal] = "#bd3531"
overall$point_color[overall$group %in% classical] = "#f79747"
overall$point_color[is.na(overall$point_color)] = "#e1468a"

overall$text_color <- ifelse(overall$group %in% basal, "gray98", "black")
overall$group <- gsub("pg_","",overall$group)
overall$group_text <- overall$group


ggplot(overall[overall$gene=="GATA6",], aes(y = -log10(padj), x = avg_log2FC)) +
    annotate("rect", xmin = -0.5, xmax = -Inf, ymin = -log10(0.05), ymax = Inf, alpha = .2) +
    geom_hline(yintercept = c(-log10(0.05)),linetype = 2) +
    geom_vline(xintercept = c(-0.5,0.5),linetype = 2) +
    geom_point(aes(color = point_color), size = 3) +
    theme_bw() +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.text.x = element_text(size = 6, color = "black"),
          text = element_text(size = 6, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          aspect.ratio = 1
    ) +
    xlim(-1.802738110, 1.802738110) + 
    geom_text(aes(label = group_text, color = text_color), nudge_y = 0,size = 2) +
    scale_color_identity() +
    xlab("Average Log2 Fold Change") + ylab("Adjusted p-value (-log10)")




ggplot(overall[overall$gene=="KRT17",], aes(y = -log10(padj), x = avg_log2FC)) +
    annotate("rect", xmin = 0.5, xmax = -Inf, ymin =-log10(0.05), ymax = -Inf, alpha = .2) +
    annotate("rect", xmin = 0.5, xmax = -Inf, ymin = -log10(0.05), ymax = Inf, alpha = .2) +
    annotate("rect", xmin = 0.5, xmax = Inf, ymin = -log10(0.05), ymax = -Inf, alpha = .2) +
    geom_hline(yintercept = c(-log10(0.05)),linetype = 2) +
    geom_vline(xintercept = c(-0.5,0.5),linetype = 2) +
    geom_point(aes(color = point_color), size = 3) +
    theme_bw() +
    theme(axis.text = element_text(size = 6, color = "black"),
          axis.text.x = element_text(size = 6, color = "black"),
          text = element_text(size = 6, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          aspect.ratio = 1
    ) +
    xlim(-2.84191530, 2.84191530) + 
    geom_text(aes(label = group_text, color = text_color), nudge_y = 0,size = 2) +
    scale_color_identity() +
    xlab("Average Log2 Fold Change") + ylab("Adjusted p-value (-log10)")

ggsave(filename = "krt17_volcanoPlot_shaded_manualLabel_PAVGreen.pdf",device = "pdf", width = 2, height = 2, units = "in")





########################################################################################################################
# Annotate epithelial/malignant cluster cell types
########################################################################################################################

annot <- data.frame(cell = proj$cellNames,
                    sample = proj$Sample_ID,
                    pg = proj$cistopic50_pg_k15)

annot$celltype <- NA
annot$celltype[annot$pg %in% c("pg_1","pg_5","pg_6","pg_10","pg_12","pg_22")] = "Subtype-Intermediate Progenitor"
annot$celltype[annot$pg %in% c("pg_2","pg_8")] = "Intestinal Epithelial"
annot$celltype[annot$pg %in% c("pg_3","pg_13","pg_17")] = "Basal Committed"
annot$celltype[annot$pg %in% c("pg_4","pg_16")] = "Malignant PAV"
annot$celltype[annot$pg %in% c("pg_7")] = "Metaplastic Gastric"
annot$celltype[annot$pg %in% c("pg_9","pg_14","pg_19")] = "Classical-NRP"
annot$celltype[annot$pg %in% c("pg_11")] = "Acinar/Ductal"
annot$celltype[annot$pg %in% c("pg_15")] = "Islet"
annot$celltype[annot$pg %in% c("pg_18","pg_20")] = "Metaplastic NRP"
annot$celltype[annot$pg %in% c("pg_21")] = "Pulmonary Epithelial"
annot$celltype[annot$pg %in% c("pg_23")] = "Other PAV"

# Separate pulmonary epithelial cells types not expected to be found in primary tumors
annot$celltype[annot$celltype %in% c("Pulmonary Epithelial") & annot$sample %in% unique(annot$sample[grep("-T",annot$sample)])] <- "Other Epithelial (Primary Tumors)"

# Arrange nonsense epithelial from metastases into a separate group
annot$celltype[annot$celltype %in% c("Metaplastic NRP","Intestinal Epithelial","Islet","Acinar/Ductal","Metaplastic Gastric") & annot$sample %in% unique(annot$sample[grep("-M",annot$sample)])] = "Other Epithelial (Metastases)"

# Add to ArchRProject
proj <- addCellColData(proj, cells = annot$cell, data = annot$celltype, force = T, name = "celltype")



############################################################
# Plot UMAP colored by cell types

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

plotEmbedding(proj, name = "celltype", embedding = "cisTopic_50", pal = colorPalette)



############################################################
# Save ArchRProject

proj <- saveArchRProject(proj)