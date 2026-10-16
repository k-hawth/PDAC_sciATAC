########################################################################################################################
# PLOT TOPIC x AUC SCORES
########################################################################################################################
# Import AUC signature x topic enrichment scores

auc <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/AUC_Epithelial_Signatures/Epithelial_cisTopic50_cisTopicAUC_cisTopicAUC_SubtypeSignatures_AUCScoresForHeatmap.txt", header = T)

# Rename signatures to be more human-readable
rownames(auc) <- gsub("_"," ",rownames(auc))
rownames(auc) <- gsub(" MalignantLineage","",rownames(auc))
rownames(auc) <- gsub("PDAC ","",rownames(auc))
rownames(auc) <- gsub("NeuralLikeProgenitor","Neural-like Progenitor",rownames(auc))
rownames(auc) <- gsub("IC hg38$","Intermediate Co-Exp.",rownames(auc))
rownames(auc) <- gsub("QuasiMesenchymal","Quasi-Mesenchymal",rownames(auc))

rownames(auc) <- gsub(" vsRest.*","",rownames(auc))
rownames(auc) <- gsub("PancreaticProgenitor","Pancreatic Progenitor",rownames(auc))
rownames(auc) <- gsub(" Like","-like",rownames(auc))

sig <- sub(".*? ","",rownames(auc))
author <- sub("* .*","",rownames(auc))

author <- gsub("([A-Z])(\\d)","\\1 \\2",author)
author <- gsub("([a-z])([A-Z])","\\1 \\2",author)
author <- paste("(",author,")",sep = "")

rownames(auc) <- paste(sig, author)

rownames(auc) <- gsub(" [A-Z]* \\d.*\\)",")", rownames(auc))
colnames(auc) <- gsub("^Topic","T",colnames(auc))



############################################################
# Limit to subtype-associated topics

topics_classicalNRP <- paste("T", c(34,41,30,37,10,3,18,44,40,46,29,45), sep="")
topics_basalCommitted <- paste("T", c(23,47,9,22), sep="")
topics_SIP <- paste("T", c(6,20,42), sep="")


auc <- auc[,colnames(auc) %in% c(topics_classicalNRP, topics_basalCommitted, topics_SIP)]


# Order columns of matrix
auc <- auc[,c(topics_basalCommitted, topics_SIP, topics_classicalNRP)]


############################################################
# ANNOTATION DF: subtype-associated topics

heatmap_annot_col <- data.frame(col = colnames(auc))
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% topics_classicalNRP] <- "Classical-NRP"
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% topics_SIP] <- "Subtype-Intermediate Progenitor"
heatmap_annot_col$Subtype[heatmap_annot_col$col %in% topics_basalCommitted] <- "Basal Committed"
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
# Heatmap

library(pheatmap)

cairo_pdf(width = 4.6, height = 2.17,
          file = "~/cistopicAUC_heatmap.pdf",
          family = 'Arial Unicode MS')

ComplexHeatmap::pheatmap(as.matrix(auc),
                         clustering_distance_rows = "correlation",
                         clustering_method = "ward",
                         
                         cluster_cols = F,
                         
                         scale = "none",
                         
                         annotation_col = heatmap_annot_col,
                         annotation_colors = annot_colors_col,
                         
                         fontsize = 6,
                         angle_col = "90",
                         
                         color = colorRampPalette(ArchR::ArchRPalettes$greenBlue)(10000),
                         
                         annotation_legend = F,
                         cellwidth = 7,
                         cellheight = 7,
                         
                         treeheight_row = 4,
                         treeheight_col = 4,
                         heatmap_legend_param = list(border = FALSE),
                         
                         border_color = NA, 
                         row_names_side="left",
                         row_dend_side = "right",
                         column_names_side = "bottom"
)
dev.off()