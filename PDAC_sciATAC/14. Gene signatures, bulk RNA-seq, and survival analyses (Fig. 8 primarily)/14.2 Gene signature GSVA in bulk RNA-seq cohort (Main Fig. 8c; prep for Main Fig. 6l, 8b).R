# Analysis of sci-ATAC-seq gene signatures in the bulk RNA-seq OPTR cohort published in Link JM, et al. (2025) (https://doi.org/10.1038/s43018-024-00881-3).

# Here, we apply the GSVA R package to score each primary tumor in the OPTR cohort (n=218) according to their mRNA expression of the sci-ATAC-seq gene signatures.
# We correlate these scores against scores for published/canonical classical/basal gene signatures as well.



############################################################
# Prepare workspace

library(rtracklayer)
library(dplyr)
library(tximport)
library(DESeq2)
library(GSVA)
library(pheatmap)

library(corrr)

set.seed(1)



############################################################
# Prepare sample space

# 218 primary tumors used in Link JM, et al. (2025)
samples <- list.dirs("~/Dataset/Primary_218/Kallisto_raw_output", recursive = F, full.names = F)

# Read in Kallisto outputs
# Obtain paths to each "abundance.h5" file for each sample
files <- file.path(
    "~/Dataset/Primary_218/Kallisto_raw_output",
    samples,
    "abundance.h5")
names(files) <- samples



############################################################
# Prepare annotation objects

# Import GENCODE genome annotation info
annotation <- import("~/reference_genomes/gencode/gencode.v24.annotation.gtf.gz")
annotation <- as.data.frame(annotation@elementMetadata@listData)

# Prepare tx2gene data for tximport
tx2gene <- data.frame(TXNAME = unique(annotation$transcript_id))
tx2gene <- data.frame(TXNAME = tx2gene[!(is.na(tx2gene$TXNAME)),])
tx2gene$GENEID <- annotation$gene_id[match(tx2gene$TXNAME, annotation$transcript_id)]
tx2gene <- tibble(tx2gene)

# Prepare annotation object for converting ensembl ID to gene name
annotation$gene_id_short <- gsub("\\..*","",annotation$gene_id)  # Remove Ensembl ID end tags
annotation <- annotation[c("gene_id","gene_id_short","gene_type","gene_name")]  # Keep only the info of interest
annotation <- unique(annotation)



############################################################
# tximport

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)



############################################################
# DESeq2

# Set up sample list
coldata <- data.frame(sample = colnames(txi$counts))
row.names(coldata) = coldata$sample

# Generate DESeq dataset from txi object (design is ~ 1 as we are not testing differential expression)
dds <- DESeqDataSetFromTximport(txi, coldata, ~ 1)

# Filter out genes that do not have at a sum of at least 10 reads across all samples 
dds <- dds[rowSums(counts(dds)) >= 10,]



############################################################
# Variance-stabilizing transformation (VST)

# Perform variance-stabilizing transformation
vst <- vst(dds)

# Extract VST counts matrix
matrix <- as.matrix(assay(vst))

# Flip matrix
matrix <- as.data.frame(t(matrix))



############################################################
# Remove intermediate objects

rm(files,
   tx2gene,
   txi,
   coldata,
   vst)



############################################################
# Collect published PDAC subtype gene signatures

signatures <- as.data.frame(readxl::read_xlsx("~/gene_sets/PDAC/PDAC_Gene_Signatures.xlsx", sheet = 2))

features <- list()

for(col in colnames(signatures)){
    features[[col]] = signatures[[col]]
    features[[col]] = features[[col]][is.na(features[[col]])==F]
    features[[col]] = annotation$gene_id_short[match(features[[col]],annotation$gene_name)]
    rm(col)
}



############################################################
# Collect sci-ATAC-seq gene signatures

signatures <- read.table("~/gene_sets/PDAC/sciATAC_Gene_Signatures.txt", header = T, sep = '\t', quote = "", comment.char = "")
signatures <- signatures[,grep("_Ensembl$",colnames(signatures)),drop=F]
colnames(signatures) <- gsub("_Ensembl$","",colnames(signatures))

for(col in colnames(signatures)){
    features[[col]] = signatures[[col]]
    features[[col]] = features[[col]][is.na(features[[col]])==F]
    rm(col)
}



############################################################
# Run GSVA

matrix_for_gsva <- as.matrix(assay(vst(dds)))
rownames(matrix_for_gsva) <- gsub("\\..*","", rownames(matrix_for_gsva))


gsva_results_genesets <- as.data.frame(t(gsva(
    matrix_for_gsva,
    features,
    method = "gsva",
    kcdf = "Gaussian",
    mx.diff = TRUE,
    verbose = TRUE
)))



############################################################
# Export GSVA results

gsva_results_genesets <- gsva_results_genesets
gsva_results_genesets$Public_Specimen_ID <- rownames(gsva_results_genesets)
gsva_results_genesets <- gsva_results_genesets[,c("Public_Specimen_ID",grep("_Signature",colnames(gsva_results_genesets),value = T),grep("_Signature",colnames(gsva_results_genesets),invert = T,value = T))]
gsva_results_genesets <- gsva_results_genesets[order(gsva_results_genesets$SIP_Signature,decreasing = T),]

write.table(source_data, file = "~/BulkRNA_GSVAScores.txt", row.names = F, col.names = T, sep = '\t', quote = F)



############################################################
# Correlation heatmap

correlation <- as.data.frame(corrr::correlate(gsva_results_genesets))
rownames(correlation) <- correlation$term

correlation <- correlation[
    grep("HwangWL2022_MalignantLineage|Moffitt|Raghavan|Bailey",rownames(correlation)),
    grep("SIP|BC|C_NRP",colnames(correlation)),
    drop=F]
correlation <- correlation[grep("_IC_hg38",rownames(correlation), invert = T),]

# Rename signatures to be more human-readable
rownames(correlation) <- gsub("_"," ",rownames(correlation))
rownames(correlation) <- gsub(" MalignantLineage","",rownames(correlation))
rownames(correlation) <- gsub(" Sig.*","",rownames(correlation))
rownames(correlation) <- gsub("PDAC ","",rownames(correlation))
rownames(correlation) <- gsub("NeuralLikeProgenitor","Neural-like Progenitor",rownames(correlation))
rownames(correlation) <- gsub("IC hg19$","Intermediate Co-Exp.",rownames(correlation))



pheatmap(correlation,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward",
         
         scale = "none",
         
         fontsize = 6,
         angle_col = 90,
         cellwidth = 10,
         cellheight = 10,
         
         color = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(100),
         
         treeheight_row = 10,
         treeheight_col = 0,
         
         border_color = NA
)



