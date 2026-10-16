library(dplyr)



############################################################
# Calculate CNV counts (gains + losses) per cell

matrix <- read.table("~/epiAneufinder/AllSamples_TSSMixed_nFrags1k_FRIS0.4/epiAneufinder/results_table.tsv")
colnames(matrix) <- gsub("^cell\\.","sears_tumors_hg38#", colnames(matrix))
rownames(matrix) <- paste(matrix$seq, matrix$start, matrix$end, sep = "_")

matrix <- as.data.frame(t(matrix))
matrix <- matrix[4:length(rownames(matrix)),]

matrix_backup <- matrix

matrix <- sapply(matrix, as.integer)
matrix <- as.data.frame(matrix)

rownames(matrix) <- rownames(matrix_backup)

matrix_backup <- matrix

matrix <- t(matrix)



matrix_gain <- matrix
matrix_gain[matrix_gain==0] <- 0
matrix_gain[matrix_gain==1] <- 0
matrix_gain[matrix_gain==2] <- 1

matrix_gain <- as.data.frame(colSums(matrix_gain))
colnames(matrix_gain) <- "Gain"



matrix_loss <- matrix
matrix_loss[matrix_loss==0] <- 3
matrix_loss[matrix_loss==1] <- 0
matrix_loss[matrix_loss==2] <- 0
matrix_loss[matrix_loss==3] <- 1

matrix_loss <- as.data.frame(colSums(matrix_loss))
colnames(matrix_loss) <- "Loss"



matrix_normal <- matrix
matrix_normal[matrix_normal==0] <- 0
matrix_normal[matrix_normal==1] <- 1
matrix_normal[matrix_normal==2] <- 0

matrix_normal <- as.data.frame(colSums(matrix_normal))
colnames(matrix_normal) <- "Normal"



output <- cbind(matrix_normal,matrix_gain,matrix_loss)
output$Cell <- rownames(output)
output <- output[,c("Cell","Normal","Gain","Loss")]

write.table(output, file = "epiAneufinder_SummaryOfEvents_PerCell.txt", row.names = F, col.names = T, sep = '\t', quote = F)
