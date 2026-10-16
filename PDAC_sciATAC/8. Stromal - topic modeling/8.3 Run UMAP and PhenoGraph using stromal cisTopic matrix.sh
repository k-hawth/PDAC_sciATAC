# Run UMAP using scaled cisTopic cell-by-topic score matrix
scitools matrix-umap cisTopic_alpha50_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics27_Matrix_Scaled.matrix

# Run PhenoGraph using scaled cisTopic cell-by-topic score matrix
scitools matrix-pg -k 50 cisTopic_alpha50_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics27_Matrix_Scaled.matrix