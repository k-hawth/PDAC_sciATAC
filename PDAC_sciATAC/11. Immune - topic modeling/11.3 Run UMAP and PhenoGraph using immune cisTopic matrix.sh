# Run UMAP using scaled cisTopic cell-by-topic score matrix
scitools matrix-umap cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Matrix_Scaled.matrix

# Run PhenoGraph using scaled cisTopic cell-by-topic score matrix
scitools matrix-pg -k 15 cisTopic_alpha50_ImmuneTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics43_Matrix_Scaled.matrix