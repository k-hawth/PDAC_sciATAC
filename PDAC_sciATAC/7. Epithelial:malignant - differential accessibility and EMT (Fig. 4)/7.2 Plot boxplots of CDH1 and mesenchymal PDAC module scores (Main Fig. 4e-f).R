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
# Calculate Impute Weights using the cisTopic model as the reduced dimensional space
########################################################################################################################

# ArchR's addImputeWeights function by default uses the iterative LSI reduced dimensions
# Below is a modified version of ArchR's addImputeWeights function that accepts other reduced dimensions (here, cisTopic matrix)

addImputeWeights_custom <- function(
    ArchRProj = NULL,
    customDims = cistopic,
    reducedDims = "IterativeLSI", 
    dimsToUse = NULL, 
    scaleDims = NULL,
    corCutOff = 0.75, 
    td = 3,
    ka = 4,
    sampleCells = 5000,
    nRep = 2,
    k = 15,
    epsilon = 1,
    useHdf5 = TRUE,
    randomSuffix = FALSE,
    threads = getArchRThreads(),
    seed = 1,
    verbose = TRUE,
    logFile = createLogFile("addImputeWeights")
){
    
    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
    ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer", "null"))
    ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean", "null"))
    ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
    ArchR:::.validInput(input = td, name = "td", valid = c("integer"))
    ArchR:::.validInput(input = ka, name = "ka", valid = c("integer"))
    ArchR:::.validInput(input = sampleCells, name = "sampleCells", valid = c("integer", "null"))
    ArchR:::.validInput(input = nRep, name = "nRep", valid = c("integer"))
    ArchR:::.validInput(input = k, name = "k", valid = c("integer"))
    ArchR:::.validInput(input = epsilon, name = "epsilon", valid = c("numeric"))
    ArchR:::.validInput(input = useHdf5, name = "useHdf5", valid = c("boolean"))
    ArchR:::.validInput(input = randomSuffix, name = "randomSuffix", valid = c("boolean"))
    ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
    ArchR:::.validInput(input = seed, name = "seed", valid = c("integer"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
    ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
    
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "addImputeWeights Input-Parameters", logFile = logFile)  
    
    #Adapted From
    #https://github.com/dpeerlab/magic/blob/master/R/R/run_magic.R
    
    set.seed(seed)
    
    tstart <- Sys.time()
    ArchR:::.logDiffTime("Computing Impute Weights Using Magic (Cell 2018)", t1 = tstart, verbose = verbose, logFile = logFile)
    
    ## OLD ## Get Reduced Dims
    # matDR <- getReducedDims(ArchRProj, reducedDims = reducedDims, dimsToUse = dimsToUse, corCutOff = corCutOff)
    # N <- nrow(matDR)
    # rn <- rownames(matDR)
    
    ## CUSTOM ## Get Reduced Dims
    matDR <- customDims
    N <- nrow(matDR)
    rn <- rownames(matDR)
    
    if(!is.null(sampleCells)){
        if(sampleCells > nrow(matDR)){
            sampleCells <- NULL
        }
    }
    if(is.null(sampleCells)){
        binSize <- N
        nRep <- 1
    }else{
        cutoffs <- lapply(seq_len(1000), function(x){
            N / x
        }) %>% unlist
        binSize <- min(cutoffs[order(abs(cutoffs - sampleCells))[1]] + 1, N)
    }
    
    if(useHdf5){
        dir.create(file.path(getOutputDirectory(ArchRProj), "ImputeWeights"), showWarnings = FALSE)
        if(randomSuffix){
            weightFiles <- .tempfile("Impute-Weights", tmpdir = file.path(gsub(paste0(getwd(),"/"),"",getOutputDirectory(ArchRProj)), "ImputeWeights"))
            weightFiles <- paste0(weightFiles, "-Rep-", seq_len(nRep))
        }else{
            weightFiles <- file.path(getOutputDirectory(ArchRProj), "ImputeWeights", paste0("Impute-Weights-Rep-", seq_len(nRep)))
        }
        o <- suppressWarnings(file.remove(weightFiles))
    }
    
    weightList <- ArchR:::.safelapply(seq_len(nRep), function(y){
        
        ArchR:::.logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s)", y, nRep), t1 = tstart, verbose = FALSE, logFile = logFile)
        
        if(!is.null(sampleCells)){
            idx <- sample(seq_len(nrow(matDR)), nrow(matDR))
            blocks <- split(rownames(matDR)[idx], ceiling(seq_along(idx)/binSize))
        }else{
            blocks <- list(rownames(matDR)) 
        }
        
        if(useHdf5){
            weightFile <- weightFiles[y]
            o <- h5createFile(weightFile)
        }
        
        blockList <- lapply(seq_along(blocks), function(x){
            
            if(x %% 10 == 0){
                ArchR:::.logDiffTime(sprintf("Computing Partial Diffusion Matrix with Magic (%s of %s, Iteration %s of %s)", y, nRep, x, length(blocks)), 
                                     t1 = tstart, verbose = FALSE, logFile = logFile)
            }
            
            ix <- blocks[[x]]
            Nx <- length(ix)
            
            #Compute KNN
            knnObj <- nabor::knn(data = matDR[ix,], query = matDR[ix, ], k = k)
            knnIdx <- knnObj$nn.idx
            knnDist <- knnObj$nn.dists
            rm(knnObj)
            
            if(ka > 0){
                knnDist <- knnDist / knnDist[,ka]
            }
            
            if(epsilon > 0){
                W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=c(knnDist), dims = c(Nx, Nx))
            } else {
                W <- Matrix::sparseMatrix(rep(seq_len(Nx), k), c(knnIdx), x=1, dims = c(Nx, Nx)) # unweighted kNN graph
            }
            W <- W + Matrix::t(W)
            
            #Compute Kernel
            if(epsilon > 0){
                W@x <- exp(-(W@x / epsilon^2))
            }
            
            #Markov normalization
            W <- W / Matrix::rowSums(W) 
            
            #Initialize Matrix
            Wt <- W
            
            #Computing Diffusion Matrix
            for(i in seq_len(td)){
                Wt <- Wt %*% W
            }
            rownames(Wt) <- ix
            colnames(Wt) <- ix
            
            rm(knnIdx)
            rm(knnDist)
            rm(W)
            gc()
            
            if(useHdf5){
                o <- ArchR:::.suppressAll(h5createGroup(file = weightFile, paste0("block", x)))      
                o <- ArchR:::.suppressAll(h5write(obj = ix, file = weightFile, name = paste0("block", x, "/Names"), level = 0))
                o <- ArchR:::.suppressAll(h5write(obj = as.matrix(Wt), file = weightFile, name = paste0("block", x, "/Weights"), level = 0))
                return(weightFile)
            }else{
                Wt
            }
            
        }) %>% SimpleList
        
        if(useHdf5){
            return(weightFile)
        }else{
            names(blockList) <- paste0("b",seq_along(blockList))
            blockList
        }
    }, threads = threads) %>% SimpleList
    names(weightList) <- paste0("w",seq_along(weightList))
    
    ArchR:::.logDiffTime(sprintf("Completed Getting Magic Weights! Object size - %s.", round(object.size(weightList) / 10^9, 3)), 
                         t1 = tstart, verbose = FALSE, logFile = logFile)
    
    ArchRProj@imputeWeights <- SimpleList(
        Weights = weightList, 
        Params = 
            list(
                reducedDims = reducedDims, 
                td = td, 
                k = k, 
                ka = ka,
                epsilon = epsilon
            )
    )
    
    ArchRProj
    
}



# Load scaled cell-by-topic score matrix
cistopic <- read.table("~/cisTopic/AllSamples_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS3_AllSamplesTSSMixednFrags1kFRIS0.4/Epithelial/cisTopic100_Alpha50/Perplexity_topics/cisTopic100Alpha50_EpithelialAllSamplesTSSMixnFrags1kFRIS0.4_AllSamplesMACS3_Topics50_Matrix_Scaled.matrix")

# Transpose: rows = cells, cols = topics
cistopic <- t(cistopic)

# Add ArchRProjct prefix
rownames(cistopic) <- paste("sears_tumors_hg38#",rownames(cistopic),sep="")

# Reorder cells (rows) to match in ArchRProject
cistopic <- cistopic[proj$cellNames,]

# Add impute weights using cisTopic matrix
proj <- addImputeWeights_custom(ArchRProj = proj,
                                customDims = cistopic,
                                sampleCells = length(proj$cellNames),
                                scaleDims = FALSE)



##################################################################
# Module score matrix

modules <- read.table("~/ArchR/Epithelial_AllSamples_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_cisTopic50/Plots/cisTopic50_EpithelialPrimary_TSSMixed_nFrags1k_FRIS0.4_PeaksMACS3_CellIdentities_ModuleScores.txt", header = T, sep = '\t', quote = "", comment.char = "")

modules_imputed <- modules
rownames(modules_imputed) <- modules_imputed$cellNames
modules_imputed <- modules_imputed[,grep("^Module",colnames(modules_imputed)),drop=F]


# Impute (columns = cells, rows = features)
modules_imputed <- imputeMatrix(
    mat = t(modules_imputed),
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = FALSE,
    logFile = createLogFile("imputeMatrix")
)

modules_imputed <- as.data.frame(t(modules_imputed))



############################################################
# GeneScoreMatrix

matrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix", logFile = NULL, binarize = F)
gene_index <- as.data.frame(matrix@elementMetadata@listData)
gene_index$index <- 1:length(rownames(gene_index))
matrix <- assay(matrix)



gene = "CDH1"
matrix_gene <- data.frame(row.names = proj$cellNames,
                          gene_imputed = matrix[as.numeric(rownames(gene_index[gene_index$name==gene,])),proj$cellNames])
matrix_gene <- as.matrix(matrix_gene)
matrix_imputed <- imputeMatrix(
    mat = t(matrix_gene),
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = FALSE,
    logFile = createLogFile("imputeMatrix")
)
matrix_imputed <- as.data.frame(matrix_imputed)
colnames(matrix_imputed) <- gene



############################################################
# Prepare plot data frame

groups_to_plot <- c("A/D",
                    "C1","C10","C12","C17",
                    "C5","C6","C22","C3","C14",
                    "C19","C13","C9")

plot <- data.frame(cell = proj$cellNames,
                   sample = proj$Sample_ID,
                   group = proj$cistopic50_pg_k15)

# Rename clusters
plot$group <- gsub("^pg_","C",plot$group)
plot$group <- gsub("^C11$","A/D",plot$group)

# Keep relevant groups
plot <- plot[plot$group %in% groups_to_plot,]

# Remove metastatic cells in acinar/ductal cluster
plot <- plot[(plot$group=="A/D" & plot$sample %in% grep("-M",plot$sample,value = T))==F,]

# Subtype groups
colorPalette[["A/D"]] <- "#174a7e" 
colorPalette[["Subtype-Intermediate Progenitor"]] <- "#e1468a"
colorPalette[["Basal Committed"]] <- "#bd3531"
colorPalette[["Classical-NRP"]] <- "#f79747"
plot$subtype[plot$group %in% c("A/D")] <- "A/D"
plot$subtype[plot$group %in% c("C1","C5","C6","C10","C12","C22")] <- "Subtype-Intermediate Progenitor"
plot$subtype[plot$group %in% c("C3","C13","C17")] <- "Basal Committed"
plot$subtype[plot$group %in% c("C9","C14","C19")] <- "Classical-NRP"

# EMT groups
plot$emt[plot$group %in% c("A/D","C1","C10","C12","C17")] <- "Epithelial"
plot$emt[plot$group %in% c("C5","C6","C22","C3","C14")] <- "Partial EMT"
plot$emt[plot$group %in% c("C19","C13","C9")] <- "Mesenchymal"

# Order of groups
plot$group <- factor(plot$group, levels = groups_to_plot)
plot$emt <- factor(plot$emt, levels = c("Epithelial","Partial EMT","Mesenchymal"))




############################################################
# Plot data

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




# Mesenchymal module score
plot$feature <- modules_imputed$Module.HwangWL2022_MalignantLineage_Mesenchymal[match(plot$cell, rownames(modules_imputed))]

pdf(height = 1.1, width = 3.24, file = "EMT_Violin_MesenchymalModuleScore.pdf")
ggplot(plot, aes(x = group, y = feature, color = subtype)) +
    geom_boxplot(outlier.size = point_size) +
    scale_color_manual(values = colorPalette) +
    facet_grid(cols = vars(plot$emt), scales = "free", space = "free") +
    default_theme +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_text(size = 6, color = "black", vjust = -1.8),
        
        plot.title = element_blank()
    ) +
    ylab("Score (Imputed)")
dev.off()



# CDH1
plot$feature <- matrix_imputed[match(plot$cell,rownames(matrix_imputed)),gene]

pdf(height = 1.1, width = 3.19, file = "EMT_Violin_CDH1GeneScore.pdf")
ggplot(plot, aes(x = group, y = feature, color = subtype)) +
    geom_boxplot(outlier.size = point_size) +
    scale_color_manual(values = colorPalette) +
    facet_grid(cols = vars(plot$emt), scales = "free", space = "free") +
    default_theme +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        
        strip.background = element_blank(),
        strip.text = element_text(size = 6, color = "black", vjust = -1.8),
        
        plot.title = element_blank()
    ) +
    ylab("Score (Imputed)")
dev.off()
