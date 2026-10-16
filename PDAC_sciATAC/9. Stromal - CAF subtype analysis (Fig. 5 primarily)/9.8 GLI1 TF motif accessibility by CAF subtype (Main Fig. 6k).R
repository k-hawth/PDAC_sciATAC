library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/")



########################################################################################################################
# Calculate IMPUTED GLI1 TF motif scores

# In this section, ArchR's implementation of MAGIC will be applied to the TF motif scores
# However, imputation weights will be calculated using the cisTopic model as the reduced dimensional space
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
cistopic <- read.table("~/cisTopic/Primary_TSSMixed_nFrags1k_FRIS0.4/peakset_MACS2ArchR_PrimaryTSSMixnFrags1kFRIS0.4/Fibroblast/cisTopic_50/Derivative_topics/cisTopic_alpha50_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_Topics27_Matrix_Scaled.matrix")

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
# Extract TF motif scores and impute

matrix <- getMatrixFromProject(proj, useMatrix = "MotifMatrix", logFile = NULL, binarize = F)
matrix <- assay(matrix, "z")
rownames(matrix) <- gsub("_.*","",rownames(matrix))
matrix <- t(matrix)

matrix_imputed <- imputeMatrix(
    mat = t(matrix),
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = FALSE,
    logFile = createLogFile("imputeMatrix")
)

matrix_imputed <- t(matrix_imputed)



############################################################
# ggplot theme

library(ggplot2)

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

colorPalette <- c(
    "LRRC15+ myCAF" = "#E23C3C",
    "NGFR+ myCAF" = "#b176ef",
    "iCAF" = "#FFB132",
    "PSC" = "#89DDFF",
    "Endothelial" = "#3ECE39",
    "Unknown" = "gray88"
)



############################################################
# Plot data

# Extract GLI1 TF motif scores
plot <- matrix_imputed[,"GLI1"]
plot <- as.data.frame(plot)
colnames(plot) <- "Score"
plot$Cell_Type <- proj$celltype[match(rownames(plot),proj$cellNames)]

# Filter for cell types
plot <- plot[grep("myCAF|iCAF|PSC", plot$Cell_Type),]





p <- ggplot(plot, aes(x = Cell_Type, y = Score)) +
    geom_violin(aes(color = Cell_Type, fill = Cell_Type),
                alpha = 0.5,
                linewidth = line_width/2) +
    geom_boxplot(aes(color = Cell_Type),
                 outlier.alpha = 0,
                 lwd = line_width*0.5,
                 fatten = line_width*2,
                 width = 0.5) +
    scale_fill_manual(values = colorPalette) +
    scale_color_manual(values = colorPalette) +
    scale_x_discrete(labels = scales::label_wrap(10)) +
    default_theme +
    theme(
        legend.position = "none"
    ) +
    ylab("Deviation Z-Score (Imputed)")
p