library(ArchR)
library(dplyr)
library(parallel)
addArchRGenome("hg38")
addArchRThreads(10)
set.seed(1)

PlotPrefix = "cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR"



############################################################
# Load ArchRProject

proj <- loadArchRProject("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/")



############################################################
# Collect stromal gene signatures

signatures <- as.data.frame(readxl::read_xlsx("~/gene_sets/PDAC/PDAC_Gene_Signatures.xlsx", sheet = 3))

features <- list()
for(col in colnames(signatures)){
    features[[col]] = signatures[[col]]
    features[[col]] = features[[col]][is.na(features[[col]])==F]
    rm(col)
}

for(module in names(features)){
    gene_list = features[[module]]
    if(length(gene_list[gene_list %in% getFeatures(proj)==F])>0){
        print(paste(length(gene_list[gene_list %in% getFeatures(proj)==F]), " features not found in ", module, sep = ""))
    }
    gene_list = gene_list[gene_list %in% getFeatures(proj)]
    features[[module]] <- gene_list
    rm(gene_list,module)
}



############################################################
# Calculate per-cell module scores

set.seed(1)
proj_temp <- addModuleScore(proj, features = features, useMatrix = "GeneScoreMatrix")


# Export module scores
moduleScores <- proj_temp@cellColData
moduleScores <- moduleScores[,grep("^Module\\.",colnames(moduleScores))]
moduleScores$cellNames <- rownames(moduleScores)
moduleScores <- moduleScores[,c("cellNames",grep("^Module\\.",colnames(moduleScores),value=T))]
write.table(moduleScores, file = paste(getOutputDirectory(proj),"/Plots/",PlotPrefix,"_CellIdentities_ModuleScores",".txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)






########################################################################################################################
# Calculate IMPUTED module scores

# In this section, ArchR's implementation of MAGIC will be applied to the module scores
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
# Module scores: imputed version

# Load module scores
modules <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/Plots/cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_CellIdentities_ModuleScores.txt", header = T, sep = '\t', quote = "", comment.char = "")

# Turn modules df into a matrix
rownames(modules) <- modules$cellNames
modules <- modules[,grep("^Module\\.",colnames(modules),value = T),drop=F]
modules <- as.matrix(modules)

# Impute matrix
modules_imputed <- imputeMatrix(
    mat = t(modules),
    imputeWeights = getImputeWeights(proj),
    threads = getArchRThreads(),
    verbose = FALSE,
    logFile = createLogFile("imputeMatrix")
)

modules_imputed <- as.data.frame(t(modules_imputed))
modules_imputed$cellNames <- rownames(modules_imputed)
modules_imputed <- modules_imputed[,c(length(colnames(modules_imputed)),(1:length(colnames(modules_imputed))-1)),]

write.table(modules_imputed, file = paste(getOutputDirectory(proj),"/Plots/",PlotPrefix,"_CellIdentities_ModuleScores_Imputed",".txt",sep=""), row.names = F, col.names = T, sep = '\t', quote = F)






########################################################################################################################
# Plot module score UMAPs
########################################################################################################################


############################################################
# Load imputed module scores

modules <- read.table("~/ArchR/Primary_TSSMixed_nFrags1k_FRIS0.4/Fibroblast_MACS2ArchR_cisTopic27/Plots/cisTopic27_FibroblastTSSMixnFrags1kFRIS0.4_MACS2ArchR_CellIdentities_ModuleScores_Imputed.txt", header = T, sep = '\t', quote = "", comment.char = "")

# Add to ArchRProject
for(module_name in grep("^Module\\.",colnames(modules),value = T)){
    proj <- addCellColData(proj, cells = modules$cellNames, data = modules[[module_name]], name = module_name, force = T)
    rm(module_name)
}



############################################################
# Set directories and plot PDFs

PlotDestination = paste(getOutputDirectory(proj), "Plots", sep = "/")



for(feature_name in grep("^Module\\.",colnames(modules),value = T)){
    
    PlotSuffix = paste("_", gsub("^Module\\.Imputed\\.","",feature_name), "_Imputed.pdf", sep = "")
    PlotName = paste(PlotPrefix, PlotSuffix, sep = "")
    PlotPath = paste(PlotDestination, PlotName, sep = "/")
    PlotPath <- gsub("//","/",PlotPath)
    plot_list <- list()
    
    # Initialize plot DF, but don't show any data points (due to point size issue)
    p <- plotEmbedding(proj,
                       embedding = "cisTopic_27",
                       name = feature_name,
                       imputeWeights = NULL,
                       baseSize = 0,
                       alpha = 0,
                       labelMeans = F,
                       rastr = F,
                       plotAs = "points")
    
    
    
    # Set breakpoints for plot legend
    # break_high = (mean(p[["data"]][["color"]]) + max(p[["data"]][["color"]]))/2
    # break_low = (mean(p[["data"]][["color"]]) + min(p[["data"]][["color"]]))/2
    
    break_mid = (max(p[["data"]][["color"]]) + min(p[["data"]][["color"]]))/2
    break_high = (max(p[["data"]][["color"]]) + break_mid)/2
    break_low = (min(p[["data"]][["color"]]) + break_mid)/2
    
    break_high <- break_high*1.1
    if(0>break_low){
        break_low <- abs(break_low)*1.1
        break_low <- break_low*(-1)
    } else{
        break_low <- break_low*0.9
    }
    
    break_high <- signif(break_high, digits = 2)
    break_low <- signif(break_low, digits = 2)
    
    break_high_display <- as.character(break_high)
    break_high_display <- gsub("^0\\.|^-0\\.",".",break_high_display)
    if(0>break_high){
        break_high_display <- paste("-",break_high_display)
    }
    
    break_low_display <- as.character(break_low)
    break_low_display <- gsub("^0\\.|^-0\\.",".",break_low_display)
    if(0>break_low){
        break_low_display <- paste("-",break_low_display,sep="")
    }
    
    
    
    
    # Add data points back manually (to circumvent point size issue)
    p <- p +
        geom_point(size = 0.8, aes(stroke = 0)) +
        
        theme(
            #legend.position = "none",
            legend.text = element_text(vjust = 6.5,
                                       color = "white",
                                       #face = "bold",
                                       size = 6),
            
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)
        ) +
        scale_color_gradientn(colors = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000000),
                              breaks = c(break_low,
                                         #break_mid,
                                         break_high), labels = c(break_low_display,break_high_display)) +
        
        guides(colour=guide_colourbar(barwidth=4*0.65, 
                                      barheight = 1*0.65,
                                      ticks.linewidth = 0.5))
    
    # Add to plot list
    plot_list[[length(plot_list)+1]] <- p
    
    # Add plot w/ legend to plot list
    plot_list[[length(plot_list)+1]] <- p +
        geom_point(size = 0.8, aes(stroke = 0)) +
        
        theme(
            legend.position = "none",
            
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            
            plot.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            
            plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)
        ) +
        scale_color_gradientn(colors = colorRampPalette(ArchR::ArchRPalettes$solarExtra)(1000000))
    
    # Save plots
    pdf(file = PlotPath,
        width = 0.9, height = 0.9)
    print(plot_list)
    dev.off()
}