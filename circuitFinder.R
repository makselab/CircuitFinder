library(tidyr)
library(dplyr)
library(stringr)
library(foreach)
library(igraph)
library(gridExtra)
library(doParallel)
library(RColorBrewer)

getCircuits <- function(adjMatrix, graph){
  # run code to find circuits
  adjMatrix <- graph_from_adjacency_matrix(adjMatrix, mode = "directed", weighted = T, diag = T)
  circuits <- subgraph_isomorphisms(adjMatrix, graph, method = "vf2")
  if(length(circuits) == 0)
    return(data.frame(matrix(vector(), 0, 1, dimnames=list(c(), c("Nodes"))), stringsAsFactors=F))
  
  circuits <- foreach(i = 1:length(circuits), .combine = rbind) %do% {
    paste(
      (
        names(
          unlist(circuits[i])
        )
      ), collapse = "; "
    )
  }
  
  # reshape circuits for output
  circuits <- as.data.frame(circuits[, 1], stringsAsFactors = F)
  colnames(circuits)[1] <- "Nodes"
  return(circuits)
}

arrangeLine <- function(lineToSort) {
  lineToSort <- unlist(strsplit(lineToSort, split = "; "))
  return(paste(sort(lineToSort), collapse = "; "))
}

twoNodeSimilarity <- function(edgeList, node, nodePrime) {
  # common k-in divided by average k-in of two nodes
  # returns 1 if both have 0 k-in
  nodeInput      <- edgeList$V1[edgeList$V2 == node      & edgeList$V1 != nodePrime]
  nodeInputPrime <- edgeList$V1[edgeList$V2 == nodePrime & edgeList$V1 != node]
  if(length(nodeInput) == 0 & length(nodeInputPrime) == 0) {return(0)}
  return(2 * length(intersect(nodeInput, nodeInputPrime)) / length(c(nodeInput, nodeInputPrime)))
}

getSimilarity <- function(data, nodeIds, edgeList) {
  if(nrow(data) == 0)
    return(data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c("Nodes", "SimilarityY", "SimilarityX"))), stringsAsFactors = F))
  
  numberOfNewColumns <- length(nodeIds) / 2
  similarity <- foreach(i = 1:nrow(data), .combine = rbind) %do% {
    circuitNodes <- unlist(strsplit(data$Nodes[i], split = "; "))
    foreach(j = 1:numberOfNewColumns, .combine = c) %do% {
      round(twoNodeSimilarity(edgeList, circuitNodes[nodeIds[2 * j - 1]], circuitNodes[nodeIds[2 * j]]), digits = 2)
    }
  }
  if(nrow(data) == 1) {
    data <- as.data.frame(c(data, similarity), stringsAsFactors = F)
  } else {
    data <- cbind(data, similarity)
  }

  # this is my best shot at trying to give names to all columns depending on how many columns we have
  # and the best shot at sorting dataframe by similarity starting from first finishing last
  colnames(data)[1 + 1:numberOfNewColumns] <- paste("Similarity", 1:numberOfNewColumns, sep = "")
  for(i in (1 + numberOfNewColumns:1)) {
    data <- data[sort(data[, c(i)], decreasing = T, index.return = T)$ix, ]
  }
  return(data)
}

makePngs <- function(dataset, circuitType, circuits, fullEdgeList) {
  if(length(circuits) == 0) {return()}
  system(paste("mkdir ", "circuits/", dataset, "_", circuitType, sep = ""))
  if(ncol(fullEdgeList) == 3) {
    fullEdgeList$color <- group_indices(fullEdgeList, V3)
    numberOfColors <- max(fullEdgeList$color)
    if(numberOfColors < 9 & numberOfColors > 2) {
      edgeColors <- brewer.pal(numberOfColors, "Set1")
    } else {
      edgeColors <- rainbow(numberOfColors)
    }
    fullEdgeList$color <- edgeColors[fullEdgeList$color]
    
    legendColors <- foreach(f = 1:numberOfColors, .combine = c) %do% {first(fullEdgeList$V3[grepl(edgeColors[f], fullEdgeList$color)])}
  }
  for(i in 1:length(circuits)) {
    circuitNodes <- unlist(strsplit(circuits[i], split = "; "))
    circuitEdges <- fullEdgeList[fullEdgeList$V1 %in% circuitNodes & fullEdgeList$V2 %in% circuitNodes, ]
    circuitEdges[, 1:2] <- apply(circuitEdges[, 1:2], 2, toupper)
    
    network <- graph_from_data_frame(d = circuitEdges, vertices = toupper(circuitNodes), directed = T)
    V(network)$label.size <- 30
    
    png(filename = paste("circuits/", dataset, "_", circuitType, "/", i, ".png", sep = ""),
        width = 640, height = 360)
    oldMargins<-par("mar")
    par(mar = c(0, 0, 0, 0))
    if(ncol(fullEdgeList) == 4) {
      plot(network, edge.color = circuitEdges$color, vertex.label.cex = 2.5)
      legend(x = 1.2, y = 1.1, legend = legendColors,
             col = edgeColors, lty = 1, lwd = 3, cex = 1,
             title = circuitType, text.font = 4, bg = 'white')
    } else {
      plot(network, vertex.label.cex = 2.5)
    }
    par(mar = oldMargins)
    dev.off()
  }
}

checkCircuitRemovalCondition <- function(circuit_node_names, circuit, graph) {
  # here we check if adjacency matrix of subgraph of reduced (without self-loops and multilinks) is equal to adjacency matrix searched for
  circuit_nodes <- unlist(strsplit(circuit_node_names, split = "; "))
  circuit_Subgraph <- induced_subgraph(graph, circuit_nodes, impl = "create_from_scratch")
  circuit_Subgraph_nodes <- vertex_attr(circuit_Subgraph, "name")

  permutation_order <- unlist(lapply(circuit_Subgraph_nodes, function(circuit_Subgraph_node) grep(paste("^", circuit_Subgraph_node, "$", sep = ""), circuit_nodes)))
  circuit_Subgraph <- permute(circuit_Subgraph, permutation_order)
  
  found_adjacency <- as.matrix(as_adjacency_matrix(circuit_Subgraph, names = T))
  return(all(found_adjacency == circuit))
}

# source("~/Dropbox/groupoid_finding_codes/fibers/R/classifier.R")
# dataset <- gsub(".*/", "", files[idx])
# dataset <- gsub("\\.txt", "", dataset)
# print(paste("Running", dataset))
# prefix <- paste("/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData/output", dataset, sep = "/")
# classifiedStructures <- getBlocks(prefix)
# regulators <- classifiedStructures$Regulators[classifiedStructures$Class == "Unsynchronized Star Fiber"]
# i <- 1
# for(i in 1:nrow(circuits)) {
#   circuits_Nodes <- unlist(strsplit(circuits$Nodes[i], split = "; "))
#   circuits$Regulating[i] <- sum(circuits_Nodes[c(4, 5)] %in% regulators)
#   circuits$ForwardLink[i] <- paste(fullEdgeList$V3[fullEdgeList$V1 %in% circuits_Nodes[4] & fullEdgeList$V2 %in% circuits_Nodes[5]], collapse = "; ")
#   circuits$BackwardLink[i] <- paste(fullEdgeList$V3[fullEdgeList$V1 %in% circuits_Nodes[5] & fullEdgeList$V2 %in% circuits_Nodes[4]], collapse = "; ")
#   circuits$upLink[i] <- paste(fullEdgeList$V3[fullEdgeList$V1 %in% circuits_Nodes[2] & fullEdgeList$V2 %in% circuits_Nodes[4]], collapse = "; ")
#   circuits$downLink[i] <- paste(fullEdgeList$V3[fullEdgeList$V1 %in% circuits_Nodes[3] & fullEdgeList$V2 %in% circuits_Nodes[5]], collapse = "; ")
#   circuits$Pair_Nodes[i] <- paste(sort(circuits_Nodes[c(4, 5)]), collapse = "; ")
# }
# circuits$BackAndForthRepression <- grepl("Repression", circuits$ForwardLink) + grepl("Repression", circuits$BackwardLink)
# circuits$updownRepression <- grepl("Repression", circuits$upLink) + grepl("Repression", circuits$downLink)
# circuits <- circuits %>%
#   arrange(desc(BackAndForthRepression), desc(updownRepression), Regulating)
# grep("BRCA1; ESR1", circuits$Nodes)
# 
# circuit_FFF
# circuit_JK
# circuit <- circuit_AR
# similarity_nodes <- c(4, 5, 2, 3)
# name <- "FFF"
findCircuits <- function(circuit, similarity_nodes, name, dataset, fullEdgeList, pvalues) {
  # returns number of circuits found and writes circuits to files
  nodes <- unique(c(fullEdgeList$V1, fullEdgeList$V2))
  edgeList <- fullEdgeList %>%
    filter(V1 != V2)
  edgeList <- edgeList[!duplicated(edgeList[, c(1, 2)]), ]
  graph <- graph.data.frame(edgeList)
  
  circuits <- getCircuits(circuit, graph)
  
  # remove circuits, which are rearrangement of one another
  duplicationTable <- sapply(circuits$Nodes, function(x) arrangeLine(x))
  circuits <- circuits[!duplicated(duplicationTable), ]
  
  if(length(circuits) != 0) {
    exactCircuits <- foreach(i = 1:length(circuits), .combine = c) %do% {checkCircuitRemovalCondition(circuits[i], circuit, graph)}
    circuits <- circuits[exactCircuits]
  }
  
  circuits <- as.data.frame(circuits, stringsAsFactors = F)
  colnames(circuits) <- "Nodes"
  
  # calculate measure of similarity in between Xs and Ys
  circuits <- getSimilarity(circuits, similarity_nodes, edgeList)
  
  if(length(circuits) != 0 & !pvalues) {
    write(unlist(unite(circuits, Output, sep = "\t")),
          paste("circuits/", dataset, "_", name, ".txt", sep = ""))
  }
  
  if(!pvalues) {
    makePngs(dataset, name, circuits$Nodes, fullEdgeList)
  }
  return(nrow(circuits))
}

###########################################################
########## End of functions, beginning of script ##########
###########################################################

networkDirectory <- "/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData"
#"/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/cleanedData"
files <- list.files(networkDirectory, full.names = T)

fileNames <- "modelSample"
  "Regulations_in_ATRM
regulations_bacilus_cleaned
Ecoli
human_network_tf_gene_weighted
trrust_rawdata.human
Micobacterium_tuberculosis
trrust_rawdata.mouse
salmonella_Typhi_SL1344_RN
SGD_full_yeast_TF-G
YTRP_TF_gene_regulatory"

fileNames <- unlist(strsplit(fileNames, split = "\n"))
idxs <- foreach(i = 1:length(fileNames), .combine = c) %do% {grep(paste("/", fileNames[i], sep = ""), files)}
idxs <- idxs[!grepl("KEGG", files[idxs])]

workingDirectory <- "~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN"
setwd(workingDirectory)
pvalues = T

# loop over all files
myCluster <- makeCluster(detectCores() - 3, outfile = "")
registerDoParallel(myCluster)
pb <- txtProgressBar(min(idxs), max(idxs), style = 3)
start.time <- Sys.time()

summary <- foreach(idx = idxs, .combine = rbind, .errorhandling = "remove") %do% {
  # setTxtProgressBar(pb, idx)
  
  library(tidyr, quietly = T)
  library(dplyr, quietly = T)
  library(stringr, quietly = T)
  library(foreach, quietly = T)
  library(igraph, quietly = T)
  
  #print(files[idx])
  dataset <- gsub(".*/", "", files[idx])
  dataset <- gsub("\\.txt", "", dataset)
  
  fullEdgeList <- read.table(files[idx], stringsAsFactors = F)
  
  circuit_AR <- matrix(
    data = c(0, 1,
             1, 0), ncol = 2)
  circuit_count_AR <- findCircuits(circuit_AR, c(1, 2), "AR", dataset, fullEdgeList, pvalues)
  
  circuit_FFF <- matrix(
    data = c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             0, 1, 0, 0, 1,
             0, 0, 1, 1, 0), ncol = 5)
  circuit_count_FFF <- findCircuits(circuit_FFF, c(4, 5, 2, 3), "FFF", dataset, fullEdgeList, pvalues)

  # circuit_JK <- matrix(
  #   data = c(0, 0, 0, 0, 0,
  #            1, 0, 0, 0, 1,
  #            1, 0, 0, 1, 0,
  #            0, 1, 0, 0, 1,
  #            0, 0, 1, 1, 0), ncol = 5)
  # circuit_count_JK <- findCircuits(circuit_JK, c(4, 5, 2, 3), "Fibonacci", dataset, fullEdgeList, pvalues)
  
  circuit_JK_no_clock <- matrix(
    data = c(0, 0, 0, 1,
             0, 0, 1, 0,
             1, 0, 0, 1,
             0, 1, 1, 0), ncol = 4)
  circuit_count_JK <- findCircuits(circuit_JK_no_clock, c(3, 4, 1, 2), "JK_no_clock", dataset, fullEdgeList, pvalues)
  
  # write data to files as a backup in case code bugs
  if(pvalues == T) {
    write(paste(dataset, circuit_count_AR, circuit_count_FFF, circuit_count_JK, sep = ";"),
          "~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/pvalue_runs_summary.txt",
          append = T)
  }
  
  c(dataset, circuit_count_AR, circuit_count_FFF, circuit_count_JK)
}
now.time <- Sys.time()
stopCluster(myCluster)
time.taken <- now.time - start.time
print(time.taken)

summary <- as.data.frame(summary, stringsAsFactors = F)
colnames(summary) <- c("Dataset", "AR", "FFF", "JK")
summary[, -1] <- sapply(summary[, -1], as.numeric)

#write.csv(summary, "pvalueSummary.csv", row.names = F, quote = F)

# pdf("output.pdf", height = 11, width = 14)
# grid.table(summary)
# dev.off()
# 
# # here we count p values
# zscoretopvalue <- function(x) {
#   return(1 - (erf(x) + 1) / 2)
# }
# 
# summary$Dataset <- gsub("(modelSample|[0-9])", "", summary$Dataset)
# 
# meanSummary <- summary %>%
#   group_by(Dataset) %>%
#   summarise_if(is.numeric, mean)
# 
# sdSummary <- summary %>%
#   group_by(Dataset) %>%
#   summarise_if(is.numeric, sd)
# 
# pvalueSummary <- cbind(meanSummary, sdSummary)
# pvalueSummary <- pvalueSummary[, -5]
# colnames(pvalueSummary)[2:7] <- apply(expand.grid(colnames(pvalueSummary)[2:4], c("Mean", "SD")), 1, paste, collapse="")

getDatasetByPrefix <- function(prefix) {
  prefix <- gsub("Regulations_in_ATRM",            "Arabidopsis Thaliana", prefix)
  prefix <- gsub("regulations_bacilus_cleaned",    "Bacillus subtilis", prefix)
  prefix <- gsub("Micobacterium_tuberculosis",     "Micobacterium tuberculosis", prefix)
  prefix <- gsub("salmonella_Typhi_SL1344_RN",     "Salmonella SL1344", prefix)
  prefix <- gsub("SGD_full_yeast_TF-G",            "Yeast SGD", prefix)
  prefix <- gsub("YTRP_TF_gene_regulatory",        "Yeast YTRP", prefix)
  prefix <- gsub("trrust_rawdata.mouse",           "Mouse TRRUST", prefix)
  prefix <- gsub("trrust_rawdata.human",           "Human TRRUST2", prefix)
  prefix <- gsub("human_network_tf_gene_weighted", "Human TRRUST", prefix)
  prefix <- gsub("combinedKEGG",                   "Human KEGG", prefix)
  return(prefix)
}

# here we do pdfs

structuresPath <- "~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/"
setwd(structuresPath)
circuitDirectories <- list.dirs(paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/"), full.names = T)
circuitDirectories <- circuitDirectories[-1]
circuitDirectories <- as.data.frame(circuitDirectories, stringsAsFactors = F)

circuitDirectories$Dataset <- gsub(".*//(.*)_[[:alpha:]]*$", "\\1", circuitDirectories$circuitDirectories)
#circuitDirectories <- circuitDirectories[!grepl("human_network_tf_gene_weighted", circuitDirectories$circuitDirectories), ]
#circuitDirectories <- circuitDirectories[foreach(i = 1:length(fileNames), .combine = c) %do% {grep(paste(fileNames[i], sep = ""), circuitDirectories$Dataset)}, ]

samplePNGs = T
# there is roughly 4 times more asymm circuits than symm, so we sample to get equal amount of both types in structures
asymmProportionToSample = 0.15
symmProportionToSample = 0.6

texFileName <- "structures.tex"

write("\\documentclass[preprint,aps,preprintnumbers,amsmath,amssymb]{revtex4}

\\usepackage[]{graphicx}

\\begin{document}", texFileName)

# first we get all asymmetric circuits

foreach(dataset = unique(circuitDirectories$Dataset), .combine = c) %do% {
  directoriesToCheck <- circuitDirectories %>%
    filter(Dataset == dataset)

  AR_circuits <- try(read.table(paste("circuits/", dataset, "_AR.txt", sep = ""), sep = "\t", stringsAsFactors = F))
  FFF_circuits <- try(read.table(paste("circuits/", dataset, "_FFF.txt", sep = ""), sep = "\t", stringsAsFactors = F))
  Fibonacci_circuits <- try(read.table(paste("circuits/", dataset, "_Fibonacci.txt", sep = ""), sep = "\t", stringsAsFactors = F))

  dataset <- getDatasetByPrefix(dataset)

  pngFileNames <- foreach(i = 1:nrow(directoriesToCheck), .combine = c) %do% {
    paste(directoriesToCheck$circuitDirectories[i], list.files(directoriesToCheck$circuitDirectories[i]), sep = "/")
  }

  if(samplePNGs) {
    fibonacciPNG <- pngFileNames[grepl("Fibonacci", pngFileNames)]
    otherPNG <- pngFileNames[!grepl("Fibonacci", pngFileNames)]
    otherPNG <- otherPNG[sample(1:length(otherPNG), as.integer(length(pngFileNames) * asymmProportionToSample))]
    pngFileNames <- c(otherPNG, fibonacciPNG)
  }

  for(i in 1:length(pngFileNames)) {
    type <- gsub(".*_([[:alpha:]]*)/[0-9]*.png$", "\\1", pngFileNames[i])
    idx <- as.integer(gsub(".*/([0-9]*).png$", "\\1", pngFileNames[i]))
    if(type == "AR")
      nodes <- toupper(AR_circuits$V1[idx])
    if(type == "FFF")
      nodes <- toupper(FFF_circuits$V1[idx])
    if(type == "Fibonacci")
      nodes <- toupper(Fibonacci_circuits$V1[idx])
    nodes <- gsub("_", "\\\\_", nodes)
    
    write(paste("\\clearpage ", dataset, " \\\\ ",
                paste(type, "Broken"),
                " \\\\ ", nodes, " \\\\", sep = ""),
          file = texFileName, append = T)
    write(paste("\\includegraphics[width=\\textwidth, height=100in, keepaspectratio]{",
                gsub(".*/(circuits)", "\\1", pngFileNames[i]),
                "} \\\\", sep = ""),
          file = texFileName, append = T)
  }
}

# now we get all symmetric circuits
writeToStructures <- function(prefix, type, subblocks, structuresPath, symmStructuresPath, texFileName, samplePNGs, proportionToSample) {
  if(nrow(subblocks) == 0) {return()}
  subblocks <- subblocks[str_count(subblocks$Regulators, ",") < 50, ]
  
  idxs = 1:nrow(subblocks)
  if(samplePNGs) {
    idxs <- sample(1:nrow(subblocks), as.integer(nrow(subblocks) * proportionToSample))
  }
  
  for(i in idxs) {
    write(paste("\\clearpage ", getDatasetByPrefix(prefix), " \\\\ ",
                type, " \\\\ ",
                "Regulators: ", gsub("_", "\\\\_", toupper(subblocks$Regulators[i])), " \\\\ ",
                "Fiber: ", gsub("_", "\\\\_", toupper(subblocks$Fiber[i])), " \\\\",
                sep = ""),
          file = paste(structuresPath, texFileName, sep = ""), append = T)
    write(paste("\\includegraphics[width=\\textwidth, height=100in, keepaspectratio]{",
                symmStructuresPath, "/output/", prefix, "buildingBlocks/", subblocks$Id[i], ".png",
                "} \\\\", sep = ""),
          file = paste(structuresPath, texFileName, sep = ""), append = T)
  }
}

symmStructuresPath <- "/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData"
source("~/Dropbox/groupoid_finding_codes/fibers/R/classifier.R")

prefixes <- c(unique(circuitDirectories$Dataset), "combinedKEGG")

for(prefix in prefixes) {
  if(prefix == "combinedKEGG") {
    blocks <- read.csv("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/blocks.txt", sep = "\t", header = T)
  } else {
    blocks <- getBlocks(paste(symmStructuresPath, "/output/", prefix, sep = ""))
  }

  ARBlocks <- blocks[grepl("(Chain|Chain-Star|Synchronized Star Fiber|Repression Chain)", blocks$Class), ]
  FFFBlocks <- blocks[grepl("(Feed-Forward Fiber|UNSAT Feed-Forward Fiber)", blocks$Class), ]
  FibonacciBlocks <- blocks[grepl("(Feedback Fiber|Fibonacci n = 1)", blocks$Class), ]
  n2Blocks <- blocks[grepl("(n > 1|Negative n = 2|Unclassified)", blocks$Class), ]

  writeToStructures(prefix, "AR", ARBlocks, structuresPath, symmStructuresPath, texFileName, samplePNGs, symmProportionToSample)
  writeToStructures(prefix, "FFF", FFFBlocks, structuresPath, symmStructuresPath, texFileName, samplePNGs, symmProportionToSample)
  writeToStructures(prefix, "Fibonacci", FibonacciBlocks, structuresPath, symmStructuresPath, texFileName, samplePNGs, symmProportionToSample)
  writeToStructures(prefix, "n = 2", n2Blocks, structuresPath, symmStructuresPath, texFileName, F, 1)
}

write("\\end{document}", texFileName, append = T)

system(paste("pdflatex ", texFileName))

# system(paste("mv ",
#              gsub(".tex", ".pdf", texFileName),
#              " SM.pdf",
#              gsub("[^/]*.tex", paste(dataset, ".pdf", sep = ""), texFileName),
#              sep = ""))



# this is old code that finds all types of circuits

# summary <- foreach(idx = idxs, .combine = rbind, .errorhandling = "remove") %do% {
#   library(tidyr)
#   library(dplyr)
#   library(stringr)
#   library(foreach)
#   library(igraph)
#   library(gridExtra)
#   print(idx)
#   fullEdgeList <- read.table(files[idx], stringsAsFactors = F)
#   nodes <- unique(c(fullEdgeList$V1, fullEdgeList$V2))
#   edgeList <- fullEdgeList %>%
#     filter(V1 != V2)
#   edgeList <- edgeList[!duplicated(edgeList[, c(1, 2)]), ]
#   graph <- graph.data.frame(edgeList)
#   
#   # AR section
#   backNForth <- matrix(
#     data = c(0, 1,
#              1, 0), ncol = 2)
#   
#   backNForth <- getCircuits(backNForth, graph, "backNForth")
#   # remove duplicated lines
#   backNForth$Nodes <- sapply(backNForth$Nodes, function(x) arrangeLine(x))
#   backNForth <- backNForth[!duplicated(backNForth$Nodes), ]
#
# source("~/Dropbox/groupoid_finding_codes/fibers/R/classifier.R")
#   dataset <- gsub(".*/", "", files[idx])
#   dataset <- gsub("\\.txt", "", dataset)
#   print(paste("Running", dataset))
#   prefix <- paste("/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData/output", dataset, sep = "/")
#   classifiedStructures <- getBlocks(prefix)
#   
#   possibleARb <- nrow(backNForth)
#   possibleARbList <- backNForth$Nodes
#   
#   possibleARa <- 0
#   possibleARaList <- NULL
#   for(i in 1:nrow(backNForth)) {
#     possibleAR <- unlist(strsplit(backNForth$Nodes[i], split = "; "))
#     if(nrow(edgeList[edgeList$V2 == possibleAR[1], ]) == 2) {
#       if(nrow(edgeList[edgeList$V2 == possibleAR[2], ]) == 2) {
#         possibleARa <- possibleARa + 1
#         possibleARaList <- c(possibleARaList, T)
#         next
#       }
#     }
#     possibleARaList <- c(possibleARaList, F)
#   }
#   possibleARaList <- backNForth$Nodes[possibleARaList]
#   
#   possibleARc <- sum(classifiedStructures$nl == "Fibonacci")
#   possibleARcList <- gsub(",", ";", classifiedStructures$Node[classifiedStructures$nl == "Fibonacci"])
#   
#   # write to output files
#   # write(possibleARaList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_ARa.txt", sep = ""))
#   # write(possibleARbList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_ARb.txt", sep = ""))
#   # write(possibleARcList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_ARc.txt", sep = ""))
#   
#   # FFF without fiber section
#   FFFnoFiberA <- 0
#   FFFnoFiberAlist <- NULL
#   FFFnoFiberB <- 0
#   FFFnoFiberBlist <- NULL
#   for(i in 1:length(possibleARaList)) {
#     if(length(possibleARaList) == 0) {break}
#     possibleY <- unlist(strsplit(possibleARaList[i], split = "; "))
#     possibleYprime <- possibleY[1]
#     possibleY <- possibleY[2]
#     possibleX      <- edgeList$V1[edgeList$V2 == possibleY      & edgeList$V1 != possibleYprime]
#     possibleXprime <- edgeList$V1[edgeList$V2 == possibleYprime & edgeList$V1 != possibleY]
#     if(possibleX != possibleXprime) {
#       FFFnoFiberA <- FFFnoFiberA + 1
#       FFFnoFiberAlist <- c(FFFnoFiberAlist, paste(possibleX, possibleXprime, possibleY, possibleYprime, sep = "; "))
#       
#       possibleXclock      <- edgeList$V1[edgeList$V2 == possibleX]
#       possibleXprimeClock <- edgeList$V1[edgeList$V2 == possibleXprime]
#       
#       possibleClock <- intersect(possibleXclock, possibleXprimeClock)
#       if(length(possibleClock) != 0) {
#         FFFnoFiberB <- FFFnoFiberB + length(possibleClock)
#         newEntries <- foreach(clock = possibleClock) %do% {
#           paste(clock, possibleX, possibleXprime, possibleY, possibleYprime, sep = "; ")
#         }
#         FFFnoFiberBlist <- c(FFFnoFiberBlist, unlist(newEntries))
#       }
#     }
#   }
#   
#   SRnoClock <- matrix(
#     data = c(0, 0, 0, 0,
#              0, 0, 0, 0,
#              1, 0, 0, 1,
#              0, 1, 1, 0), ncol = 4)
#   SRnoClock <- getCircuits(SRnoClock, graph, "NoClock")
#   # remove circuits, which are rearrangement of one another
#   duplicationTable <- sapply(SRnoClock$Nodes, function(x) arrangeLine(x))
#   SRnoClock <- SRnoClock[!duplicated(duplicationTable), ]
#   
#   # calculate measure of similarity in between Xs and Ys
#   SRnoClock <- getSimilarity(SRnoClock, c(3, 4, 1, 2), edgeList)
#   
#   SRclock <- matrix(
#     data = c(0, 0, 0, 0, 0,
#              1, 0, 0, 0, 0,
#              1, 0, 0, 0, 0,
#              0, 1, 0, 0, 1,
#              0, 0, 1, 1, 0), ncol = 5)
#   SRclock <- getCircuits(SRclock, graph, "Clock")
#   # remove circuits, which are rearrangement of one another
#   duplicationTable <- sapply(SRclock$Nodes, function(x) arrangeLine(x))
#   SRclock <- SRclock[!duplicated(duplicationTable), ]
#   
#   # calculate measure of similarity in between Xs and Ys
#   SRclock <- getSimilarity(SRclock, c(4, 5, 2, 3), edgeList)
#   
#   # oldList <- JKclock
#   # newList <- oldList[0, ]
#   # nodesInBlocks <- unique(unlist(strsplit(oldList$Nodes, split = "; ")))
#   # nodesInBlocks
#   # pointer <- 1
#   # while(pointer <= length(nodesInBlocks)) {
#   #   # take node on pointer in a list of available nodes
#   #   blockIds <- grep(paste("(^| )", nodesInBlocks[pointer], "($|;)", sep = ""), SRclock$Nodes)
#   #   blockIds <- sample(blockIds, replace = F)
#   #   
#   #   # for loop in blockIds until find block, all nodes of which are available
#   #   increasePointer <- T
#   #   for(blockId in blockIds) {
#   #     block <- unlist(strsplit(oldList$Nodes[blockId], split = "; "))
#   #     
#   #     # if all nodes from this block are still in the list
#   #     #   add block with this id to new list
#   #     #   remove all nodes from this block from list of available nodes
#   #     #   break cycle
#   #     if(all(block %in% nodesInBlocks)) {
#   #       increasePointer <- F
#   #       
#   #       newList <- rbind(newList, oldList[blockId, ])
#   #       oldList <- oldList[c(-blockIds), ]
#   #       
#   #       nodesInBlocks <- nodesInBlocks[!nodesInBlocks %in% block]
#   #       pointer <- 1
#   #       break
#   #     }
#   #   }
#   #   if(increasePointer)
#   #     pointer <- pointer + 1
#   # }
#   
#   # write to output files
#   # write(FFFnoFiberAlist, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_FFFa.txt", sep = ""))
#   # write(FFFnoFiberBlist, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_FFFb.txt", sep = ""))
#   # write(unlist(unite(SRclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
#   #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_FFFc.txt", sep = ""))
#   # write(unlist(unite(SRclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
#   #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_FFFd.txt", sep = ""))
#   
#   # JK section
#   JKnoClock <- matrix(
#     data = c(0, 0, 0, 1,
#              0, 0, 1, 0,
#              1, 0, 0, 1,
#              0, 1, 1, 0), ncol = 4)
#   JKnoClock <- getCircuits(JKnoClock, graph, "NoClock")
#   # remove circuits, which are rearrangement of one another
#   duplicationTable <- sapply(JKnoClock$Nodes, function(x) arrangeLine(x))
#   JKnoClock <- JKnoClock[!duplicated(duplicationTable), ]
#   
#   # calculate measure of similarity in between Xs and Ys
#   JKnoClock <- getSimilarity(JKnoClock, c(3, 4, 1, 2), edgeList)
#   
#   JKclock <- matrix(
#     data = c(0, 0, 0, 0, 0,
#              1, 0, 0, 0, 1,
#              1, 0, 0, 1, 0,
#              0, 1, 0, 0, 1,
#              0, 0, 1, 1, 0), ncol = 5)
#   JKclock <- getCircuits(JKclock, graph, "Clock")
#   # remove circuits, which are rearrangement of one another
#   duplicationTable <- sapply(JKclock$Nodes, function(x) arrangeLine(x))
#   JKclock <- JKclock[!duplicated(duplicationTable), ]
#   
#   # calculate measure of similarity in between Xs and Ys
#   JKclock <- getSimilarity(JKclock, c(4, 5, 2, 3), edgeList)
#   
#   # write to output files
#   # write(unlist(unite(JKclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
#   #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_JKa.txt", sep = ""))
#   # write(unlist(unite(JKnoClock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
#   #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits/", dataset, "_JKb.txt", sep = ""))
#   
#   # makePngs(dataset, "ARa", possibleARaList)
#   # makePngs(dataset, "ARb", possibleARbList)
#   # makePngs(dataset, "ARc", possibleARcList)
#   # 
#   # makePngs(dataset, "FFFa", FFFnoFiberAlist)
#   # makePngs(dataset, "FFFb", FFFnoFiberBlist)
#   # makePngs(dataset, "FFFc", SRclock$Nodes)
#   # makePngs(dataset, "FFFd", SRnoClock$Nodes)
#   # 
#   # makePngs(dataset, "JKa", JKclock$Nodes)
#   # makePngs(dataset, "JKb", JKnoClock$Nodes)
#   c(dataset, possibleARa, possibleARb, possibleARc, FFFnoFiberA, FFFnoFiberB, nrow(SRclock), nrow(SRnoClock), nrow(JKclock), nrow(JKnoClock))
# }