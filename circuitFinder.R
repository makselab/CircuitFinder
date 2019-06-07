library(tidyr)
library(dplyr)
library(stringr)
library(foreach)
library(igraph)
library(gridExtra)
library(doParallel)
library(RColorBrewer)

source("~/Dropbox/groupoid_finding_codes/fibers/R/classifier.R")

getCircuits <- function(adjMatrix, graph, type){
  # run code to find circuits
  adjMatrix <- graph_from_adjacency_matrix(adjMatrix, mode = "directed", weighted = T, diag = T)
  circuits <- subgraph_isomorphisms(adjMatrix, graph, method = "vf2")
  if(length(circuits) == 0)
    return(data.frame(matrix(vector(), 0, 2, dimnames=list(c(), c("Nodes", "CircuitType"))), stringsAsFactors=F))
  
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
  circuits$CircuitType = type
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
    return(data.frame(matrix(vector(), 0, 4, dimnames=list(c(), c("Nodes", "CircuitType", "SimilarityY", "SimilarityX"))), stringsAsFactors = F))
  
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
  colnames(data)[c(3, 4)] <- c("SimilarityY", "SimilarityX")
  data <- arrange(data, desc(SimilarityY), desc(SimilarityX))
  return(data)
}

makePngs <- function(dataset, circuitType, circuits) {
  if(length(circuits) == 0) {return()}
  system(paste("mkdir ", "/home/ian/Dropbox/Research/PhD\\ work/shared\\ folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_", circuitType, sep = ""))
  if(ncol(fullEdgeList) == 3) {
    fullEdgeList$color <- group_indices(fullEdgeList, V3)
    numberOfColors <- max(fullEdgeList$color)
    if(numberOfColors < 9) {
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
    
    network <- graph_from_data_frame(d = circuitEdges, vertices = circuitNodes, directed = T)
    V(network)$label.size <- 30
    
    png(filename = paste("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_", circuitType, "/", i, ".png", sep = ""),
        width = 1280, height = 720)
    oldMargins<-par("mar")
    par(mar = c(0, 0, 0, 0))
    if(ncol(fullEdgeList) == 3) {
      plot(network, edge.color = circuitEdges$color, vertex.label.cex = 2.5)
      legend(x = 1.5, y = 1.1, legend = legendColors,
             col = edgeColors, lty = 1, lwd = 3, cex = 1,
             title = circuitType, text.font = 4, bg = 'white')
    } else {
      plot(network, vertex.label.cex = 2.5)
    }
    par(mar = oldMargins)
    dev.off()
  }
}

###########################################################
########## End of functions, beginning of script ##########
###########################################################

networkDirectory <- "/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData"
#"/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData"
files <- list.files(networkDirectory, full.names = T)

fileNames <- 
  "trrust_rawdata.human"
  "Regulations_in_ATRM
B_subtilis_quantitative_transcription_network
regulations_bacilus
Ecoli
human_network_tf_gene
trrust_rawdata.human
Micobacterium_tuberculosis
trrust_rawdata.mouse
salmonella_Typhi_SL1344_RN
salmonella_Typhi_T2_RN
evidence_yeast_TF-G
full_yeast_TF-G
SGD_full_yeast_TF-G
Yeast_TRN
YTRP_TF_gene_binding
YTRP_TF_gene_regulatory"
fileNames <- unlist(strsplit(fileNames, split = "\n"))
idxs <- foreach(i = 1:length(fileNames), .combine = c) %do% {grep(paste("/", fileNames[i], sep = ""), files)}

# loop over all files
#idx = idxs[1]
#files
#idxs = c(3, 4, 7, 8)
#myCluster <- makeCluster(detectCores() - 2, outfile = "")
#registerDoParallel(myCluster)
summary <- foreach(idx = idxs, .combine = rbind, .errorhandling = "remove") %do% {
  library(tidyr)
  library(dplyr)
  library(stringr)
  library(foreach)
  library(igraph)
  library(gridExtra)
  print(idx)
  fullEdgeList <- read.table(files[idx], stringsAsFactors = F)
  nodes <- unique(c(fullEdgeList$V1, fullEdgeList$V2))
  edgeList <- fullEdgeList %>%
    filter(V1 != V2)
  edgeList <- edgeList[!duplicated(edgeList[, c(1, 2)]), ]
  graph <- graph.data.frame(edgeList)

  # AR section
  backNForth <- matrix(
    data = c(0, 1,
             1, 0), ncol = 2)
  
  backNForth <- getCircuits(backNForth, graph, "backNForth")
  # remove duplicated lines
  backNForth$Nodes <- sapply(backNForth$Nodes, function(x) arrangeLine(x))
  backNForth <- backNForth[!duplicated(backNForth$Nodes), ]
  
  dataset <- gsub(".*/", "", files[idx])
  dataset <- gsub("\\.txt", "", dataset)
  print(paste("Running", dataset))
  prefix <- paste("/home/ian/Dropbox/groupoid_finding_codes/naturePhysRuns/automaticRunData/output", dataset, sep = "/")
  classifiedStructures <- getBlocks(prefix)
  
  possibleARb <- nrow(backNForth)
  possibleARbList <- backNForth$Nodes
  
  possibleARa <- 0
  possibleARaList <- NULL
  for(i in 1:nrow(backNForth)) {
    possibleAR <- unlist(strsplit(backNForth$Nodes[i], split = "; "))
    if(nrow(edgeList[edgeList$V2 == possibleAR[1], ]) == 2) {
      if(nrow(edgeList[edgeList$V2 == possibleAR[2], ]) == 2) {
        possibleARa <- possibleARa + 1
        possibleARaList <- c(possibleARaList, T)
        next
      }
    }
    possibleARaList <- c(possibleARaList, F)
  }
  possibleARaList <- backNForth$Nodes[possibleARaList]
  
  possibleARc <- sum(classifiedStructures$nl == "Fibonacci")
  possibleARcList <- gsub(",", ";", classifiedStructures$Node[classifiedStructures$nl == "Fibonacci"])
  
  # write to output files
  # write(possibleARaList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_ARa.txt", sep = ""))
  # write(possibleARbList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_ARb.txt", sep = ""))
  # write(possibleARcList, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_ARc.txt", sep = ""))
  
  # FFF without fiber section
  FFFnoFiberA <- 0
  FFFnoFiberAlist <- NULL
  FFFnoFiberB <- 0
  FFFnoFiberBlist <- NULL
  for(i in 1:length(possibleARaList)) {
    if(length(possibleARaList) == 0) {break}
    possibleY <- unlist(strsplit(possibleARaList[i], split = "; "))
    possibleYprime <- possibleY[1]
    possibleY <- possibleY[2]
    possibleX      <- edgeList$V1[edgeList$V2 == possibleY      & edgeList$V1 != possibleYprime]
    possibleXprime <- edgeList$V1[edgeList$V2 == possibleYprime & edgeList$V1 != possibleY]
    if(possibleX != possibleXprime) {
      FFFnoFiberA <- FFFnoFiberA + 1
      FFFnoFiberAlist <- c(FFFnoFiberAlist, paste(possibleX, possibleXprime, possibleY, possibleYprime, sep = "; "))
    
      possibleXclock      <- edgeList$V1[edgeList$V2 == possibleX]
      possibleXprimeClock <- edgeList$V1[edgeList$V2 == possibleXprime]
    
      possibleClock <- intersect(possibleXclock, possibleXprimeClock)
      if(length(possibleClock) != 0) {
        FFFnoFiberB <- FFFnoFiberB + length(possibleClock)
        newEntries <- foreach(clock = possibleClock) %do% {
          paste(clock, possibleX, possibleXprime, possibleY, possibleYprime, sep = "; ")
        }
        FFFnoFiberBlist <- c(FFFnoFiberBlist, unlist(newEntries))
      }
    }
  }

  SRnoClock <- matrix(
    data = c(0, 0, 0, 0,
             0, 0, 0, 0,
             1, 0, 0, 1,
             0, 1, 1, 0), ncol = 4)
  SRnoClock <- getCircuits(SRnoClock, graph, "NoClock")
  # remove circuits, which are rearrangement of one another
  duplicationTable <- sapply(SRnoClock$Nodes, function(x) arrangeLine(x))
  SRnoClock <- SRnoClock[!duplicated(duplicationTable), ]
  
  # calculate measure of similarity in between Xs and Ys
  SRnoClock <- getSimilarity(SRnoClock, c(3, 4, 1, 2), edgeList)
  
  SRclock <- matrix(
    data = c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             1, 0, 0, 0, 0,
             0, 1, 0, 0, 1,
             0, 0, 1, 1, 0), ncol = 5)
  SRclock <- getCircuits(SRclock, graph, "Clock")
  # remove circuits, which are rearrangement of one another
  duplicationTable <- sapply(SRclock$Nodes, function(x) arrangeLine(x))
  SRclock <- SRclock[!duplicated(duplicationTable), ]
  
  # calculate measure of similarity in between Xs and Ys
  SRclock <- getSimilarity(SRclock, c(4, 5, 2, 3), edgeList)

  newList <- SRclock[0, ]
  nodesInBlocks <- unique(unlist(strsplit(SRclock$Nodes, split = "; ")))
  # take first element and find random block, which it belongs to
  newList <- rbind(newList, SRclock[sample(grep(paste("(^| )", nodesInBlocks[1], "($|;)", sep = ""), SRclock$Nodes), 1), ])
  block <- newList$Nodes[nrow(newList)]
  block <- unlist(strsplit(block, split = "; "))
  
  
  # write to output files
  # write(FFFnoFiberAlist, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_FFFa.txt", sep = ""))
  # write(FFFnoFiberBlist, paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_FFFb.txt", sep = ""))
  # write(unlist(unite(SRclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
  #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_FFFc.txt", sep = ""))
  # write(unlist(unite(SRclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
  #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_FFFd.txt", sep = ""))
  
  # JK section
  JKnoClock <- matrix(
    data = c(0, 0, 0, 1,
             0, 0, 1, 0,
             1, 0, 0, 1,
             0, 1, 1, 0), ncol = 4)
  JKnoClock <- getCircuits(JKnoClock, graph, "NoClock")
  # remove circuits, which are rearrangement of one another
  duplicationTable <- sapply(JKnoClock$Nodes, function(x) arrangeLine(x))
  JKnoClock <- JKnoClock[!duplicated(duplicationTable), ]
  
  # calculate measure of similarity in between Xs and Ys
  JKnoClock <- getSimilarity(JKnoClock, c(3, 4, 1, 2), edgeList)
  
  JKclock <- matrix(
    data = c(0, 0, 0, 0, 0,
             1, 0, 0, 0, 1,
             1, 0, 0, 1, 0,
             0, 1, 0, 0, 1,
             0, 0, 1, 1, 0), ncol = 5)
  JKclock <- getCircuits(JKclock, graph, "Clock")
  # remove circuits, which are rearrangement of one another
  duplicationTable <- sapply(JKclock$Nodes, function(x) arrangeLine(x))
  JKclock <- JKclock[!duplicated(duplicationTable), ]
  
  # calculate measure of similarity in between Xs and Ys
  JKclock <- getSimilarity(JKclock, c(4, 5, 2, 3), edgeList)
  
  # write to output files
  # write(unlist(unite(JKclock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
  #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_JKa.txt", sep = ""))
  # write(unlist(unite(JKnoClock[, c("Nodes", "SimilarityY", "SimilarityX")], Output, sep = "\t")),
  #       paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/", dataset, "_JKb.txt", sep = ""))
  
  # makePngs(dataset, "ARa", possibleARaList)
  # makePngs(dataset, "ARb", possibleARbList)
  # makePngs(dataset, "ARc", possibleARcList)
  # 
  # makePngs(dataset, "FFFa", FFFnoFiberAlist)
  # makePngs(dataset, "FFFb", FFFnoFiberBlist)
  # makePngs(dataset, "FFFc", SRclock$Nodes)
  # makePngs(dataset, "FFFd", SRnoClock$Nodes)
  # 
  # makePngs(dataset, "JKa", JKclock$Nodes)
  # makePngs(dataset, "JKb", JKnoClock$Nodes)
  c(dataset, possibleARa, possibleARb, possibleARc, FFFnoFiberA, FFFnoFiberB, nrow(SRclock), nrow(SRnoClock), nrow(JKclock), nrow(JKnoClock))
}
#stopCluster(myCluster)

summary1 <- as.data.frame(summary, stringsAsFactors = F)
colnames(summary1) <- c("Dataset", "AR A", "AR B", "AR C", "FFF A", "FFF B", "FFF C", "FFF D", "JK A", "JK B")
summary1[, -1] <- sapply(summary1[, -1], as.numeric)

write.csv(summary1, "/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/pvalueSummary.csv", row.names = F, quote = F)

pdf("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/output.pdf", height = 11, width = 14)
grid.table(summary1)
dev.off()


# here we do pdfs

circuitDirectories <- list.dirs(paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/circuits/"), full.names = T)
circuitDirectories <- circuitDirectories[-1]
circuitDirectories <- as.data.frame(circuitDirectories, stringsAsFactors = F)

circuitDirectories$Dataset <- gsub(".*//([A-z]*)_[[:alpha:]]*$", "\\1", circuitDirectories$circuitDirectories)

foreach(dataset = unique(circuitDirectories$Dataset), .combine = c) %do% {
  directoriesToCheck <- circuitDirectories %>%
    filter(Dataset == dataset)
  
  pngFileNames <- foreach(i = 1:nrow(directoriesToCheck), .combine = c) %do% {
    paste(directoriesToCheck$circuitDirectories[i], list.files(directoriesToCheck$circuitDirectories[i]), sep = "/")
  }
  
  pngFileNames
}

setwd("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/")
texFileName <- "structures.tex"

write("\\documentclass[preprint,aps,preprintnumbers,amsmath,amssymb]{revtex4}

\\usepackage[]{graphicx}

\\begin{document}", texFileName)

for(i in 1:length(pngFileNames)) {
  write(paste("\\clearpage ",
              gsub(".*_([[:alpha:]]*)/[^/]*$", "\\1", pngFileNames[i]),
              " \\\\", sep = ""),
        file = texFileName, append = T)
  write(paste("\\includegraphics[width=\\textwidth, height=100in, keepaspectratio]{",
              gsub(".*/(circuits)", "\\1", pngFileNames[i]),
              "} \\\\", sep = ""),
        file = texFileName, append = T)
}

write("\\end{document}", texFileName, append = T)

system(paste("pdflatex", texFileName))

system(paste("mv ",
             gsub(".tex", ".pdf", texFileName),
             " ",
             gsub("[^/]*.tex", paste(dataset, ".pdf", sep = ""), texFileName),
             sep = ""))

# 
# # here we clean 4 networks
# # 1
# trrust_mouse <- read.table("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/trrust_rawdata.mouse.tsv", stringsAsFactors = F)
# trrust_mouse <- trrust_mouse %>%
#   filter(V3 != "Unknown")
# 
# write.table(trrust_mouse[, 1:3], "/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/trust_mouse.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# 
# # 2
# trrust_human <- read.table("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/trrust_rawdata.human.tsv", stringsAsFactors = F)
# trrust_human <- trrust_human %>%
#   filter(V3 != "Unknown")
# 
# write.table(trrust_human[, 1:3], "/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/trust_human.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# 
# # 3
# regnetwork_mouse <- read.table("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/new_kegg.mouse.reg.direction.txt", stringsAsFactors = F)
# unique(regnetwork_mouse$V5)
# regnetwork_mouse$V5 <- gsub("\\-\\->", "Activation", regnetwork_mouse$V5)
# regnetwork_mouse$V5 <- gsub("\\-\\-\\|", "Repression", regnetwork_mouse$V5)
# regnetwork_mouse <- regnetwork_mouse[grepl("(Activation|Repression)", regnetwork_mouse$V5), ]
# 
# write.table(regnetwork_mouse[, c(1, 3, 5)], "/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/regnetwork_mouse.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# 
# # 4
# regnetwork_human <- read.table("/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/new_kegg.human.reg.direction.txt", stringsAsFactors = F)
# unique(regnetwork_human$V5)
# regnetwork_human$V5 <- gsub("\\-\\->", "Activation", regnetwork_human$V5)
# regnetwork_human$V5 <- gsub("\\-\\-\\|", "Repression", regnetwork_human$V5)
# regnetwork_human <- regnetwork_human[grepl("(Activation|Repression)", regnetwork_human$V5), ]
# 
# write.table(regnetwork_human[, c(1, 3, 5)], "/home/ian/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/NetworksAnalysis/cleanedData/regnetwork_human.txt", sep = "\t", quote = F, row.names = F, col.names = F)
