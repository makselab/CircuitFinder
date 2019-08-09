library(tidyr)
library(dplyr)

getDatasetByPrefix <- function(prefix) {
  prefix <- gsub("Regulations_in_ATRM",            "Arabidopsis Thaliana", prefix)
  prefix <- gsub("regulations_bacilus_cleaned",    "Bacillus subtilis", prefix)
  prefix <- gsub("Micobacterium_tuberculosis",     "Micobacterium tuberculosis", prefix)
  prefix <- gsub("salmonella_Typhi_SL1344_RN",     "Salmonella SL1344", prefix)
  prefix <- gsub("SGD_full_yeast_TF-G",            "Yeast SGD", prefix)
  prefix <- gsub("YTRP_TF_gene_regulatory",        "Yeast YTRP regulatory", prefix)
  prefix <- gsub("YTRP_TF_gene_binding",           "Yeast YTRP binding", prefix)
  prefix <- gsub("trrust_rawdata.mouse",           "Mouse TRRUST", prefix)
  prefix <- gsub("trrust_rawdata.human",           "Human TRRUST2", prefix)
  prefix <- gsub("human_network_tf_gene_weighted", "Human TRRUST", prefix)
  prefix <- gsub("combinedKEGG",                   "Human KEGG", prefix)
  return(prefix)
}

getTypeName <- function(type) {
  type <- gsub("AR",        "SR flip-flop", type)
  type <- gsub("FFF",       "Clocked SR flip-flop", type)
  type <- gsub("Fibonacci", "Clocked JK flip-flop", type)
  return(type)
}

structuresPath <- "~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/"
setwd(structuresPath)
circuitDirectories <- list.dirs(paste("~/Dropbox/Research/PhD work/shared folders/MARIANO-SERIES/TRANSISTOR/IAN/circuits3/"), full.names = T)
circuitDirectories <- circuitDirectories[-1]
circuitDirectories <- as.data.frame(circuitDirectories, stringsAsFactors = F)

circuitDirectories$Dataset <- gsub(".*//(.*)_[[:alpha:]]*$", "\\1", circuitDirectories$circuitDirectories)
#circuitDirectories <- circuitDirectories[!grepl("human_network_tf_gene_weighted", circuitDirectories$circuitDirectories), ]
#circuitDirectories <- circuitDirectories[foreach(i = 1:length(fileNames), .combine = c) %do% {grep(paste(fileNames[i], sep = ""), circuitDirectories$Dataset)}, ]

samplePNGs = T
# there is roughly 4 times more asymm circuits than symm, so we sample to get equal amount of both types in structures

proportionToSample = 10/68.4 / 3
asymmProportionToSample = (proportionToSample / 5) ^ (1/2)
symmProportionToSample = asymmProportionToSample * 4

texFileName <- "structures.tex"

write("\\documentclass[preprint,aps,preprintnumbers,amsmath,amssymb]{revtex4}
      
      \\usepackage[]{graphicx}
      
      \\begin{document}", texFileName)

# first we get all asymmetric circuits

foreach(dataset = unique(circuitDirectories$Dataset), .combine = c) %do% {
  directoriesToCheck <- circuitDirectories %>%
    filter(Dataset == dataset)
  
  AR_circuits <- try(read.table(paste("circuits3/", dataset, "_AR.txt", sep = ""), sep = "\t", stringsAsFactors = F))
  FFF_circuits <- try(read.table(paste("circuits3/", dataset, "_FFF.txt", sep = ""), sep = "\t", stringsAsFactors = F))
  Fibonacci_circuits <- try(read.table(paste("circuits3/", dataset, "_Fibonacci.txt", sep = ""), sep = "\t", stringsAsFactors = F))
  
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
  
  if(length(pngFileNames) == 0) {return()}
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
                #paste(type, "Broken"),
                getTypeName(type),
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