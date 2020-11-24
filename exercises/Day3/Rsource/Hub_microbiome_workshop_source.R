# 6/3/2020: source-script for setting up environment to run the USF Omics Hub Microbiome Data-Analysis Workshop pipeline(s).
# Author: J. Oberstaller

# install necessary general packages/dependencies if they aren't already installed (combining AS, JG, JO, and TEK package-requirements)
packages <- c("BiocManager",
              "knitr",
              "caTools",
              "ggpubr",
              "dplyr",
              "RColorBrewer",
              "reshape2",
              "ggplot2",
              "hrbrthemes",
              "gcookbook",
              "RCurl",
              "fs",
              "readr",
              "GUniFrac",
              "simEd",
              "tictoc",
              "vegan"
)
if (length(setdiff(packages,
                   rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages,
                           rownames(installed.packages())),
                   repos = 'https://cloud.r-project.org',
                   type = "binary",
                   update = TRUE,
                   #ask = FALSE,
                   verbose = TRUE
  ) 
}

# install our primary analysis-packages from bioconductor if they aren't already installed
bioconductor.packages <- c("microbiome",
                           "phyloseq",
                           "SIAMCAT",
                           "metagenomeSeq",
                           "dada2",
                           "DESeq2")
if (length(setdiff(bioconductor.packages, rownames(installed.packages()))) > 0){
  BiocManager::install(c("microbiome",
                         "phyloseq",
                         "SIAMCAT",
                         "metagenomeSeq",
                         "dada2"),
                       type = "binary",
                       update = FALSE,
                       #ask = FALSE,
                       verbose = TRUE
  )
}

# now install this one separately as there is no binary available for this version of R yet (GenomeInfoDb is a dependency of ALDEx2)
bioconductor.packages2 <- c("GenomeInfoDb","ALDEx2")

if (length(setdiff(bioconductor.packages2, rownames(installed.packages()))) > 0){
  BiocManager::install(c("GenomeInfoDb", "ALDEx2"),
                       type = "both",
                       update = FALSE,
                       #ask = FALSE,
                       verbose = TRUE
  )
}

# clean up leftover variables from above
rm(list=ls())
#rm("packages", "bioconductor.packages")

# make required directories if they don't already exist
makeRdirs <- function () {
  # make ./Routput directory if it doesn't already exist
  mainDir = getwd()
  subDir = "/Routput"
  
  dir.create(paste(mainDir,
                   subDir,
                   sep = ""),
             showWarnings = FALSE)
  
  # make ./Rfigs directory if it doesn't exist
  subDir2 = "/Rfigs"
  dir.create(paste(mainDir,
                   subDir2,
                   sep = ""),
             showWarnings = FALSE)
}

makeRdirs()

# add these load-package parts to the markdown for the appropriate day of the course instead of this source-document
# load required packages (AS, Day2)
# library(dada2)
# library(RCurl)
# library(fs)
# library(dplyr)
# library(GUniFrac)
# library(simEd)
# library(tictoc)

# load required packages (JG, Day3)
library(microbiome)
library(phyloseq)
library(vegan)
library(metagenomeSeq)
library(knitr)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(hrbrthemes)
library(gcookbook)

# load required packages (TEK, Day3)
library(ALDEx2)
library(readr)
library(dplyr)
library(SIAMCAT)
